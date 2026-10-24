"""Periodic restrained minimization of the actual solvated six-chain core/PE model.
No heating, dynamics, pressure equilibration or adsorption/free-energy claim.
"""
import argparse,csv,hashlib,json,math,os,resource,time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import unit,XmlSerializer,Platform,VerletIntegrator,CustomExternalForce,NonbondedForce,MinimizationReporter,LocalEnergyMinimizer
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,Simulation,PME,HBonds,PDBxFile,element
try:
    from .extend_pe_npt import atomic_json
except ImportError:
    from extend_pe_npt import atomic_json


def geometry(xyz,box,heavy,water_o,solute_heavy,exceptions):
    n=len(xyz);wrapped=np.mod(xyz,box)
    heavy_pairs=cKDTree(wrapped[heavy],boxsize=box).query_pairs(.12,output_type='ndarray')
    severe_heavy=[]
    for a,b in heavy_pairs:
        i,j=sorted((int(heavy[a]),int(heavy[b])))
        if i*n+j not in exceptions:severe_heavy.append([i,j])
    short=cKDTree(wrapped,boxsize=box).query_pairs(.08,output_type='ndarray')
    severe_all=[]
    for a,b in short:
        i,j=sorted((int(a),int(b)))
        if i*n+j not in exceptions:severe_all.append([i,j])
    dist=cKDTree(wrapped[solute_heavy],boxsize=box).query(wrapped[water_o])[0]
    return {'nonexcluded_heavy_pairs_under_0_12nm':len(severe_heavy),'nonexcluded_allatom_pairs_under_0_08nm':len(severe_all),'example_heavy_atom_index_pairs':severe_heavy[:30],'example_allatom_index_pairs':severe_all[:30],'minimum_water_O_solute_heavy_distance_nm':float(dist.min()),'policy':'Minimum-image periodic distances. All NonbondedForce exception pairs (including bonded 1-2/1-3/1-4 pairs) excluded from severe-contact counts; thresholds are hard-overlap diagnostics, not all steric clashes.'}


def run(bundle,output,max_iterations=200,start_state=None):
    bundle,out=map(lambda x:Path(x).resolve(),(bundle,output));out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    atomic_json(out/'progress.json',{'status':'loading_and_checking_input','elapsed_seconds':0})
    manifest=json.loads((bundle/'manifest.json').read_text())
    if manifest['status']!='solvated_and_neutralized_input_not_equilibrated':raise ValueError('Unexpected source state')
    for name in ['solvated.prmtop','solvated.inpcrd']:
        if hashlib.sha256((bundle/name).read_bytes()).hexdigest()!=manifest['sha256'][name]:raise ValueError('Source checksum differs')
    config={'atoms':577890,'task':'restrained energy minimization only','source_bundle':str(bundle),'source_sha256':manifest['sha256'],'protein_k_kJ_mol_nm2':1000,'all_PE_carbon_k_kJ_mol_nm2':1000,'Ca_and_Na':'unrestrained','water':'TIP3P rigid water','constraints':'HBonds; all physical masses retained','electrostatics':'PME, 1 nm cutoff, Ewald error tolerance 0.0005','box_policy':'Fixed periodic box; no barostat','threads':4,'max_iterations_per_constraint_repair_attempt':max_iterations,'target_minimizer_tolerance_kJ_mol_nm':10,'wall_budget_seconds':1800,'MD_steps':0,'equilibrium_established':False,'production_ready':False}
    atomic_json(out/'config.json',config)
    prm=AmberPrmtopFile(str(bundle/'solvated.prmtop'));crd=AmberInpcrdFile(str(bundle/'solvated.inpcrd'))
    xyz=np.asarray(crd.positions.value_in_unit(unit.nanometer));vectors=np.array([v.value_in_unit(unit.nanometer) for v in crd.boxVectors]);box=np.diag(vectors)
    if not np.allclose(vectors,np.diag(box)) or len(xyz)!=577890:raise ValueError('Unsupported model/cell')
    start_xyz=xyz
    if start_state is not None:
        restart=XmlSerializer.deserialize(Path(start_state).read_text())
        start_xyz=restart.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        restart_box=np.array([v.value_in_unit(unit.nanometer) for v in restart.getPeriodicBoxVectors()])
        if len(start_xyz)!=len(xyz) or not np.allclose(restart_box,vectors):raise ValueError('Restart particle count/box mismatch')
        config['restart_state']=str(Path(start_state).resolve())
        config['restart_sha256']=hashlib.sha256(Path(start_state).read_bytes()).hexdigest()
        config['restart_policy']='Use preserved diagnostic coordinates from prior budget-limited attempt; reapply exact bond constraints before minimization. Position anchors remain the original solvated coordinates.'
        atomic_json(out/'config.json',config)
    atoms=list(prm.topology.atoms());dry_n=manifest['dry_atoms'];heavy=np.array([a.index for a in atoms if a.element!=element.hydrogen]);water_o=np.array([a.index for a in atoms if a.residue.name in ('HOH','WAT') and a.element==element.oxygen]);solute_heavy=heavy[heavy<dry_n]
    ss=sum(a.name==b.name=='SG' for a,b in prm.topology.bonds())
    if ss!=25:raise ValueError('Source disulfide count differs')
    print('Building actual periodic system',len(xyz),'atoms',flush=True)
    system=prm.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*unit.nanometer,constraints=HBonds,rigidWater=True,ewaldErrorTolerance=.0005)
    system.setDefaultPeriodicBoxVectors(*crd.boxVectors)
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce));n=system.getNumParticles()
    charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(n))
    if abs(charge)>1e-4:raise ValueError('Non-neutral system')
    exclusions=set()
    for k in range(nb.getNumExceptions()):
        a,b,*_=nb.getExceptionParameters(k);i,j=sorted((int(a),int(b)));exclusions.add(i*n+j)
    print('Checking periodic close contacts',flush=True)
    before=geometry(start_xyz,box,heavy,water_o,solute_heavy,exclusions);atomic_json(out/'geometry-before.json',before)
    # Soft restraints preserve the observed scaffold and finite patch while waters/H/ions relax.
    restraint=CustomExternalForce('0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2')
    for name in ['k','x0','y0','z0']:restraint.addPerParticleParameter(name)
    restrained=[];protein_heavy=[];pe_carbon=[]
    for a in atoms:
        if a.index>=dry_n or a.element==element.hydrogen:continue
        if a.residue.name=='PEO':pe_carbon.append(a.index)
        elif a.element==element.calcium:continue
        else:protein_heavy.append(a.index)
        restraint.addParticle(a.index,[1000,*xyz[a.index]]);restrained.append(a.index)
    restraint.setForceGroup(31);system.addForce(restraint)
    integrator=VerletIntegrator(1*unit.femtosecond) # integration never stepped
    sim=Simulation(prm.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'})
    sim.context.setPeriodicBoxVectors(*crd.boxVectors);sim.context.setPositions(start_xyz*unit.nanometer)
    raw=sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    sim.context.applyConstraints(1e-6)
    initial=sim.context.getState(getPositions=True,getEnergy=True);initial_e=initial.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    (out/'system.xml').write_text(XmlSerializer.serialize(system));(out/'integrator.xml').write_text(XmlSerializer.serialize(integrator))
    (out/'initial-constrained-state.xml').write_text(XmlSerializer.serialize(initial))
    atomic_json(out/'progress.json',{'status':'minimizing','elapsed_seconds':time.monotonic()-started,'raw_input_energy_kJ_mol':raw,'initial_constrained_total_energy_kJ_mol':initial_e,'geometry_before':before,'restrained_protein_heavy':len(protein_heavy),'restrained_PE_carbons':len(pe_carbon)})
    count=0;last_gradient=None
    with (out/'minimization.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['callback','optimizer_iteration','elapsed_seconds','physical_plus_position_restraint_kJ_mol','constraint_penalty_kJ_mol','max_relative_constraint_error','objective_gradient_RMS_kJ_mol_nm'])
        class Progress(MinimizationReporter):
            def report(self,iteration,x,grad,args):
                nonlocal count,last_gradient
                count+=1;last_gradient=float(np.sqrt(np.mean(np.asarray(grad)**2)));elapsed=time.monotonic()-started
                row=[count,iteration,elapsed,args['system energy'],args['restraint energy'],args['max constraint error'],last_gradient]
                if not all(math.isfinite(float(v)) for v in row):raise ValueError('Non-finite minimizer state')
                writer.writerow(row);f.flush()
                if count%10==0 or count==1:
                    atomic_json(out/'progress.json',{'status':'minimizing','callbacks':count,'optimizer_iteration':iteration,'elapsed_seconds':elapsed,'total_energy_kJ_mol':args['system energy'],'constraint_penalty_kJ_mol':args['restraint energy'],'max_relative_constraint_error':args['max constraint error'],'objective_gradient_RMS_kJ_mol_nm':last_gradient})
                    print('Minimization callback',count,'iteration',iteration,'energy',args['system energy'],'elapsed',round(elapsed,1),'s',flush=True)
                if count%50==0:
                    s=sim.context.getState(getPositions=True,getEnergy=True);(out/'inspection-state-next.xml').write_text(XmlSerializer.serialize(s));(out/'inspection-state-next.xml').replace(out/'inspection-state.xml')
                if elapsed>1800:raise TimeoutError('Minimization wall budget exceeded; review partial diagnostic state before proceeding')
                return False
        LocalEnergyMinimizer.minimize(sim.context,10*unit.kilojoule_per_mole/unit.nanometer,max_iterations,Progress())
    final=sim.context.getState(getPositions=True,getEnergy=True,getForces=True)
    p=final.getPositions(asNumpy=True).value_in_unit(unit.nanometer);energy=final.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    if not np.isfinite(p).all() or not math.isfinite(energy):raise ValueError('Non-finite final state')
    after=geometry(p,box,heavy,water_o,solute_heavy,exclusions)
    maxerr=0
    for k in range(system.getNumConstraints()):
        a,b,d=system.getConstraintParameters(k);target=d.value_in_unit(unit.nanometer);delta=p[a]-p[b];delta-=box*np.round(delta/box);maxerr=max(maxerr,abs(np.linalg.norm(delta)-target)/target)
    delta=p-xyz;delta-=box*np.round(delta/box)
    valid=energy<initial_e and maxerr<1e-3 and after['nonexcluded_heavy_pairs_under_0_12nm']==0 and after['nonexcluded_allatom_pairs_under_0_08nm']==0
    (out/'state.xml').write_text(XmlSerializer.serialize(final));sim.saveCheckpoint(str(out/'checkpoint.chk'))
    prm.topology.setPeriodicBoxVectors(crd.boxVectors)
    with (out/'minimized.cif').open('w') as f:PDBxFile.writeFile(prm.topology,final.getPositions(),f)
    flat=(p*10).ravel()
    with (out/'minimized.inpcrd').open('w') as f:
        f.write('Restrained periodic core/PE/TIP3P minimization only; not equilibrated\n');f.write(f'{n:6d}\n')
        for i in range(0,len(flat),6):f.write(''.join(f'{v:12.7f}' for v in flat[i:i+6])+'\n')
        f.write(''.join(f'{v:12.7f}' for v in [*(box*10),90,90,90])+'\n')
    physical=sim.context.getState(getEnergy=True,groups={0}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    result={'status':'restrained_minimization_completed_preliminary_checks_passed' if valid else 'minimization_completed_requires_review','atoms':n,'MD_steps':0,'disulfides':ss,'net_charge_e':charge,'before':before,'after':after,'raw_input_energy_kJ_mol':raw,'initial_constrained_total_energy_kJ_mol':initial_e,'final_total_energy_kJ_mol':energy,'final_forcefield_energy_kJ_mol':physical,'final_position_restraint_energy_kJ_mol':energy-physical,'max_relative_constraint_error':float(maxerr),'max_protein_heavy_displacement_nm':float(np.linalg.norm(delta[protein_heavy],axis=1).max()),'max_PE_carbon_displacement_nm':float(np.linalg.norm(delta[pe_carbon],axis=1).max()),'recorded_minimizer_callbacks':count,'last_reported_objective_gradient_RMS_kJ_mol_nm':last_gradient,'convergence':'Not certified; iteration-limited minimization with constraint repair. Gradient includes internal constraint penalty and is not an independent physical-force convergence test.','elapsed_seconds':time.monotonic()-started,'peak_process_RSS_GiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'box_vectors_nm':vectors.tolist(),'initial_total_density_g_cm3':manifest['initial_total_density_g_cm3'],'equilibrium_established':False,'production_ready':False,'allow_heating_preparation':bool(valid),'limitations':['Observed six-chain core, not full-length fibrinogen; glycans and long unresolved tails omitted','Finite oligomer PE patch, not validated bulk crystal','Soft restraints on protein heavy atoms and all PE carbons; water/H/Ca/Na mobile subject to bond constraints','Fixed box; minimization cannot establish liquid density, temperature, pressure or interface equilibrium','No adsorption, binding/free energy, PTFE or PVDF result'],'sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['system.xml','state.xml','minimized.cif','minimized.inpcrd']}}
    atomic_json(out/'result.json',result);atomic_json(out/'progress.json',{'status':'completed' if valid else 'requires_review','elapsed_seconds':time.monotonic()-started,'final_total_energy_kJ_mol':energy,'allow_heating_preparation':bool(valid)})
    print(json.dumps({k:v for k,v in result.items() if k not in ('before','after','limitations','sha256')},indent=2),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('bundle');p.add_argument('output');p.add_argument('--iterations',type=int,default=200);p.add_argument('--start-state');a=p.parse_args()
    if not 1<=a.iterations<=500:p.error('iterations must be 1..500')
    try:run(a.bundle,a.output,a.iterations,a.start_state)
    except BaseException as e:
        out=Path(a.output)
        if out.exists():atomic_json(out/'failure.json',{'status':'failed','error':f'{type(e).__name__}: {e}'})
        raise
