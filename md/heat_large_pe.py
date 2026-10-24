"""One-ps restrained low-temperature NVT integrity pilot of the actual large model.
A 10→50→100 K diagnostic, not completed heating to 310 K or equilibration.
"""
import argparse,csv,hashlib,json,math,resource,struct,time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import unit,XmlSerializer,Platform,LangevinMiddleIntegrator,CustomExternalForce,MonteCarloBarostat,NonbondedForce
from openmm.app import AmberPrmtopFile,Simulation,DCDReporter,PDBxFile,element
try:
    from .extend_pe_npt import atomic_json
    from .protein_pilot import aligned_rmsd
    from .minimize_large_pe import geometry
except ImportError:
    from extend_pe_npt import atomic_json
    from protein_pilot import aligned_rmsd
    from minimize_large_pe import geometry


def run(source,bundle,output):
    source,bundle,out=map(lambda p:Path(p).resolve(),(source,bundle,output));out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    atomic_json(out/'progress.json',{'status':'validating_source','executed_steps':0,'planned_steps':1000,'time_ps':0})
    result=json.loads((source/'result.json').read_text());prior_config=json.loads((source/'config.json').read_text());validation=json.loads((source/'validation.json').read_text())
    if not result['allow_heating_preparation'] or validation['status']!='independent_saved_output_checks_passed':raise ValueError('Reviewed source minimization not accepted')
    for name in ['state.xml','system.xml']:
        if hashlib.sha256((source/name).read_bytes()).hexdigest()!=result['sha256'][name]:raise ValueError('Minimization source checksum mismatch')
    prm_path=bundle/'solvated.prmtop'
    if hashlib.sha256(prm_path.read_bytes()).hexdigest()!=prior_config['source_sha256']['solvated.prmtop']:raise ValueError('Topology checksum mismatch')
    prm=AmberPrmtopFile(str(prm_path));system=XmlSerializer.deserialize((source/'system.xml').read_text());start=XmlSerializer.deserialize((source/'state.xml').read_text())
    n=system.getNumParticles()
    if n!=577890 or prm.topology.getNumAtoms()!=n:raise ValueError('Unexpected atom count')
    if any(isinstance(f,MonteCarloBarostat) for f in system.getForces()):raise ValueError('Fixed-box heating must have no barostat')
    pos=next(f for f in system.getForces() if isinstance(f,CustomExternalForce))
    if pos.getNumParticles()!=31974 or any(abs(pos.getParticleParameters(i)[1][0]-1000)>1e-9 for i in range(pos.getNumParticles())):raise ValueError('Position support changed')
    if any(system.getParticleMass(i).value_in_unit(unit.dalton)<=0 for i in range(n)):raise ValueError('Massless particle')
    schedule=[{'target_K':10,'steps':200},{'target_K':50,'steps':200},{'target_K':100,'steps':600}]
    cfg={'atoms':n,'planned_steps':1000,'time_ps':1,'dt_fs':1,'ensemble':'NVT fixed box','schedule':schedule,'friction_per_ps':5,'seed':45,'threads':4,'constraint_tolerance':1e-6,'position_restraint_k_kJ_mol_nm2':1000,'restrained_heavy_atoms':31974,'restraint_anchors':'Unchanged original solvated anchors from reviewed minimization; do not reset anchors','source_minimization':str(source),'source_sha256':result['sha256'],'topology':str(prm_path),'topology_sha256':prior_config['source_sha256']['solvated.prmtop'],'restart_policy':'Keep reviewed minimization coordinates and periodic box; project bond lengths at 1e-6 and initialize new constrained Maxwell velocities at 10 K, seed45. This is a new MD run, not a checkpoint continuation.','trajectory_interval_steps':100,'observations_interval_steps':25,'checkpoint_interval_steps':100,'equilibrium_established':False,'production_ready':False,'limitations':['One-ps low-temperature diagnostic, not heating completed to 310 K','Fixed initially underdense box; no pressure or density adjustment','Observed six-chain core omits unresolved tails and glycans','Finite C24H50 PE patch, not validated bulk crystal','No binding/free energy, spontaneous adsorption or PTFE/PVDF comparison']}
    atomic_json(out/'config.json',cfg);(out/'system.xml').write_text(XmlSerializer.serialize(system))
    integrator=LangevinMiddleIntegrator(10*unit.kelvin,5/unit.picosecond,1*unit.femtoseconds);integrator.setConstraintTolerance(1e-6);integrator.setRandomNumberSeed(45)
    (out/'integrator-initial.xml').write_text(XmlSerializer.serialize(integrator))
    sim=Simulation(prm.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'})
    sim.context.setPeriodicBoxVectors(*start.getPeriodicBoxVectors());sim.context.setPositions(start.getPositions());sim.context.applyConstraints(1e-6);sim.context.setVelocitiesToTemperature(10*unit.kelvin,45);sim.context.applyVelocityConstraints(1e-6)
    initial=sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);reference=initial.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    original=start.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    box_vectors=np.array([v.value_in_unit(unit.nanometer) for v in start.getPeriodicBoxVectors()]);box=np.diag(box_vectors);volume=float(np.linalg.det(box_vectors))
    if not np.allclose(box_vectors,np.diag(box)):raise ValueError('Metrics require orthorhombic box')
    atomlist=list(prm.topology.atoms());maps=json.loads((Path(__file__).resolve().parents[1]/'data/fibrinogen-core03/manifest.json').read_text())['residue_mapping']
    native=lambda a:a.residue.index<len(maps) and maps[a.residue.index]['source_label_seq'] is not None
    ph=np.array([a.index for a in atomlist if native(a) and a.element!=element.hydrogen]);ca=np.array([a.index for a in atomlist if native(a) and a.name=='CA']);pe=np.array([a.index for a in atomlist if a.residue.name=='PEO' and a.element==element.carbon])
    if len(ca)!=1913 or len(pe)!=16560:raise ValueError('Wrong core or polymer identities')
    protein_rest=np.array([int(pos.getParticleParameters(i)[0]) for i in range(pos.getNumParticles()) if atomlist[int(pos.getParticleParameters(i)[0])].residue.name!='PEO']);pe_rest=np.array([int(pos.getParticleParameters(i)[0]) for i in range(pos.getNumParticles()) if atomlist[int(pos.getParticleParameters(i)[0])].residue.name=='PEO'])
    constraints=[system.getConstraintParameters(i) for i in range(system.getNumConstraints())];ia=np.array([int(x[0]) for x in constraints]);ib=np.array([int(x[1]) for x in constraints]);targets=np.array([x[2].value_in_unit(unit.nanometer) for x in constraints]);del constraints
    dof=3*n-system.getNumConstraints()-3;mass=sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(n));density=mass*1.66053906660e-3/volume
    sim.reporters.append(DCDReporter(str(out/'trajectory.dcd'),100))
    fields=['step','time_ps','target_K','temperature_K','potential_including_position_restraints_kJ_mol','kinetic_kJ_mol','total_energy_kJ_mol','total_density_g_cm3','volume_nm3','max_relative_constraint_error','protein_CA_RMSD_from_heating_start_A','max_protein_heavy_shift_from_heating_start_nm','max_PE_carbon_shift_from_heating_start_nm','minimum_native_protein_PE_distance_nm','contact_native_residue_count']
    observations=[];contacts=[]
    def observe(writer,file,target):
        nonlocal contacts
        state=sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);p=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        e=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);ke=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole);temp=2*ke/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
        if not np.isfinite(p).all() or not all(math.isfinite(x) for x in (e,ke,temp)) or temp>350 or temp<0:raise ValueError('Nonfinite or runaway temperature/energy')
        cell=np.array([v.value_in_unit(unit.nanometer) for v in state.getPeriodicBoxVectors()])
        if not np.allclose(cell,box_vectors,atol=1e-10):raise ValueError('Fixed periodic cell changed')
        delta=p[ia]-p[ib];delta-=box*np.round(delta/box);error=float(np.max(np.abs(np.linalg.norm(delta,axis=1)-targets)/targets))
        if error>1e-4:raise ValueError('Excessive bond constraint error')
        displacement=p-reference;displacement-=box*np.round(displacement/box)
        p_shift=float(np.linalg.norm(displacement[protein_rest],axis=1).max());pe_shift=float(np.linalg.norm(displacement[pe_rest],axis=1).max())
        if max(p_shift,pe_shift)>.3:raise ValueError('Excessive restrained heavy-atom movement in short low-T test')
        d=cKDTree(np.mod(p[pe],box),boxsize=box).query(np.mod(p[ph],box))[0]
        if d.min()<.15:raise ValueError('Severe native protein/PE overlap')
        labels={(maps[atomlist[int(ph[i])].residue.index]['source_chain'],maps[atomlist[int(ph[i])].residue.index]['source_label_seq']) for i in np.flatnonzero(d<.45)}
        contacts=[{'source_chain':c,'source_label_seq':int(r)} for c,r in sorted(labels)]
        row=[sim.currentStep,sim.currentStep*.001,target,temp,e,ke,e+ke,density,volume,error,aligned_rmsd(p[ca],reference[ca])*10,p_shift,pe_shift,float(d.min()),len(contacts)]
        writer.writerow(row);file.flush();observations.append(dict(zip(fields,row)))
        atomic_json(out/'progress.json',{'status':'running','executed_steps':sim.currentStep,'planned_steps':1000,'time_ps':sim.currentStep*.001,'elapsed_seconds':time.monotonic()-started,'last':observations[-1],'contact_residues':contacts,'equilibrium_established':False})
        print(f'Low-T test {sim.currentStep}/1000, target {target} K, T {temp:.2f} K, constraint {error:.2g}, elapsed {time.monotonic()-started:.1f}s',flush=True)
        return state
    with (out/'observables.csv').open('w',newline='') as file:
        writer=csv.writer(file);writer.writerow(fields);last=observe(writer,file,10)
        (out/'initial-state.xml').write_text(XmlSerializer.serialize(initial))
        for spec in schedule:
            integrator.setTemperature(spec['target_K']*unit.kelvin)
            for j in range(spec['steps']//25):
                sim.step(25);last=observe(writer,file,spec['target_K'])
                if sim.currentStep%100==0:
                    sim.saveState(str(out/'state-next.xml'));(out/'state-next.xml').replace(out/'state.xml');sim.saveCheckpoint(str(out/'checkpoint-next.chk'));(out/'checkpoint-next.chk').replace(out/'checkpoint.chk')
    # Full final spatial checks; exception pairs are explicitly excluded from hard-overlap counts.
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce));exclusions=set()
    for i in range(nb.getNumExceptions()):
        a,b,*_=nb.getExceptionParameters(i);a,b=sorted((int(a),int(b)));exclusions.add(a*n+b)
    heavy=np.array([a.index for a in atomlist if a.element!=element.hydrogen]);water_o=np.array([a.index for a in atomlist if a.residue.name in ('HOH','WAT') and a.element==element.oxygen]);dry_n=81275
    geom=geometry(last.getPositions(asNumpy=True).value_in_unit(unit.nanometer),box,heavy,water_o,heavy[heavy<dry_n],exclusions)
    if geom['nonexcluded_heavy_pairs_under_0_12nm'] or geom['nonexcluded_allatom_pairs_under_0_08nm']:raise ValueError('Final hard-overlap check failed')
    with (out/'trajectory.dcd').open('rb') as f:header=f.read(12)
    count=struct.unpack('<i',header[8:12])[0]
    if header[4:8]!=b'CORD' or count!=10:raise ValueError('Incomplete DCD')
    sim.saveState(str(out/'state.xml'));sim.saveCheckpoint(str(out/'checkpoint.chk'));(out/'integrator-final.xml').write_text(XmlSerializer.serialize(integrator))
    prm.topology.setPeriodicBoxVectors(last.getPeriodicBoxVectors())
    with (out/'final.cif').open('w') as file:PDBxFile.writeFile(prm.topology,last.getPositions(),file)
    r={'status':'low_temperature_MD_completed_numerical_checks_passed','atoms':n,'steps':1000,'time_ps':1,'dt_fs':1,'initial':observations[0],'final':observations[-1],'maximum_recorded_temperature_K':max(x['temperature_K'] for x in observations),'maximum_recorded_constraint_relative_error':max(x['max_relative_constraint_error'] for x in observations),'geometry_after':geom,'DCD_frames':count,'final_contact_residues':contacts,'initial_constraint_projection_max_atom_shift_nm':float(np.linalg.norm(reference-original,axis=1).max()),'elapsed_seconds':time.monotonic()-started,'peak_process_RSS_GiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'equilibrium_established':False,'heating_to_310K_completed':False,'production_ready':False,'limitations':cfg['limitations'],'sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['system.xml','state.xml','final.cif','observables.csv']}}
    atomic_json(out/'result.json',r);atomic_json(out/'progress.json',{'status':'completed','executed_steps':1000,'planned_steps':1000,'time_ps':1,'elapsed_seconds':time.monotonic()-started,'last':observations[-1],'equilibrium_established':False});print(json.dumps({k:v for k,v in r.items() if k not in ('limitations','sha256')},indent=2),flush=True);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source_minimization');p.add_argument('bundle');p.add_argument('output');a=p.parse_args()
    try:run(a.source_minimization,a.bundle,a.output)
    except BaseException as e:
        out=Path(a.output)
        if out.exists():atomic_json(out/'failure.json',{'status':'failed','error':f'{type(e).__name__}: {e}'})
        raise
