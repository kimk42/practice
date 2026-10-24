"""Run an explicit-solvent protein + finite PE-patch pilot. Not an adsorption-free-energy calculation."""
import argparse,csv,hashlib,json,math,os,time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from openmm import unit,LangevinMiddleIntegrator,Platform,CustomExternalForce,NonbondedForce,XmlSerializer,MinimizationReporter
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,Simulation,PDBFile,DCDReporter,PME,HBonds,element
try:
    from .protein_pilot import aligned_rmsd
except ImportError:
    from protein_pilot import aligned_rmsd


def surface_metrics(xyz,protein_heavy,polymer_heavy,atom_residues,box_lengths):
    # Minimum-image pair distances for an orthorhombic periodic cell.
    tree=cKDTree(np.mod(xyz[polymer_heavy],box_lengths),boxsize=box_lengths)
    distance,_=tree.query(np.mod(xyz[protein_heavy],box_lengths),k=1)
    contact_residues=sorted({atom_residues[protein_heavy[i]] for i,d in enumerate(distance) if d<.45})
    return float(distance.min()),contact_residues


def run(bundle,output,steps=1000):
    root=Path(bundle);out=Path(output);out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((root/'manifest.json').read_text())
    for name in ('complex.prmtop','complex.inpcrd'):
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=manifest['sha256'][name]:raise ValueError('Model checksum mismatch')
    prm=AmberPrmtopFile(str(root/'complex.prmtop'));crd=AmberInpcrdFile(str(root/'complex.inpcrd'))
    xyz=crd.positions.value_in_unit(unit.nanometer)
    system=prm.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*unit.nanometer,constraints=HBonds,rigidWater=True,ewaldErrorTolerance=.0005)
    system.setDefaultPeriodicBoxVectors(*crd.boxVectors)
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
    charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
    if abs(charge)>1e-4:raise ValueError(f'Combined solvated system not neutral: {charge}')
    atoms=list(prm.topology.atoms());residues=list(prm.topology.residues())
    protein_heavy=np.array([a.index for a in atoms if a.residue.index<252 and a.element!=element.hydrogen])
    polymer_heavy=np.array([a.index for a in atoms if a.residue.name=='PEO' and a.element!=element.hydrogen])
    ca=np.array([a.index for a in atoms if a.residue.index<252 and a.name=='CA'])
    if len(ca)!=252 or len(polymer_heavy)!=24*manifest['polymer_chains']:raise ValueError('Unexpected protein/polymer identity counts')
    atom_residues={a.index:int(a.residue.index)+143 for a in atoms if a.residue.index<252}
    vectors=np.array([v.value_in_unit(unit.nanometer) for v in crd.boxVectors])
    if not np.allclose(vectors,np.diag(np.diag(vectors))):raise ValueError('Metrics require orthorhombic box')
    box=np.diag(vectors)
    # The first half of the oligomers form the lower layer; restrain only their carbons.
    pe_res=[r for r in residues if r.name=='PEO'];lower=set(pe_res[:len(pe_res)//2])
    restraint=CustomExternalForce('0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2');restraint.addGlobalParameter('k',1000)
    for name in ('x0','y0','z0'):restraint.addPerParticleParameter(name)
    for a in atoms:
        if a.residue in lower and a.element==element.carbon:restraint.addParticle(a.index,xyz[a.index])
    system.addForce(restraint)
    (out/'system.xml').write_text(XmlSerializer.serialize(system))
    integrator=LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,1*unit.femtoseconds);integrator.setRandomNumberSeed(42)
    sim=Simulation(prm.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'});sim.context.setPositions(crd.positions)
    def get():
        s=sim.context.getState(getPositions=True,getEnergy=True);e=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        p=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        if not math.isfinite(e) or not np.isfinite(p).all():raise ValueError('Non-finite state')
        return s,p,e
    s,original,energy=get();initial_gap,initial_contacts=surface_metrics(original,protein_heavy,polymer_heavy,atom_residues,box)
    with (out/'initial.pdb').open('w') as f:PDBFile.writeFile(prm.topology,s.getPositions(),f)
    print('Prepared mixed system:',system.getNumParticles(),'atoms; net charge',charge,'e; initial min distance',initial_gap,'nm',flush=True)
    class Progress(MinimizationReporter):
        def report(self,iteration,x,grad,args):
            if iteration%10==0:print('Minimizing',iteration,'energy',args['system energy'],flush=True)
            return False
    sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=50,reporter=Progress())
    s,reference,min_energy=get()
    with (out/'minimized.pdb').open('w') as f:PDBFile.writeFile(prm.topology,s.getPositions(),f)
    sim.context.setVelocitiesToTemperature(310*unit.kelvin,42)
    sim.reporters.append(DCDReporter(str(out/'trajectory.dcd'),100))
    started=time.monotonic();rows=[];dof=3*system.getNumParticles()-system.getNumConstraints()-3
    with (out/'observables.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['step','time_ps','potential_kJ_mol','temperature_K','min_protein_PE_heavy_distance_nm','contact_residues_within_0.45nm','protein_CA_RMSD_from_minimized_A'])
        for start in range(0,steps,100):
            sim.step(min(100,steps-start));s,p,e=get();distance,contacts=surface_metrics(p,protein_heavy,polymer_heavy,atom_residues,box)
            temp=2*s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
            row=[sim.currentStep,sim.currentStep*.001,e,temp,distance,len(contacts),aligned_rmsd(p[ca],reference[ca])*10]
            rows.append(row);writer.writerow(row);f.flush()
            print('Step',sim.currentStep,'distance',round(distance,4),'nm; contacts',len(contacts),'T',round(temp,1),'K',flush=True)
    with (out/'final.pdb').open('w') as f:PDBFile.writeFile(prm.topology,s.getPositions(),f)
    sim.saveCheckpoint(str(out/'checkpoint.chk'));sim.saveState(str(out/'state.xml'))
    result={'status':'completed','system':'actual fibrinogen gamma fragment + finite PE patch + explicit TIP3P water and neutralizing ions','model_manifest':manifest,
            'atoms':system.getNumParticles(),'steps':steps,'timestep_fs':1,'temperature_target_K':310,'platform':'CPU','threads':4,'seed':42,'net_charge_e':charge,
            'restraint':'1000 kJ/mol/nm^2 on lower-layer PE carbons; upper layer unrestrained','minimization_max_iterations':50,'minimization_convergence':'not assessed',
            'initial_min_distance_nm':initial_gap,'initial_contact_residues':initial_contacts,'final_min_distance_nm':rows[-1][4],'final_contact_residues':contacts,
            'final_protein_CA_RMSD_A':rows[-1][-1],'final_temperature_K':rows[-1][3],'MD_elapsed_seconds':time.monotonic()-started,
            'limitations':['Short preparation pilot, not equilibrated adsorption','Pre-positioned near surface; contact formation does not establish spontaneous adsorption','Finite PE oligomer model not validated against experimental surface properties','No binding energy or adsorption free energy computed','No PTFE/PVDF comparison has been performed']}
    result['contact_details']=[{'source_label_sequence_position':i,'residue':residues[i-143].name,'distance_in_residues_from_fragment_boundary':min(i-143,394-i)} for i in contacts]
    if any(min(i-143,394-i)<5 for i in contacts):
        result['boundary_warning']='Contact occurs within five residues of an artificial fragment boundary; do not generalize it to intact fibrinogen.'
    (out/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='model_manifest'},indent=2))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bundle');p.add_argument('output');p.add_argument('--steps',type=int,default=1000);a=p.parse_args()
    if not 100<=a.steps<=100000:p.error('steps must be 100..100000')
    run(a.bundle,a.output,a.steps)
