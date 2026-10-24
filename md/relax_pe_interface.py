"""Short post-NPT, fixed-volume interface-preparation diagnostic with explicit support.
Protein restraints are temporarily reintroduced at the NPT endpoint and then released.
Bottom PE carbon support remains; this is not an unrestrained adsorption production run.
"""
import argparse,csv,json,math,struct,time
from pathlib import Path
import numpy as np
from openmm import unit,XmlSerializer,Platform,CustomExternalForce,MonteCarloBarostat,LangevinMiddleIntegrator
from openmm.app import AmberPrmtopFile,Simulation,DCDReporter,PDBFile,element
try:
    from .extend_pe_npt import atomic_json
    from .run_pe_complex import surface_metrics
    from .protein_pilot import aligned_rmsd
    from .analyze_pe_followup import structural_snapshot
except ImportError:
    from extend_pe_npt import atomic_json
    from run_pe_complex import surface_metrics
    from protein_pilot import aligned_rmsd
    from analyze_pe_followup import structural_snapshot


def relax(npt,bundle,analysis,output):
    npt,bundle,analysis,out=map(Path,(npt,bundle,analysis,output))
    check=json.loads((analysis/'analysis.json').read_text())
    if not check['allow_next_preparation']:raise ValueError('Prior structural checks require review')
    prior=json.loads((npt/'result.json').read_text())
    if prior['steps']!=50000:raise ValueError('100 ps NPT stage incomplete')
    out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    prm=AmberPrmtopFile(str(bundle/'complex.prmtop'));system=XmlSerializer.deserialize((npt/'system.xml').read_text());state=XmlSerializer.deserialize((npt/'state.xml').read_text())
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i),(MonteCarloBarostat,CustomExternalForce)):system.removeForce(i)
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    atoms=list(prm.topology.atoms());pe_res=[r for r in prm.topology.residues() if r.name=='PEO'];lower=set(pe_res[:24]);xyz=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    force=CustomExternalForce('0.5*(is_protein*protein_k+is_lower*lower_k)*periodicdistance(x,y,z,x0,y0,z0)^2');force.addGlobalParameter('protein_k',100);force.addGlobalParameter('lower_k',1000)
    for name in ['is_protein','is_lower','x0','y0','z0']:force.addPerParticleParameter(name)
    for a in atoms:
        if a.element==element.hydrogen:continue
        prot=a.residue.index<252;bottom=a.residue in lower
        if prot or bottom:force.addParticle(a.index,[int(prot),int(bottom),*xyz[a.index]])
    system.addForce(force)
    integrator=LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,2*unit.femtoseconds);integrator.setRandomNumberSeed(44)
    sim=Simulation(prm.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'})
    sim.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors());sim.context.setPositions(state.getPositions());sim.context.setVelocities(state.getVelocities());sim.context.setTime(0)
    (out/'system-initial.xml').write_text(XmlSerializer.serialize(system));(out/'system.xml').write_text(XmlSerializer.serialize(system));(out/'integrator.xml').write_text(XmlSerializer.serialize(integrator))
    schedule=[{'protein_k':100,'steps':1000,'time_ps':2},{'protein_k':10,'steps':1000,'time_ps':2},{'protein_k':0,'steps':3000,'time_ps':6}]
    config={'ensemble':'NVT','atoms':43288,'dt_fs':2,'temperature_K':310,'lower_PE_carbon_k_kJ_mol_nm2':1000,'upper_PE_k':0,'schedule':schedule,'seed':44,'source':str(npt),'policy':'Temporary protein restraints anchored at the post-NPT coordinates, then released; lower PE support retained throughout. Not a continuation of progressively releasing the already unrestrained NPT stage.','equilibrium_established':False}
    atomic_json(out/'config.json',config)
    ph=np.array([a.index for a in atoms if a.residue.index<252 and a.element!=element.hydrogen]);pe=np.array([a.index for a in atoms if a.residue.name=='PEO' and a.element==element.carbon]);ca=np.array([a.index for a in atoms if a.residue.index<252 and a.name=='CA']);mapping={a.index:a.residue.index+143 for a in atoms if a.residue.index<252};dof=3*system.getNumParticles()-system.getNumConstraints()-3
    sim.reporters.append(DCDReporter(str(out/'trajectory.dcd'),500));rows=[]
    fields=['step','interface_time_ps','protein_k_kJ_mol_nm2','temperature_K','potential_kJ_mol','min_protein_PE_distance_nm','contact_residues','CA_RMSD_from_NPT_endpoint_A']
    try:
        with (out/'observables.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(fields)
            for i,stage in enumerate(schedule):
                sim.context.setParameter('protein_k',stage['protein_k'])
                for j in range(stage['steps']//100):
                    sim.step(100);s=sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);p=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer);box=np.array([v.value_in_unit(unit.nanometer) for v in s.getPeriodicBoxVectors()]);e=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);t=2*s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin));dist,contacts=surface_metrics(p,ph,pe,mapping,np.diag(box))
                    if not np.isfinite(p).all() or not all(math.isfinite(x) for x in [e,t,dist]) or t>1000 or dist<.15:raise ValueError('Unstable interface state')
                    row=[sim.currentStep,sim.currentStep*.002,stage['protein_k'],t,e,dist,len(contacts),aligned_rmsd(p[ca],xyz[ca])*10];rows.append(row);w.writerow(row);f.flush()
                    if sim.currentStep%500==0:
                        sim.saveState(str(out/'state-next.xml'));(out/'state-next.xml').replace(out/'state.xml');sim.saveCheckpoint(str(out/'checkpoint.chk'))
                        atomic_json(out/'progress.json',{'status':'running','interface_time_ps':sim.currentStep*.002,'planned_time_ps':10,'protein_k':stage['protein_k'],'elapsed_seconds':time.monotonic()-started,'equilibrium_established':False});print('Interface',sim.currentStep*.002,'ps, protein k',stage['protein_k'],flush=True)
                with (out/f'stage-{i+1}-final.pdb').open('w') as handle:
                    prm.topology.setPeriodicBoxVectors(s.getPeriodicBoxVectors());PDBFile.writeFile(prm.topology,s.getPositions(),handle)
        snap=structural_snapshot(prm,s,system)
        if snap['max_relative_constraint_error']>1e-3:raise ValueError('Constraint validation failed')
        force.setGlobalParameterDefaultValue(0,0)  # Serialized final system agrees with final Context protein_k=0.
        (out/'system.xml').write_text(XmlSerializer.serialize(system))
        sim.saveState(str(out/'state.xml'));sim.saveCheckpoint(str(out/'checkpoint.chk'))
        with (out/'final.pdb').open('w') as handle:PDBFile.writeFile(prm.topology,s.getPositions(),handle)
        with (out/'trajectory.dcd').open('rb') as f:frames=struct.unpack('<i',f.read(12)[8:12])[0]
        if frames!=10:raise ValueError('Unexpected interface DCD frame count')
        result={'status':'interface_preparation_completed_equilibrium_not_established','steps':5000,'time_ps':10,'atoms':43288,'final_protein_k':0,'lower_PE_carbon_k':1000,'final_contact_residue_source_labels':contacts,'final':dict(zip(fields,rows[-1])),'final_structure':snap,'DCD_frames':frames,'elapsed_seconds':time.monotonic()-started,'equilibrium_established':False,'production_ready':False,'limitations':['10 ps fixed-volume supported-patch preparation, not equilibrium adsorption','Protein restraints released; lower PE carbons remain supported','Uncapped isolated gamma fragment; artificial terminal contacts remain a concern','No affinity, free energy, spontaneous adsorption, PTFE or PVDF result']}
        atomic_json(out/'result.json',result);atomic_json(out/'progress.json',{'status':'completed','interface_time_ps':10,'elapsed_seconds':time.monotonic()-started});return result
    except BaseException as e:
        atomic_json(out/'failure.json',{'status':'failed','executed_steps':sim.currentStep,'error':f'{type(e).__name__}: {e}'});raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('npt');p.add_argument('bundle');p.add_argument('analysis');p.add_argument('output');a=p.parse_args();relax(a.npt,a.bundle,a.analysis,a.output)
