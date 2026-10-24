"""Extend a reviewed small finite PE patch in NPT, retaining coordinates and velocities.
Ensemble changes restart stochastic streams; this is not a bitwise checkpoint continuation.
"""
import argparse,csv,hashlib,json,math,os,struct,time
from pathlib import Path
import numpy as np
from openmm import unit,XmlSerializer,Platform,NonbondedForce,CustomExternalForce,MonteCarloBarostat,LangevinMiddleIntegrator
from openmm.app import AmberPrmtopFile,Simulation,DCDReporter,PDBFile,element
try:
    from .run_pe_complex import surface_metrics
    from .protein_pilot import aligned_rmsd
    from .equilibrate_small_pe import density_g_cm3
except ImportError:
    from run_pe_complex import surface_metrics
    from protein_pilot import aligned_rmsd
    from equilibrate_small_pe import density_g_cm3


def atomic_json(path,data):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,ensure_ascii=False,indent=2));temp.replace(path)


def run(bundle,source,output,ps=100,seed=43,threads=4):
    bundle,source,out=map(lambda x:Path(x).resolve(),(bundle,source,output))
    if ps<=0 or abs(ps/.002-round(ps/.002))>1e-6:raise ValueError('Duration must be positive integer multiples of 2 fs')
    manifest=json.loads((bundle/'manifest.json').read_text())
    for name in ['complex.prmtop','complex.inpcrd']:
        if hashlib.sha256((bundle/name).read_bytes()).hexdigest()!=manifest['sha256'][name]:raise ValueError('Source bundle checksum mismatch')
    # Validate that the reviewed staged run and the model bundle share identical parameters.
    parent_config=json.loads((source.parent/'config.json').read_text())
    if parent_config['input_sha256']['complex.prmtop']!=manifest['sha256']['complex.prmtop']:raise ValueError('Starting state model identity mismatch')
    checks=json.loads((source.parent/'validation.json').read_text())
    expected=next(x['state_sha256'] for x in checks['stages'] if x['stage']==source.name)
    if hashlib.sha256((source/'state.xml').read_bytes()).hexdigest()!=expected:raise ValueError('Starting state differs from validated state')
    prm=AmberPrmtopFile(str(bundle/'complex.prmtop'))
    state=XmlSerializer.deserialize((source/'state.xml').read_text())
    system=XmlSerializer.deserialize((source/'system.xml').read_text())
    if system.getNumParticles()!=43288 or prm.topology.getNumAtoms()!=43288:raise ValueError('Unexpected model')
    for i in reversed(range(system.getNumForces())):
        if isinstance(system.getForce(i),(CustomExternalForce,MonteCarloBarostat)):system.removeForce(i)
    barostat=MonteCarloBarostat(1*unit.bar,310*unit.kelvin,25);barostat.setRandomNumberSeed(seed);system.addForce(barostat)
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
    charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
    if abs(charge)>1e-4:raise ValueError('Non-neutral source system')
    out.mkdir(parents=True,exist_ok=False)
    steps=int(round(ps/.002));start_wall=time.monotonic()
    cfg={'atoms':43288,'planned_steps':steps,'planned_time_ps':ps,'dt_fs':2,'temperature_K':310,'pressure_bar':1,'threads':threads,'seed':seed,'barostat_frequency':25,'position_restraints':0,'input_bundle':str(bundle),'input_sha256':manifest['sha256'],'source_stage':str(source),'starting_state_sha256':expected,'source_state_time_ps':state.getTime().value_in_unit(unit.picosecond),'restart_policy':'Retain source coordinates, box and velocities; change NVT to isotropic NPT, remove position restraints, restart thermostat/barostat random streams. No reminimization or velocity reinitialization. Extension time starts at zero.','limitations':parent_config['limitations']+['Isotropic NPT is a finite-patch preparation protocol, not a validated macroscopic slab protocol.']}
    atomic_json(out/'config.json',cfg);(out/'system.xml').write_text(XmlSerializer.serialize(system))
    integ=LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,2*unit.femtoseconds);integ.setRandomNumberSeed(seed)
    (out/'integrator.xml').write_text(XmlSerializer.serialize(integ))
    sim=Simulation(prm.topology,system,integ,Platform.getPlatformByName('CPU'),{'Threads':str(threads)})
    sim.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors());sim.context.setPositions(state.getPositions());sim.context.setVelocities(state.getVelocities());sim.context.setTime(0)
    atoms=list(prm.topology.atoms());ph=np.array([a.index for a in atoms if a.residue.index<252 and a.element!=element.hydrogen]);pe=np.array([a.index for a in atoms if a.residue.name=='PEO' and a.element!=element.hydrogen]);ca=np.array([a.index for a in atoms if a.residue.index<252 and a.name=='CA']);mapping={a.index:a.residue.index+143 for a in atoms if a.residue.index<252}
    ref=state.getPositions(asNumpy=True).value_in_unit(unit.nanometer);mass=sum(prm._prmtop.getMasses());dof=3*system.getNumParticles()-system.getNumConstraints()-3
    sim.reporters.append(DCDReporter(str(out/'trajectory.dcd'),500)) # 1 ps per frame
    fields=['step','extension_time_ps','temperature_K','total_density_g_cm3','volume_nm3','potential_kJ_mol','min_protein_PE_distance_nm','contact_residues','CA_RMSD_from_extension_start_A']
    rows=[];contacts=[]
    def sample(writer,file):
        nonlocal contacts
        s=sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True)
        xyz=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer);vec=np.array([v.value_in_unit(unit.nanometer) for v in s.getPeriodicBoxVectors()]);vol=float(np.linalg.det(vec))
        if not np.allclose(vec,np.diag(np.diag(vec))):raise ValueError('Orthorhombic metrics require orthorhombic cell')
        e=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);temp=2*s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
        if not np.isfinite(xyz).all() or not all(math.isfinite(x) for x in (temp,e,vol)) or temp>1000 or vol<=0 or min(np.diag(vec))<=2.1:raise ValueError('Unstable coordinates, energy, temperature or box')
        dist,contacts=surface_metrics(xyz,ph,pe,mapping,np.diag(vec))
        if dist<.15:raise ValueError('Severe protein/PE overlap')
        row=[sim.currentStep,sim.currentStep*.002,temp,density_g_cm3(mass,vol),vol,e,dist,len(contacts),aligned_rmsd(xyz[ca],ref[ca])*10]
        writer.writerow(row);file.flush();rows.append(row)
        return s
    def save(s,status):
        # Alternating checkpoint files preserve at least one previous checkpoint if interrupted.
        sim.saveCheckpoint(str(out/'checkpoint-next.chk'));(out/'checkpoint-next.chk').replace(out/'checkpoint.chk')
        (out/'state-next.xml').write_text(XmlSerializer.serialize(s));(out/'state-next.xml').replace(out/'state.xml')
        elapsed=time.monotonic()-start_wall
        atomic_json(out/'progress.json',{'status':status,'executed_steps':sim.currentStep,'planned_steps':steps,'extension_time_ps':sim.currentStep*.002,'planned_time_ps':ps,'elapsed_seconds':elapsed,'estimated_remaining_seconds':elapsed*(steps-sim.currentStep)/sim.currentStep if sim.currentStep else None,'last_observation':dict(zip(fields,rows[-1])),'contact_residue_source_labels':contacts,'updated_unix':time.time(),'equilibration_established':False})
    with (out/'observables.csv').open('w',newline='') as file:
        writer=csv.writer(file);writer.writerow(fields)
        try:
            last=sample(writer,file);save(last,'running')
            while sim.currentStep<steps:
                sim.step(min(100,steps-sim.currentStep));last=sample(writer,file)
                if sim.currentStep%500==0 or sim.currentStep==steps:
                    save(last,'running');print(f'NPT extension {sim.currentStep*.002:.1f}/{ps:g} ps, T={rows[-1][2]:.1f} K, density={rows[-1][3]:.4f}, elapsed={time.monotonic()-start_wall:.1f}s',flush=True)
            prm.topology.setPeriodicBoxVectors(last.getPeriodicBoxVectors())
            with (out/'final.pdb').open('w') as f:PDBFile.writeFile(prm.topology,last.getPositions(),f)
            sxyz=last.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            maxerr=max(abs(np.linalg.norm(sxyz[a]-sxyz[b])-d.value_in_unit(unit.nanometer))/d.value_in_unit(unit.nanometer) for a,b,d in (system.getConstraintParameters(i) for i in range(system.getNumConstraints())))
            with (out/'trajectory.dcd').open('rb') as f:header=f.read(12)
            frames=struct.unpack('<i',header[8:12])[0]
            if frames!=steps//500 or maxerr>1e-3:raise ValueError('Trajectory frame or final constraint validation failed')
            atomic_json(out/'validation.json',{'status':'execution_integrity_checks_passed','atoms':system.getNumParticles(),'steps':sim.currentStep,'DCD_frames':frames,'barostats':1,'restrained_particles':0,'max_relative_constraint_error':float(maxerr),'net_charge_e':charge,'scientific_equilibrium_established':False})
            result={'status':'completed_NPT_extension_equilibrium_not_established','atoms':43288,'steps':steps,'time_ps':ps,'net_charge_e':charge,'initial':dict(zip(fields,rows[0])),'final':dict(zip(fields,rows[-1])),'final_contact_residue_source_labels':contacts,'elapsed_seconds':time.monotonic()-start_wall,'equilibration_established':False,'production_ready':False,'limitations':cfg['limitations']}
            atomic_json(out/'result.json',result);save(last,'completed');print('NPT extension completed',flush=True)
            return result
        except BaseException as e:
            atomic_json(out/'failure.json',{'status':'failed','executed_steps':sim.currentStep,'error':f'{type(e).__name__}: {e}'})
            if rows:
                try:save(sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True),'failed')
                except Exception:pass
            raise

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bundle');p.add_argument('source_stage');p.add_argument('output');p.add_argument('--ps',type=float,default=100);a=p.parse_args();run(a.bundle,a.source_stage,a.output,a.ps)
