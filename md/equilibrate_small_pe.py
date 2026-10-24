"""Staged preparation of the 43,288-atom finite PE/protein/water pilot.
Short restrained heating and NPT checks are not proof of interface equilibration.
"""
import argparse,csv,hashlib,json,math,time
from pathlib import Path
import numpy as np
from openmm import unit,XmlSerializer,Platform,LangevinMiddleIntegrator,CustomExternalForce,NonbondedForce,MonteCarloBarostat,MinimizationReporter
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,Simulation,PME,HBonds,PDBFile,DCDReporter,element
try:
    from .run_pe_complex import surface_metrics
    from .protein_pilot import aligned_rmsd
except ImportError:
    from run_pe_complex import surface_metrics
    from protein_pilot import aligned_rmsd


def density_g_cm3(mass_daltons,volume_nm3):return mass_daltons*1.66053906660e-3/volume_nm3

def run(bundle,start_run,output,npt_ps=8):
    root=Path(bundle).resolve();previous=Path(start_run).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((root/'manifest.json').read_text())
    for name in ('complex.prmtop','complex.inpcrd'):
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=manifest['sha256'][name]:raise ValueError('Input checksum mismatch')
    previous_result=json.loads((previous/'result.json').read_text())
    if previous_result['status']!='completed':raise ValueError('Previous run incomplete')
    if previous_result['model_manifest']['sha256']['complex.prmtop']!=manifest['sha256']['complex.prmtop']:raise ValueError('Starting state model identity differs from input bundle')
    prm=AmberPrmtopFile(str(root/'complex.prmtop'));crd=AmberInpcrdFile(str(root/'complex.inpcrd'))
    start=XmlSerializer.deserialize((previous/'state.xml').read_text());positions=start.getPositions();box=start.getPeriodicBoxVectors()
    if len(positions)!=43288 or prm.topology.getNumAtoms()!=43288:raise ValueError('This staged pilot expects the reviewed small model')
    atoms=list(prm.topology.atoms());pe_res=[r for r in prm.topology.residues() if r.name=='PEO'];lower=set(pe_res[:24])
    protein_heavy=np.array([a.index for a in atoms if a.residue.index<252 and a.element!=element.hydrogen]);pe_heavy=np.array([a.index for a in atoms if a.residue.name=='PEO' and a.element!=element.hydrogen]);ca=np.array([a.index for a in atoms if a.residue.index<252 and a.name=='CA'])
    mapping={a.index:a.residue.index+143 for a in atoms if a.residue.index<252}
    reference=np.array(positions.value_in_unit(unit.nanometer));mass=sum(prm._prmtop.getMasses());velocities=None;total_steps=0;started=time.monotonic()
    stages=[{'name':'01-restrained-NVT-heating','ensemble':'NVT','schedule':[(100,250),(200,250),(310,500)],'k_protein':1000,'k_lower_PE':1000,'k_upper_PE':1000},
            {'name':'02-unrestrained-NPT-density','ensemble':'NPT','schedule':[(310,int(round(npt_ps/.002)))],'k_protein':0,'k_lower_PE':0,'k_upper_PE':0},
            {'name':'03-restrained-interface-NVT-check','ensemble':'NVT','schedule':[(310,1000)],'k_protein':100,'k_lower_PE':1000,'k_upper_PE':0}]
    config={'input_bundle':str(root),'input_sha256':manifest['sha256'],'friction_per_ps':1,'atoms':43288,'dt_fs':2,'threads':4,'seed':42,'pressure_bar':1,'npt_barostat_frequency':25,'stages':stages,
            'starting_state_sha256':hashlib.sha256((previous/'state.xml').read_bytes()).hexdigest(),
            'restart_policy':'Start a new staged preparation from prior coordinates, minimize and initialize velocities at 100 K; do not represent this as continuation of the prior ensemble.',
            'NPT_policy':'Isotropic NPT for a finite oligomer patch with all Cartesian position restraints disabled; do not apply this unchanged to an infinite slab or membrane.',
            'minimization_policy':'Restrained, at most 100 iterations per constraint-repair attempt; convergence is not asserted.',
            'limitations':['Protein is an uncapped gamma fragment, not the newly prepared capped six-chain core','Short preparation stages, not established interface equilibrium','Finite PE oligomer model, no bulk crystal validation','Neutralizing ions only; no prescribed 150 mM salt','No binding/adsorption free energy or spontaneous-adsorption claim']}
    (out/'config.json').write_text(json.dumps(config,indent=2))
    rows=[];stage_reports=[]
    csv_file=(out/'observables.csv').open('w',newline='');writer=csv.writer(csv_file)
    fields=['stage','cumulative_step','time_ps','target_temperature_K','temperature_K','total_density_g_cm3','volume_nm3','potential_kJ_mol','min_protein_PE_distance_nm','contact_residues','CA_RMSD_A']
    writer.writerow(fields);csv_file.flush()
    try:
        for stage in stages:
            sub=out/stage['name'];sub.mkdir()
            system=prm.createSystem(nonbondedMethod=PME,nonbondedCutoff=1*unit.nanometer,constraints=HBonds,rigidWater=True,ewaldErrorTolerance=.0005)
            system.setDefaultPeriodicBoxVectors(*box)
            nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
            charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
            if abs(charge)>1e-4:raise ValueError('Non-neutral system')
            pose=CustomExternalForce('0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2')
            for name in ('k','x0','y0','z0'):pose.addPerParticleParameter(name)
            anchor=np.array(positions.value_in_unit(unit.nanometer))
            for a in atoms:
                if a.element==element.hydrogen:continue
                k=stage['k_protein'] if a.residue.index<252 else stage['k_lower_PE'] if a.residue in lower else stage['k_upper_PE'] if a.residue.name=='PEO' else 0
                if k:pose.addParticle(a.index,[k,*anchor[a.index]])
            system.addForce(pose)
            temp=stage['schedule'][0][0]
            if stage['ensemble']=='NPT':
                pressure=MonteCarloBarostat(1*unit.bar,temp*unit.kelvin,25);pressure.setRandomNumberSeed(42);system.addForce(pressure)
            integrator=LangevinMiddleIntegrator(temp*unit.kelvin,1/unit.picosecond,2*unit.femtoseconds);integrator.setRandomNumberSeed(42)
            sim=Simulation(prm.topology,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'});sim.context.setPositions(positions);sim.context.setPeriodicBoxVectors(*box)
            if velocities is None:
                class Progress(MinimizationReporter):
                    def report(self,iteration,x,grad,args):
                        if iteration%20==0:print('Restrained minimization',iteration,args['system energy'],flush=True)
                        return False
                sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=100,reporter=Progress())
                minimized=sim.context.getState(getPositions=True)
                with (sub/'minimized.pdb').open('w') as f:PDBFile.writeFile(prm.topology,minimized.getPositions(),f)
                sim.context.setVelocitiesToTemperature(temp*unit.kelvin,42)
            else:sim.context.setVelocities(velocities)
            (sub/'system.xml').write_text(XmlSerializer.serialize(system));sim.reporters.append(DCDReporter(str(sub/'trajectory.dcd'),100))
            dof=3*system.getNumParticles()-system.getNumConstraints()-3;stage_rows=[]
            print('Starting',stage['name'],sum(n for _,n in stage['schedule']),'steps',flush=True)
            for temp,steps in stage['schedule']:
                integrator.setTemperature(temp*unit.kelvin)
                for offset in range(0,steps,100):
                    n=min(100,steps-offset);sim.step(n);total_steps+=n
                    s=sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);p=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
                    e=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);kinetic=s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                    vectors=np.array([v.value_in_unit(unit.nanometer) for v in s.getPeriodicBoxVectors()]);volume=float(np.linalg.det(vectors));lengths=np.diag(vectors)
                    temperature=2*kinetic/(dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
                    if not np.isfinite(p).all() or not all(math.isfinite(v) for v in (e,volume,temperature)) or temperature>1000 or volume<=0:raise ValueError('Unstable/non-finite state')
                    if min(lengths)<=2.1:raise ValueError('Periodic box too small for cutoff')
                    distance,contacts=surface_metrics(p,protein_heavy,pe_heavy,mapping,lengths)
                    if distance<.15:raise ValueError('Severe protein-PE heavy-atom overlap')
                    row=[stage['name'],total_steps,total_steps*.002,temp,temperature,density_g_cm3(mass,volume),volume,e,distance,len(contacts),aligned_rmsd(p[ca],reference[ca])*10]
                    rows.append(row);stage_rows.append(row);writer.writerow(row);csv_file.flush()
                    if sim.currentStep%500==0 or sim.currentStep==sum(n for _,n in stage['schedule']):print(stage['name'],'time',round(row[2],2),'ps T',round(temperature,1),'density',round(row[5],4),'volume',round(volume,2),flush=True)
            s=sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);positions=s.getPositions();velocities=s.getVelocities();box=s.getPeriodicBoxVectors()
            prm.topology.setPeriodicBoxVectors(box)
            sim.saveState(str(sub/'state.xml'));sim.saveCheckpoint(str(sub/'checkpoint.chk'))
            with (sub/'final.pdb').open('w') as f:PDBFile.writeFile(prm.topology,positions,f)
            stage_report={'name':stage['name'],'ensemble':stage['ensemble'],'executed_steps':sim.currentStep,'end_time_ps':total_steps*.002,'final_temperature_K':stage_rows[-1][4],'initial_recorded_density_g_cm3':stage_rows[0][5],'final_density_g_cm3':stage_rows[-1][5],'initial_recorded_volume_nm3':stage_rows[0][6],'final_volume_nm3':stage_rows[-1][6],'observed_volume_changes':sum(abs(b[6]-a[6])>1e-8 for a,b in zip(stage_rows,stage_rows[1:])),'position_restraints_k':{k:stage[k] for k in ('k_protein','k_lower_PE','k_upper_PE')},'equilibration_established':False}
            (sub/'stage-result.json').write_text(json.dumps(stage_report,indent=2));stage_reports.append(stage_report)
            del sim,integrator,system
        result={'status':'staged_preparation_completed_equilibrium_not_established','atoms':43288,'steps':total_steps,'time_ps':total_steps*.002,'timestep_fs':2,'net_charge_e':charge,'stages':stage_reports,
                'final_temperature_K':rows[-1][4],'final_total_density_g_cm3':rows[-1][5],'final_volume_nm3':rows[-1][6],'final_min_distance_nm':rows[-1][8],'final_contact_residues':contacts,
                'elapsed_seconds':time.monotonic()-started,'equilibration_established':False,'production_ready':False,
                'next_required':['Continue NPT until volume/density and energy are stable over adequately long windows','Check bulk-water and interfacial water structure, PE packing and protein structure','Extend restrained interface relaxation and release restraints in stages','Repeat independent orientations/seeds; assess finite patch and artificial-terminal effects'],
                'limitations':config['limitations']}
        (out/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result
    except Exception as e:
        (out/'failure.json').write_text(json.dumps({'status':'failed','executed_steps':total_steps,'error':f'{type(e).__name__}: {e}','completed_stages':stage_reports},indent=2));raise
    finally:csv_file.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bundle');p.add_argument('previous_run');p.add_argument('output');p.add_argument('--npt-ps',type=float,default=8);a=p.parse_args()
    if not 1<=a.npt_ps<=100:p.error('npt-ps must be between 1 and 100')
    run(a.bundle,a.previous_run,a.output,a.npt_ps)
