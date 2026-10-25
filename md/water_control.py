"""Small pure-water implementation controls matched to the actual PE source.

These tests do not validate PE interfaces, fibrinogen adsorption or equilibrium.
"""
import json
import math
import time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from openmm import System,XmlSerializer,unit,NonbondedForce,MonteCarloBarostat,LangevinMiddleIntegrator,Platform,Vec3
from openmm.app import Topology,ForceField,Simulation,PME,HBonds,element,DCDReporter
from workbench.config import ROOT,MODELS,fingerprint
from workbench.manager import atomic


def initial_box(water_mass=18.016,density=.86,seed=701):
 length=(512*water_mass*.00166053906660/density)**(1/3)
 centers=np.stack(np.meshgrid(*([(np.arange(8)+.5)*length/8]*3),indexing='ij'),axis=-1).reshape(-1,3)
 angle=np.deg2rad(104.52);template=np.array([[0.,0,0],[.09572,0,0],[.09572*np.cos(angle),.09572*np.sin(angle),0]])
 rotations=Rotation.random(512,random_state=np.random.default_rng(seed)).as_matrix();xyz=(np.einsum('nij,aj->nai',rotations,template)+centers[:,None,:]).reshape(-1,3)
 return length,xyz


def summarize(rows):
 npt=[r for r in rows if r['ensemble']=='NPT'];last=[r for r in npt if r['additional_NPT_ps']>=80];a=[r['density_g_cm3'] for r in last if r['additional_NPT_ps']<90];b=[r['density_g_cm3'] for r in last if r['additional_NPT_ps']>=90];mean=float(np.mean([r['density_g_cm3'] for r in last]));drift=float(abs(np.mean(b)/np.mean(a)-1)*100);temp=float(np.mean([r['temperature_K'] for r in last]));changes=int(np.count_nonzero(np.abs(np.diff([r['volume_nm3'] for r in npt]))>1e-8))
 return {'status':'technical_control_passed' if .90<mean<1.10 and drift<1 and abs(temp-310)<10 and changes>=10 else 'needs_review','last_20ps_mean_density_g_cm3':mean,'last_adjacent_10ps_mean_density_change_percent':drift,'last_20ps_mean_temperature_K':temp,'observed_volume_changes':changes,'equilibrium_established':False,'interpretation':'Wide implementation/sampler checks only. The liquid-like density band is not a calibrated TIP3P equation-of-state benchmark or a thermodynamic convergence test.'}


def run(source_job,out,root=ROOT):
 root=Path(root);source=root/'runs/workbench'/source_job;out=Path(out);out.mkdir(parents=True,exist_ok=False);cursor=json.loads((source/'checkpoint.json').read_text());config=json.loads((source/'config.json').read_text())
 if config['model_id']!='pe-core' or json.loads((source/'validation.json').read_text())['status']!='saved_output_cross_checks_passed':raise ValueError('Validated PE source required')
 path=source/cursor['system']
 if fingerprint(path)!=path.stem.split('-')[-1]:raise ValueError('Source System identity mismatch')
 original=XmlSerializer.deserialize(path.read_text());nb=next(f for f in original.getForces() if isinstance(f,NonbondedForce));source_water_atoms=[81275,81276,81277]
 # The source manifest fixes dry-solute count and water atom order. Verify against topology.
 from openmm.app import AmberPrmtopFile
 prm=AmberPrmtopFile(str(root/MODELS['pe-core']['prmtop']));first_water=next(r for r in prm.topology.residues() if r.name in ('HOH','WAT'));source_water_atoms=[a.index for a in first_water.atoms()]
 if [a.element.symbol for a in first_water.atoms()]!=['O','H','H']:raise ValueError('Unexpected source water order')
 params=[nb.getParticleParameters(i) for i in source_water_atoms];masses=[original.getParticleMass(i) for i in source_water_atoms];mass=sum(x.value_in_unit(unit.dalton) for x in masses);bar=next(f for f in original.getForces() if isinstance(f,MonteCarloBarostat))
 source_conditions={'source_job':source_job,'source_system_sha256':fingerprint(path),'barostat_frequency_steps':bar.getFrequency(),'applied_pressure_bar':bar.getDefaultPressure().value_in_unit(unit.bar),'barostat_temperature_K':bar.getDefaultTemperature().value_in_unit(unit.kelvin),'nonbonded_method':int(nb.getNonbondedMethod()),'cutoff_nm':nb.getCutoffDistance().value_in_unit(unit.nanometer),'dispersion_correction':nb.getUseDispersionCorrection(),'Ewald_error_tolerance':nb.getEwaldErrorTolerance(),'water_OHH_mass_Da':[m.value_in_unit(unit.dalton) for m in masses],'water_OHH_charge_e':[p[0].value_in_unit(unit.elementary_charge) for p in params],'water_OHH_sigma_nm':[p[1].value_in_unit(unit.nanometer) for p in params],'water_OHH_epsilon_kJ_mol':[p[2].value_in_unit(unit.kilojoule_per_mole) for p in params]};atomic(out/'source-conditions.json',source_conditions)
 results=[]
 for seed in [701,709]:
  run_dir=out/f'seed-{seed}';run_dir.mkdir();top=Topology();chain=top.addChain()
  for _ in range(512):
   r=top.addResidue('HOH',chain);o=top.addAtom('O',element.oxygen,r);h1=top.addAtom('H1',element.hydrogen,r);h2=top.addAtom('H2',element.hydrogen,r);top.addBond(o,h1);top.addBond(o,h2)
  length,xyz=initial_box(mass,seed=seed);box=[Vec3(length,0,0),Vec3(0,length,0),Vec3(0,0,length)]*unit.nanometer;top.setPeriodicBoxVectors(box);system=ForceField('amber14/tip3p.xml').createSystem(top,nonbondedMethod=PME,nonbondedCutoff=nb.getCutoffDistance(),constraints=HBonds,rigidWater=True,ewaldErrorTolerance=nb.getEwaldErrorTolerance());control_nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce));control_nb.setUseDispersionCorrection(nb.getUseDispersionCorrection())
  for atom in range(1536):
   q,sigma,eps=control_nb.getParticleParameters(atom);expected=params[atom%3]
   if abs(q.value_in_unit(unit.elementary_charge)-expected[0].value_in_unit(unit.elementary_charge))>1e-7 or abs(eps.value_in_unit(unit.kilojoule_per_mole)-expected[2].value_in_unit(unit.kilojoule_per_mole))>1e-6 or (eps.value_in_unit(unit.kilojoule_per_mole)>0 and abs(sigma.value_in_unit(unit.nanometer)-expected[1].value_in_unit(unit.nanometer))>1e-6):raise ValueError('Water XML nonbonded parameters differ from actual source')
   # Use exact source masses and NB values, including harmless zero-epsilon H sigma.
   system.setParticleMass(atom,masses[atom%3]);control_nb.setParticleParameters(atom,*expected)
  integrator=LangevinMiddleIntegrator(310,5,.002);integrator.setConstraintTolerance(1e-6);integrator.setRandomNumberSeed(seed);sim=Simulation(top,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'2'});sim.context.setPositions(xyz*unit.nanometer);sim.minimizeEnergy(maxIterations=500);sim.context.setVelocitiesToTemperature(310,seed);rows=[];dof=1536*3-system.getNumConstraints()-3
  for ensemble,steps in [('NVT',2500),('NPT',50000)]:
   if ensemble=='NPT':
    state=sim.context.getState(getPositions=True,getVelocities=True);step=sim.currentStep;clock=state.getTime();barostat=MonteCarloBarostat(1,310,25);barostat.setRandomNumberSeed(seed+100);system.addForce(barostat);integrator=LangevinMiddleIntegrator(310,5,.002);integrator.setConstraintTolerance(1e-6);integrator.setRandomNumberSeed(seed+1);sim=Simulation(top,system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'2'});sim.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors());sim.context.setPositions(state.getPositions());sim.context.setVelocities(state.getVelocities());sim.context.setTime(clock);sim.currentStep=step;sim.reporters.append(DCDReporter(str(run_dir/'npt.dcd'),2500))
   for done in range(250,steps+1,250):
    sim.step(250);s=sim.context.getState(getEnergy=True,getPositions=True);v=np.array([a.value_in_unit(unit.nanometer) for a in s.getPeriodicBoxVectors()]);volume=float(np.linalg.det(v));temperature=2*s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(dof*.00831446261815324);energy=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);density=512*mass*.00166053906660/volume
    if not np.isfinite(s.getPositions(asNumpy=True)._value).all() or not all(math.isfinite(x) for x in [temperature,energy,density]) or not 0<temperature<600 or np.diag(v).min()<=2*source_conditions['cutoff_nm']:raise ValueError('Water control numerical/cell safety failure')
    row={'seed':seed,'ensemble':ensemble,'context_time_ps':s.getTime().value_in_unit(unit.picosecond),'additional_NPT_ps':done*.002 if ensemble=='NPT' else None,'temperature_K':temperature,'density_g_cm3':density,'volume_nm3':volume,'potential_kJ_mol':energy,'pressure_bar':barostat.computeCurrentPressure(sim.context).value_in_unit(unit.bar) if ensemble=='NPT' and done%1000==0 else None};rows.append(row)
    with (run_dir/'observations.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    if done%2500==0:sim.saveCheckpoint(str(run_dir/'checkpoint.chk'));atomic(out/'progress.json',{'status':'running','seed':seed,'last':row})
  result=summarize(rows);result.update(seed=seed,water_molecules=512,atoms=1536,NVT_ps=5,NPT_ps=100,source_nonbonded_parameters_matched=True,RNG_note='Independent control seeds. NVT-to-NPT uses a new Context with documented thermostat/barostat seeds; protein trajectory is unaffected.');atomic(run_dir/'result.json',result);(run_dir/'system.xml').write_text(XmlSerializer.serialize(system));(run_dir/'final-state.xml').write_text(XmlSerializer.serialize(sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True)));results.append(result)
 atomic(out/'result.json',{'status':'technical_controls_passed' if all(x['status']=='technical_control_passed' for x in results) else 'needs_review','source_conditions':source_conditions,'replicas':results,'equilibrium_established':False,'scope':'Pure-water implementation controls only. Does not validate PE/PVDF/PTFE surfaces, whole-system pressure equilibration, adsorption or experimental retention.'});return json.loads((out/'result.json').read_text())
