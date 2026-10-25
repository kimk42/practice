"""Prepared-model NVT/NPT engine with checkpointed controls and explicit screens."""
import csv
import hashlib
import json
import math
from pathlib import Path
import resource
import time
import numpy as np
from scipy.spatial import cKDTree
from openmm import unit,XmlSerializer,Platform,LangevinMiddleIntegrator,MonteCarloBarostat,CustomExternalForce,NonbondedForce
from openmm.app import AmberPrmtopFile,Simulation,DCDReporter,PDBxFile,element
from md.protein_pilot import aligned_rmsd
from md.minimize_large_pe import geometry
from .config import ROOT,MODELS,fingerprint
from .manager import atomic


class Cancelled(Exception):pass
class ReviewRequired(Exception):pass


def density_screen(rows,spec):
    stages=[]
    for name in dict.fromkeys(x['stage'] for x in rows if x['ensemble']=='NPT'):
        sample=[x for x in rows if x['stage']==name];t=np.array([x['context_time_ps'] for x in sample]);rho=np.array([x['density_g_cm3'] for x in sample]);u=np.array([x['potential_kJ_mol'] for x in sample])
        stages.append({'stage':name,'observations':len(sample),'observed_volume_changes':int(np.count_nonzero(np.abs(np.diff([x['volume_nm3'] for x in sample]))>1e-6)) if 'volume_nm3' in sample[0] else None,'mean_density_g_cm3':float(rho.mean()),'density_endpoint_change_percent':float((rho[-1]/rho[0]-1)*100),'mean_temperature_K':float(np.mean([x['temperature_K'] for x in sample])),'relative_potential_slope_percent_per_ps':float(np.polyfit(t-t[0],u,1)[0]/max(abs(u.mean()),1)*100)})
    check={'blocks':stages,'screen_passed':False,'equilibrium_established':False,'interpretation':spec['interpretation']}
    if len(stages)<spec['minimum_blocks']:return check
    a,b=stages[-2:];change=abs(b['mean_density_g_cm3']/a['mean_density_g_cm3']-1)*100
    check['last_adjacent_density_mean_change_percent']=change
    check['screen_passed']=bool(change<spec['maximum_adjacent_mean_density_change_percent'] and abs(b['density_endpoint_change_percent'])<spec['maximum_last_block_density_drift_percent'] and abs(b['mean_temperature_K']-310)<spec['maximum_temperature_mean_deviation_K'] and abs(b['relative_potential_slope_percent_per_ps'])<spec['maximum_relative_potential_slope_percent_per_ps'])
    if spec.get('minimum_last_block_observed_volume_changes',0):
        check['screen_passed']=check['screen_passed'] and b['observed_volume_changes'] is not None and b['observed_volume_changes']>=spec['minimum_last_block_observed_volume_changes']
    return check


def inspect_runtime(sim,system,integrator):
 parameters={str(k):float(v) for k,v in sim.context.getParameters().items()};bars=[f for f in system.getForces() if isinstance(f,MonteCarloBarostat)];positive=sum(abs(f.getParticleParameters(i)[1][0])>1e-12 for f in system.getForces() if isinstance(f,CustomExternalForce) for i in range(f.getNumParticles()))
 return {'effective_context_parameters':parameters,'thermostat_temperature_K':integrator.getTemperature().value_in_unit(unit.kelvin),'friction_per_ps':integrator.getFriction().value_in_unit(unit.picosecond**-1),'timestep_fs':integrator.getStepSize().value_in_unit(unit.femtosecond),'thermostat_seed_in_serialized_integrator':integrator.getRandomNumberSeed(),'nonzero_positional_support_particles':positive,'barostat_frequencies_steps':[b.getFrequency() for b in bars],'barostat_seeds_in_serialized_system':[b.getRandomNumberSeed() for b in bars],'integration_force_groups':int(integrator.getIntegrationForceGroups()),'note':'Seeds identify serialized objects; continued stochastic streams are loaded from the binary checkpoint, not initialized here.'}


def validate_continuation_runtime(conditions):
 p=conditions['effective_context_parameters']
 if conditions['barostat_frequencies_steps']!=[25] or conditions['nonzero_positional_support_particles'] or abs(p.get(MonteCarloBarostat.Pressure(),float('inf'))-1)>1e-9 or abs(p.get(MonteCarloBarostat.Temperature(),float('inf'))-310)>1e-9:raise ReviewRequired('Effective checkpoint pressure/temperature/support settings differ from continuation protocol')


class Engine:
 def __init__(self,job,root=ROOT):
  self.job=Path(job);self.root=Path(root);self.cfg=json.loads((self.job/'config.json').read_text());self.model=MODELS[self.cfg['model_id']];self.started=time.monotonic();self.rows=[];self.segments=[];self.sim=None;self.stage_index=0;self.stage_step=0;self.slot=0;self.system_path=None;self.base_k=None
  self.objects=self.job/'objects';self.objects.mkdir(exist_ok=True);(self.job/'trajectories').mkdir(exist_ok=True)
 def state(self,status,**data):atomic(self.job/'state.json',{'status':status,'updated_unix':time.time(),'elapsed_wall_seconds':time.monotonic()-self.started,**data})
 def object(self,kind,serialized):
  data=serialized.encode();digest=hashlib.sha256(data).hexdigest();path=self.objects/f'{kind}-{digest}.xml'
  if not path.exists():path.write_bytes(data)
  return path
 def start(self):
  self.state('validating',message='Checking prepared source identity and runtime')
  m=self.model;source=self.root/m['source'];prm=self.root/m['prmtop'];manifest=json.loads((self.root/m['manifest']).read_text())
  expected=manifest['sha256']['solvated.prmtop' if self.cfg['model_id']=='pe-core' else 'complex.prmtop']
  if fingerprint(prm)!=expected:raise ReviewRequired('Prepared topology checksum differs')
  self.prm=AmberPrmtopFile(str(prm));original_system=XmlSerializer.deserialize((source/'system.xml').read_text());original_state=XmlSerializer.deserialize((source/'state.xml').read_text())
  if self.cfg['model_id']=='pe-core':
   r=json.loads((source/'result.json').read_text());v=json.loads((source/'validation.json').read_text())
   if r['status']!='100K_hold_MD_completed_numerical_checks_passed' or v['status']!='independent_100K_hold_output_checks_passed':raise ReviewRequired('Source core preparation is not verified')
   for name in ['state.xml','system.xml']:
    if fingerprint(source/name)!=r['sha256'][name]:raise ReviewRequired('Reviewed starting state/system checksum differs')
   self.integrator=XmlSerializer.deserialize((source/'integrator-final.xml').read_text())
  else:
   v=json.loads((self.root/m['source_parent']/'validation.json').read_text());s=next(x for x in v['stages'] if x['stage']==source.name)
   if fingerprint(source/'state.xml')!=s['state_sha256']:raise ReviewRequired('Reviewed fragment state differs')
   self.integrator=LangevinMiddleIntegrator(310*unit.kelvin,1/unit.picosecond,2*unit.femtoseconds);self.integrator.setRandomNumberSeed(42)
  self.atoms=list(self.prm.topology.atoms());self.original_support=next(f for f in original_system.getForces() if isinstance(f,CustomExternalForce));self.base_k=[self.original_support.getParticleParameters(i)[1][0] for i in range(self.original_support.getNumParticles())]
  if any(original_system.getParticleMass(i).value_in_unit(unit.dalton)<=0 for i in range(original_system.getNumParticles())):raise ReviewRequired('Massless particles are not supported')
  if self.cfg['model_id']=='pe-core':
   maps=json.loads((self.root/'data/fibrinogen-core03/manifest.json').read_text())['residue_mapping'];native=lambda a:a.residue.index<len(maps) and maps[a.residue.index]['source_label_seq'] is not None
  else:native=lambda a:a.residue.index<252
  self.ph=np.array([a.index for a in self.atoms if native(a) and a.element!=element.hydrogen]);self.ca=np.array([a.index for a in self.atoms if native(a) and a.name=='CA']);self.pe=np.array([a.index for a in self.atoms if a.residue.name=='PEO' and a.element==element.carbon]);self.heavy=np.array([a.index for a in self.atoms if a.element!=element.hydrogen]);self.water=np.array([a.index for a in self.atoms if a.residue.name in ('HOH','WAT') and a.element==element.oxygen]);self.dry_n=81275 if self.cfg['model_id']=='pe-core' else 7449
  if len(self.ca)!=m['native_residues'] or original_system.getNumParticles()!=m['atoms']:raise ReviewRequired('Prepared model identity counts differ')
  c=[original_system.getConstraintParameters(i) for i in range(original_system.getNumConstraints())];self.ia=np.array([int(x[0]) for x in c]);self.ib=np.array([int(x[1]) for x in c]);self.distances=np.array([x[2].value_in_unit(unit.nanometer) for x in c]);del c
  self.mass=sum(original_system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(m['atoms']));self.dof=3*m['atoms']-len(self.ia)-3
  provenance={'model':m,'source_system_sha256':fingerprint(source/'system.xml'),'source_state_sha256':fingerprint(source/'state.xml'),'source_checkpoint_sha256':fingerprint(source/'checkpoint.chk'),'prmtop_sha256':expected,'engine_sha256':fingerprint(Path(__file__)),'openmm_platform':'CPU','threads':4,'plan':self.cfg['plan']}
  atomic(self.job/'provenance.json',provenance);(self.job/'engine-source.py').write_bytes(Path(__file__).read_bytes())
  resume=self.cfg.get('resume_from_job');continuation=self.cfg.get('continuation_from_job') if not resume else None
  if resume or continuation:
   parent=self.job.parent/(resume or continuation);cursor=json.loads((parent/'checkpoint.json').read_text());self.stage_index=cursor['stage_index'];self.stage_step=cursor['stage_step']
   for key in ['system','integrator']:
    obj=parent/cursor[key]
    if obj.parent!=parent/'objects' or fingerprint(obj)!=obj.stem.split('-')[-1]:raise ReviewRequired('Resume object checksum/path mismatch')
   self.system=XmlSerializer.deserialize((parent/cursor['system']).read_text());self.integrator=XmlSerializer.deserialize((parent/cursor['integrator']).read_text());cp=parent/cursor['checkpoint']
   if fingerprint(cp)!=cursor['sha256']:raise ReviewRequired('Resume checkpoint checksum mismatch')
   self.reference=np.load(parent/'reference.npy');self.origin_time=cursor['origin_time_ps'];self.origin_step=cursor['origin_context_step']
   if not continuation and (parent/'observations.jsonl').exists():
    for line in (parent/'observations.jsonl').read_text().splitlines():
     try:row=json.loads(line)
     except ValueError:break
     if row['context_step']<=cursor['context_step']:self.rows.append(row)
   atomic(self.job/'resume-provenance.json',{'parent_job':resume or continuation,'valid_parent_context_step':cursor['context_step'],'parent_trajectories':str(parent/'trajectories'),'note':'Parent files are preserved; observations after the saved checkpoint are not carried into this resumed lineage.'})
   if continuation:
    parent_cfg=json.loads((parent/'config.json').read_text());v=json.loads((parent/'validation.json').read_text());parent_status=json.loads((parent/'state.json').read_text())['status']
    if parent_status not in {'completed','needs_review'} or parent_cfg['model_id']!=self.cfg['model_id'] or v['status']!='saved_output_cross_checks_passed':raise ReviewRequired('Continuation source validation differs')
    support=next(f for f in self.system.getForces() if isinstance(f,CustomExternalForce));bars=[f for f in self.system.getForces() if isinstance(f,MonteCarloBarostat)]
    if any(abs(support.getParticleParameters(i)[1][0])>1e-12 for i in range(support.getNumParticles())) or len(bars)!=1:raise ReviewRequired('Continuation requires unrestrained NPT source')
    if abs(self.integrator.getStepSize().value_in_unit(unit.femtosecond)-2)>1e-9 or abs(self.integrator.getTemperature().value_in_unit(unit.kelvin)-310)>1e-9 or abs(bars[0].getDefaultPressure().value_in_unit(unit.bar)-1)>1e-9:raise ReviewRequired('Continuation ensemble settings differ')
    original_state=XmlSerializer.deserialize((parent/'final-state.xml').read_text());self.reference=original_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer);self.origin_time=cursor['context_time_ps'];self.origin_step=cursor['context_step'];self.stage_index=0;self.stage_step=0
    atomic(self.job/'continuation-provenance.json',{'parent_job':continuation,'parent_checkpoint_sha256':cursor['sha256'],'parent_final_state_sha256':fingerprint(parent/'final-state.xml'),'starting_context_step':self.origin_step,'starting_context_time_ps':self.origin_time,'velocities_reinitialized':False,'thermostat_or_barostat_RNG_reseeded':False,'reference_note':'New diagnostic reference at parent final state; original preparation reference retained in parent.'})
  else:
   self.system=original_system;cp=source/'checkpoint.chk';self.reference=original_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer);self.origin_time=original_state.getTime().value_in_unit(unit.picosecond);self.origin_step=original_state.getStepCount()
  self.sim=Simulation(self.prm.topology,self.system,self.integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'});self.sim.loadCheckpoint(str(cp))
  loaded=self.sim.context.getState(getPositions=True,getVelocities=True)
  if (resume or continuation) and (loaded.getStepCount()!=cursor['context_step'] or abs(loaded.getTime().value_in_unit(unit.picosecond)-cursor['context_time_ps'])>1e-8):raise ReviewRequired('Resume checkpoint clock differs from saved cursor')
  if not resume:
   if not np.allclose(loaded.getPositions(asNumpy=True).value_in_unit(unit.nanometer),original_state.getPositions(asNumpy=True).value_in_unit(unit.nanometer),atol=1e-10,rtol=0) or not np.allclose(loaded.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond),original_state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer/unit.picosecond),atol=1e-10,rtol=0):raise ReviewRequired('Checkpoint and reviewed state disagree')
  conditions=inspect_runtime(self.sim,self.system,self.integrator);atomic(self.job/'runtime-conditions.json',conditions)
  if self.cfg.get('continuation_from_job'):validate_continuation_runtime(conditions)
  np.save(self.job/'reference.npy',self.reference)
  with (self.job/'observations.jsonl').open('w') as f:
   for row in self.rows:f.write(json.dumps(row)+'\n')
  self.system_path=self.object('system',XmlSerializer.serialize(self.system));self.open_trajectory('source-continuation')
 def open_trajectory(self,reason):
  self.sim.reporters.clear();interval=self.cfg['trajectory_interval_steps'];step=self.sim.currentStep;clock=self.sim.context.getState().getTime().value_in_unit(unit.picosecond);dt=self.integrator.getStepSize().value_in_unit(unit.picosecond)
  name=f'trajectories/segment-{len(self.segments)+1:02d}.dcd';self.sim.reporters.append(DCDReporter(str(self.job/name),interval))
  self.segments.append({'path':name,'dt_ps':dt,'context_step_offset':step-step%interval,'context_time_offset_ps':clock-(step%interval)*dt,'reason':reason});atomic(self.job/'trajectory-index.json',self.segments)
 def checkpoint(self):
  self.slot=1-self.slot;name=f'checkpoint-{self.slot}.chk';self.sim.saveCheckpoint(str(self.job/(name+'.tmp')));(self.job/(name+'.tmp')).replace(self.job/name)
  ip=self.object('integrator',XmlSerializer.serialize(self.integrator));clock=self.sim.context.getState().getTime().value_in_unit(unit.picosecond)
  atomic(self.job/'checkpoint.json',{'checkpoint':name,'sha256':fingerprint(self.job/name),'system':str(self.system_path.relative_to(self.job)),'integrator':str(ip.relative_to(self.job)),'context_step':self.sim.currentStep,'context_time_ps':clock,'origin_time_ps':self.origin_time,'origin_context_step':self.origin_step,'stage_index':self.stage_index,'stage_step':self.stage_step})
 def controls(self):
  path=self.job/'control.json'
  if not path.exists():return
  action=json.loads(path.read_text()).get('action')
  if action=='cancel':self.checkpoint();raise Cancelled('Cancelled at a verified MD chunk boundary')
  if action=='pause':
   self.checkpoint();self.state('paused',message='Checkpoint saved; waiting for resume or cancel',last=self.rows[-1] if self.rows else None)
   while True:
    time.sleep(.25);action=json.loads(path.read_text()).get('action')
    if action=='cancel':raise Cancelled('Cancelled while paused')
    if action=='resume':break
 def change_context(self,stage):
  old=self.sim.context.getState(getPositions=True,getVelocities=True);step=self.sim.currentStep;clock=old.getTime();box=old.getPeriodicBoxVectors()
  self.sim.reporters.clear()
  for i in reversed(range(self.system.getNumForces())):
   if isinstance(self.system.getForce(i),MonteCarloBarostat):self.system.removeForce(i)
  support=next(f for f in self.system.getForces() if isinstance(f,CustomExternalForce))
  if any(abs(support.getParticleParameters(i)[1][0])>1e-12 for i in range(support.getNumParticles())):raise ReviewRequired('Cartesian position restraints must be zero before isotropic NPT')
  seed=self.cfg['seed']+self.stage_index;bar=MonteCarloBarostat(stage['pressure_bar']*unit.bar,310*unit.kelvin,25);bar.setRandomNumberSeed(seed+100);self.system.addForce(bar);self.system.setDefaultPeriodicBoxVectors(*box)
  self.integrator=LangevinMiddleIntegrator(310*unit.kelvin,self.model['friction_per_ps']/unit.picosecond,stage['dt_fs']*unit.femtoseconds);self.integrator.setConstraintTolerance(1e-6);self.integrator.setRandomNumberSeed(seed)
  self.sim=Simulation(self.prm.topology,self.system,self.integrator,Platform.getPlatformByName('CPU'),{'Threads':'4'});self.sim.context.setPeriodicBoxVectors(*box);self.sim.context.setPositions(old.getPositions());self.sim.context.setVelocities(old.getVelocities());self.sim.context.setTime(clock);self.sim.currentStep=step
  self.system_path=self.object('system',XmlSerializer.serialize(self.system));self.open_trajectory('Ensemble/timestep change; positions and velocities continuous, thermostat/barostat RNG streams restarted')
 def apply(self,spec):
  self.integrator.setTemperature(spec['target_K']*unit.kelvin);force=next(f for f in self.system.getForces() if isinstance(f,CustomExternalForce));changed=False
  for i,k in enumerate(self.base_k):
   atom,params=force.getParticleParameters(i);new=k*spec['restraint_scale']
   if abs(params[0]-new)>1e-12:
    params=list(params);params[0]=new;force.setParticleParameters(i,atom,params);changed=True
  if changed:
   force.updateParametersInContext(self.sim.context);self.system_path=self.object('system',XmlSerializer.serialize(self.system))
 def sample(self,stage,spec):
  s=self.sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);xyz=s.getPositions(asNumpy=True).value_in_unit(unit.nanometer);box=np.array([v.value_in_unit(unit.nanometer) for v in s.getPeriodicBoxVectors()]);lengths=np.diag(box);volume=float(np.linalg.det(box));e=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);ke=s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole);temp=2*ke/(self.dof*unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole/unit.kelvin))
  if not np.allclose(box,np.diag(lengths)) or min(lengths)<=2.1 or not np.isfinite(xyz).all() or not all(math.isfinite(v) for v in (e,ke,temp,volume)) or not 0<=temp<=600:raise ReviewRequired('Nonfinite state, runaway temperature or unsupported cell')
  d=xyz[self.ia]-xyz[self.ib];d-=lengths*np.round(d/lengths);constraint=float(np.max(np.abs(np.linalg.norm(d,axis=1)-self.distances)/self.distances));rms=aligned_rmsd(xyz[self.ca],self.reference[self.ca])*10
  distances=cKDTree(np.mod(xyz[self.pe],lengths),boxsize=lengths).query(np.mod(xyz[self.ph],lengths))[0]
  if constraint>1e-4 or rms>20 or distances.min()<.15:raise ReviewRequired('Constraint, backbone displacement or protein/PE hard-overlap safety screen failed')
  row={'stage':stage['name'],'ensemble':stage['ensemble'],'stage_index':self.stage_index,'stage_step':self.stage_step,'context_step':self.sim.currentStep,'context_time_ps':s.getTime().value_in_unit(unit.picosecond),'additional_time_ps':s.getTime().value_in_unit(unit.picosecond)-self.origin_time,'target_K':spec['target_K'],'temperature_K':temp,'potential_kJ_mol':e,'kinetic_kJ_mol':ke,'volume_nm3':volume,'density_g_cm3':self.mass*1.66053906660e-3/volume,'restraint_scale':spec['restraint_scale'],'CA_RMSD_from_pipeline_start_A':rms,'max_relative_constraint_error':constraint,'minimum_observed_protein_PE_nm':float(distances.min()),'contact_native_residue_count':len({self.atoms[int(self.ph[i])].residue.index for i in np.flatnonzero(distances<.45)})}
  if self.cfg.get('pressure_interval_steps'):
   row['pressure_bar']=None
   if not self.rows or self.sim.currentStep%self.cfg['pressure_interval_steps']==0:
    bar=next(f for f in self.system.getForces() if isinstance(f,MonteCarloBarostat));pressure=bar.computeCurrentPressure(self.sim.context).value_in_unit(unit.bar)
    if not math.isfinite(pressure):raise ReviewRequired('Nonfinite instantaneous pressure')
    row['pressure_bar']=pressure
  self.rows.append(row)
  with (self.job/'observations.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
  self.state('running',stage=stage['name'],stage_index=self.stage_index,stage_step=self.stage_step,last=row,maximum_planned_steps=self.cfg['maximum_planned_steps'],message='Executed MD; equilibrium and adsorption are not established')
  print(json.dumps(row),flush=True);return s
 def run(self):
  self.start();last=None;screen=None;stages=self.cfg['plan']
  while self.stage_index<len(stages):
   stage=stages[self.stage_index]
   if stage['ensemble']=='NPT' and not any(isinstance(f,MonteCarloBarostat) for f in self.system.getForces()):self.change_context(stage)
   total=sum(s['steps'] for s in stage['schedule'])
   if self.stage_step==0:
    self.apply(stage['schedule'][0]);self.checkpoint()
   for spec in stage['schedule']:
    start=sum(s['steps'] for s in stage['schedule'][:stage['schedule'].index(spec)]);end=start+spec['steps']
    if self.stage_step>=end:continue
    self.apply(spec)
    if not self.rows:last=self.sample(stage,spec)
    while self.stage_step<end:
     self.controls();n=min(self.cfg['observation_interval_steps'],end-self.stage_step);self.sim.step(n);self.stage_step+=n
     try:last=self.sample(stage,spec)
     except ReviewRequired:
      bad=self.sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);(self.job/'rejected-state.xml').write_text(XmlSerializer.serialize(bad));raise
     if self.sim.currentStep%self.cfg['checkpoint_interval_steps']==0:self.checkpoint()
   self.checkpoint()
   if stage['name']=='310K-supported-check':
    tail=[x for x in self.rows if x['stage']==stage['name']][-20:]
    if abs(np.mean([x['temperature_K'] for x in tail])-310)>10:raise ReviewRequired('310 K temperature preparation screen failed; no restraint release or NPT started')
   if stage['ensemble']=='NPT':
    screen=density_screen(self.rows,self.cfg['density_screen']);atomic(self.job/'density-screen.json',screen)
    if screen['screen_passed']:break
   self.stage_index+=1;self.stage_step=0
  if last is None:last=self.sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True)
  self.checkpoint();final=self.sim.context.getState(getPositions=True,getVelocities=True,getEnergy=True);(self.job/'final-state.xml').write_text(XmlSerializer.serialize(final));self.prm.topology.setPeriodicBoxVectors(final.getPeriodicBoxVectors())
  with (self.job/'final.cif').open('w') as f:PDBxFile.writeFile(self.prm.topology,final.getPositions(),f)
  nb=next(f for f in self.system.getForces() if isinstance(f,NonbondedForce));exceptions=set();n=self.model['atoms']
  for i in range(nb.getNumExceptions()):
   a,b,*_=nb.getExceptionParameters(i);a,b=sorted((int(a),int(b)));exceptions.add(a*n+b)
  xyz=final.getPositions(asNumpy=True).value_in_unit(unit.nanometer);box=np.array([v.value_in_unit(unit.nanometer) for v in final.getPeriodicBoxVectors()]).diagonal();geo=geometry(xyz,box,self.heavy,self.water,self.heavy[self.heavy<self.dry_n],exceptions)
  if geo['nonexcluded_heavy_pairs_under_0_12nm'] or geo['nonexcluded_allatom_pairs_under_0_08nm']:raise ReviewRequired('Final nonexception hard-overlap screen failed')
  status='completed' if screen is None or screen['screen_passed'] else 'needs_review'
  result={'status':status,'model_id':self.cfg['model_id'],'atoms':n,'final':self.rows[-1],'observations':len(self.rows),'density_screen':screen,'geometry_after':geo,'additional_time_ps':self.rows[-1]['additional_time_ps'],'context_time_ps':final.getTime().value_in_unit(unit.picosecond),'context_step':self.sim.currentStep,'equilibrium_established':False,'production_ready':False,'preparation_screen_passed':status=='completed','elapsed_wall_seconds':time.monotonic()-self.started,'peak_process_RSS_GiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,'scope':'Finite PE exploratory preparation. No 310 K/interface equilibrium, binding affinity or material ranking is inferred.'}
  pressure=[r['pressure_bar'] for r in self.rows if r.get('pressure_bar') is not None]
  if pressure:
   result['pressure_diagnostics']={'samples':len(pressure),'mean_bar':float(np.mean(pressure)),'std_bar':float(np.std(pressure)),'minimum_bar':min(pressure),'maximum_bar':max(pressure),'method':'OpenMM MonteCarloBarostat.computeCurrentPressure molecular virial, sparse instantaneous samples','equilibrium_established':False,'interpretation':'Correlated, sparse short-preparation samples; no converged mean pressure or uncertainty inferred.'}
  atomic(self.job/'result.json',result)
  with (self.job/'observables.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(self.rows[0]));w.writeheader();w.writerows(self.rows)
  return result
