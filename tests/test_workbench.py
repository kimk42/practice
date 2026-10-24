import json
from pathlib import Path
import tempfile
import threading
import unittest
from http.server import ThreadingHTTPServer
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from workbench.config import make_config
from workbench.manager import Manager,atomic
from workbench.server import application

class ConfigTests(unittest.TestCase):
 def test_preparation_releases_before_pressure(self):
  c=make_config({'model_id':'pe-core','protocol_id':'prepare-310'})
  self.assertEqual(c['maximum_planned_time_ps'],19.1)
  self.assertEqual(c['plan'][3]['schedule'][-1]['restraint_scale'],0)
  self.assertTrue(all(s['restraint_scale']==0 for stage in c['plan'] if stage['ensemble']=='NPT' for s in stage['schedule']))
  self.assertFalse(c['equilibrium_established'])
 def test_rejects_unsupported_inputs(self):
  for request in [{'model_id':'PVDF','protocol_id':'prepare-310'},{'model_id':'pe-core','protocol_id':'prepare-310','steps':100},{'model_id':'pe-fragment','protocol_id':'nvt-check','steps':True},{'model_id':'pe-core','protocol_id':'prepare-310','max_npt_blocks':99},{'model_id':'pe-core','protocol_id':'prepare-310','command':'echo hi'}]:
   with self.subTest(request=request),self.assertRaises(ValueError):make_config(request)

class ApiTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);(self.root/'index.html').write_text('public');(self.root/'.env').write_text('private')
  self.manager=Manager(self.root,check_runtime=False);self.server=ThreadingHTTPServer(('127.0.0.1',0),application(self.manager));self.thread=threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start();self.base=f'http://127.0.0.1:{self.server.server_port}'
  self.token=self.get('/api/capabilities')['token']
 def tearDown(self):self.server.shutdown();self.server.server_close();self.thread.join();self.tmp.cleanup()
 def get(self,path):
  with urlopen(self.base+path) as f:return json.load(f)
 def post(self,path,data,token=True,origin=None):
  headers={'Content-Type':'application/json'}
  if token:headers['X-Workbench-Token']=self.token
  if origin:headers['Origin']=origin
  return urlopen(Request(self.base+path,data=json.dumps(data).encode(),headers=headers))
 def test_post_protection_and_validation(self):
  req={'model_id':'pe-fragment','protocol_id':'nvt-check'}
  for token,origin in [(False,None),(True,'https://unrelated.example')]:
   with self.assertRaises(HTTPError) as e:self.post('/api/preview',req,token,origin)
   self.assertEqual(e.exception.code,403)
  with self.post('/api/preview',req) as f:self.assertEqual(json.load(f)['maximum_planned_steps'],100)
  with self.assertRaises(HTTPError) as e:self.post('/api/jobs',req)
  self.assertEqual(e.exception.code,400)
 def test_private_paths_and_traversal(self):
  for path in ['/.env','/%2eenv','/runs/anything','/../.env','/api/jobs/../../.env']:
   with self.subTest(path=path),self.assertRaises(HTTPError):urlopen(self.base+path)
 def test_persistent_cancelled_job_resume_metadata(self):
  id='a'*32;p=self.manager.jobs/id;p.mkdir();cfg=make_config({'model_id':'pe-fragment','protocol_id':'nvt-check'});cfg.update(id=id,created_unix=1);atomic(p/'config.json',cfg);atomic(p/'state.json',{'status':'cancelled'})
  self.assertEqual(self.get('/api/jobs')[0]['status'],'cancelled')
  with self.assertRaises(ValueError):self.manager.directory('../other')

class ScreenTests(unittest.TestCase):
 def test_screen_requires_complete_blocks_and_retains_caveat(self):
  from workbench.engine import density_screen
  spec=make_config({'model_id':'pe-core','protocol_id':'prepare-310'})['density_screen'];rows=[]
  for b in range(3):
   for i in range(4):rows.append({'stage':str(b),'ensemble':'NPT','context_time_ps':b*2.5+i*.5,'density_g_cm3':.99,'temperature_K':310,'potential_kJ_mol':-10000})
   result=density_screen(rows,spec);self.assertEqual(result['screen_passed'],b==2);self.assertFalse(result['equilibrium_established'])
  for r in rows[-4:]:r['temperature_K']=330
  self.assertFalse(density_screen(rows,spec)['screen_passed'])

class EngineTransitionTests(unittest.TestCase):
 def test_restraint_release_and_npt_context_preserve_state(self):
  import numpy as np
  from types import SimpleNamespace
  from unittest.mock import Mock
  from openmm import System,CustomExternalForce,NonbondedForce,LangevinMiddleIntegrator,Platform,unit,MonteCarloBarostat
  from openmm.app import Topology,Simulation,element
  from workbench.engine import Engine,ReviewRequired
  with tempfile.TemporaryDirectory() as tmp:
   e=Engine.__new__(Engine);e.objects=Path(tmp);e.cfg={'seed':45};e.stage_index=4;e.model={'friction_per_ps':5};e.base_k=[1000];e.open_trajectory=Mock()
   system=System();top=Topology();chain=top.addChain();res=top.addResidue('TST',chain)
   nb=NonbondedForce();nb.setNonbondedMethod(NonbondedForce.CutoffPeriodic);nb.setCutoffDistance(1)
   for _ in range(2):system.addParticle(12);top.addAtom('C',element.carbon,res);nb.addParticle(0,.3,.1)
   system.addForce(nb);support=CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
   for name in ['k','x0','y0','z0']:support.addPerParticleParameter(name)
   support.addParticle(0,[1000,0,0,0]);system.addForce(support);box=[(3,0,0),(0,3,0),(0,0,3)]*unit.nanometer;system.setDefaultPeriodicBoxVectors(*box);top.setPeriodicBoxVectors(box)
   e.system=system;e.prm=SimpleNamespace(topology=top);e.integrator=LangevinMiddleIntegrator(310,5,.001);e.sim=Simulation(top,system,e.integrator,Platform.getPlatformByName('CPU'),{'Threads':'1'})
   e.sim.context.setPositions([[.1,.1,.1],[.5,.1,.1]]);e.sim.context.setVelocities([[.01,0,0],[0,.01,0]]);e.sim.context.setTime(7.1);e.sim.currentStep=7100
   stage={'pressure_bar':1,'dt_fs':2}
   with self.assertRaises(ReviewRequired):e.change_context(stage)
   e.apply({'target_K':310,'restraint_scale':.1});self.assertEqual(support.getParticleParameters(0)[1][0],100)
   e.apply({'target_K':310,'restraint_scale':0});old=e.sim.context.getState(getPositions=True,getVelocities=True)
   e.change_context(stage);now=e.sim.context.getState(getPositions=True,getVelocities=True)
   self.assertTrue(np.allclose(old.getPositions(asNumpy=True)._value,now.getPositions(asNumpy=True)._value));self.assertTrue(np.allclose(old.getVelocities(asNumpy=True)._value,now.getVelocities(asNumpy=True)._value))
   self.assertEqual(e.sim.currentStep,7100);self.assertAlmostEqual(now.getTime().value_in_unit(unit.picosecond),7.1);self.assertTrue(any(isinstance(f,MonteCarloBarostat) for f in e.system.getForces()));self.assertAlmostEqual(e.integrator.getStepSize().value_in_unit(unit.femtosecond),2)

class ContinuationTests(unittest.TestCase):
 def test_bounded_continuation_has_no_heating_or_restraints(self):
  c=make_config({'model_id':'pe-core','protocol_id':'continue-npt','source_job_id':'a'*32})
  self.assertEqual(c['maximum_planned_time_ps'],15)
  self.assertEqual(c['pressure_interval_steps'],250)
  self.assertTrue(all(s['ensemble']=='NPT' and s['schedule'][0]['restraint_scale']==0 for s in c['plan']))
  for field,value in [('source_job_id','../other'),('blocks',4),('blocks',True)]:
   req={'model_id':'pe-core','protocol_id':'continue-npt','source_job_id':'a'*32};req[field]=value
   with self.assertRaises(ValueError):make_config(req)
 def test_stationary_box_is_not_a_successful_pressure_sampling_screen(self):
  from workbench.engine import density_screen
  spec=make_config({'model_id':'pe-core','protocol_id':'continue-npt','source_job_id':'a'*32})['density_screen'];rows=[]
  for block in range(3):
   for i in range(20):rows.append({'stage':str(block),'ensemble':'NPT','context_time_ps':block*5+i*.25,'density_g_cm3':.9,'volume_nm3':1000.,'temperature_K':310.,'potential_kJ_mol':-10000.})
  self.assertFalse(density_screen(rows,spec)['screen_passed'])
  for i,r in enumerate(rows[-20:]):r['volume_nm3']=1000.+i*.01
  self.assertTrue(density_screen(rows,spec)['screen_passed'])
 def test_pressure_diagnostic_units_and_state_restoration(self):
  import numpy as np
  from openmm import System,NonbondedForce,MonteCarloBarostat,VerletIntegrator,Context,Platform,unit
  system=System();nb=NonbondedForce();nb.setNonbondedMethod(nb.CutoffPeriodic);nb.setCutoffDistance(1)
  for _ in range(2):system.addParticle(12);nb.addParticle(0,.3,0)
  system.addForce(nb);bar=MonteCarloBarostat(1,310,25);system.addForce(bar);system.setDefaultPeriodicBoxVectors((3,0,0),(0,3,0),(0,0,3));integrator=VerletIntegrator(.002);context=Context(system,integrator,Platform.getPlatformByName('Reference'));context.setPositions([[.1,.1,.1],[.5,.1,.1]]);context.setVelocities([[.1,.1,.1],[.1,.1,.1]]);before=context.getState(getPositions=True,getVelocities=True,getEnergy=True);pressure=bar.computeCurrentPressure(context).value_in_unit(unit.bar);after=context.getState(getPositions=True,getVelocities=True)
  # Two independent molecules, zero interaction: instantaneous P=2 KE/(3 V).
  expected=2*before.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)/(3*27)*16.60539067173847
  self.assertAlmostEqual(pressure,expected,places=5)
  np.testing.assert_array_equal(before.getPositions(asNumpy=True)._value,after.getPositions(asNumpy=True)._value);np.testing.assert_array_equal(before.getVelocities(asNumpy=True)._value,after.getVelocities(asNumpy=True)._value)
