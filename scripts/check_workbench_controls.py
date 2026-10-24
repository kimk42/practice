import json,os,time,threading
from pathlib import Path
from urllib.request import Request,urlopen
from http.server import ThreadingHTTPServer
from workbench.config import ROOT
from workbench.manager import Manager
from workbench.server import application
root=ROOT/'runs/control-integration';root.mkdir(exist_ok=True);(root/'runs').mkdir(exist_ok=True)
for name,target in [('data',ROOT/'data'),('runs/pe-staged-preparation01',ROOT/'runs/pe-staged-preparation01')]:
 p=root/name
 if not p.exists():p.symlink_to(target,target_is_directory=True)
os.environ['FIBRINOGEN_PYTHON']=str(ROOT/'.venv/bin/python');os.environ['PYTHONPATH']=str(ROOT)
m=Manager(root);server=ThreadingHTTPServer(('127.0.0.1',0),application(m));threading.Thread(target=server.serve_forever,daemon=True).start();base=f'http://127.0.0.1:{server.server_port}'
def get(path):
 with urlopen(base+path) as r:return json.load(r)
token=get('/api/capabilities')['token']
def post(path,body):
 with urlopen(Request(base+path,data=json.dumps(body).encode(),headers={'Content-Type':'application/json','X-Workbench-Token':token})) as r:return json.load(r)
def until(id,predicate,timeout=180):
 end=time.monotonic()+timeout
 while time.monotonic()<end:
  s=get('/api/jobs/'+id)
  if predicate(s):return s
  if s['status'] in {'failed','needs_review'}:raise RuntimeError(s)
  time.sleep(.5)
 raise TimeoutError(id)
try:
 j=post('/api/jobs',{'model_id':'pe-fragment','protocol_id':'nvt-check','steps':1000});id=j['id'];until(id,lambda s:s.get('last',{}).get('stage_step',0)>=25)
 post(f'/api/jobs/{id}/pause',{});s=until(id,lambda s:s['status']=='paused');step=s['last']['context_step'];time.sleep(1);assert get('/api/jobs/'+id)['last']['context_step']==step
 post(f'/api/jobs/{id}/resume',{});until(id,lambda s:s.get('last',{}).get('context_step',0)>step)
 post(f'/api/jobs/{id}/cancel',{});s=until(id,lambda s:s['status']=='cancelled');cursor=json.loads((m.directory(id)/'checkpoint.json').read_text());
 while (m.jobs/'.active').exists():time.sleep(.25)
 new=post(f'/api/jobs/{id}/resume',{});newid=new['id'];done=until(newid,lambda s:s['status']=='completed',240)
 proof=json.loads((m.directory(newid)/'resume-provenance.json').read_text());assert proof['valid_parent_context_step']==cursor['context_step'];assert done['result']['context_step']==2000;assert done['result']['additional_time_ps']>1.99
 result={'status':'real_MD_API_pause_resume_cancel_checkpoint_restart_passed','parent_job':id,'resumed_job':newid,'paused_context_step':step,'restart_context_step':cursor['context_step'],'final_context_step':done['result']['context_step'],'additional_time_ps':done['result']['additional_time_ps'],'final_validation':json.loads((m.directory(newid)/'validation.json').read_text())}
 (ROOT/'artifacts/workbench').mkdir(exist_ok=True);(ROOT/'artifacts/workbench/control-integration.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
finally:server.shutdown();server.server_close()
