"""Durable worker jobs with one compute slot and explicit control messages."""
import json
import os
from pathlib import Path
import re
import subprocess
import threading
import time
import uuid
from .config import ROOT,MODELS,TERMINAL,make_config,availability,runtime,engine_python


def atomic(path,value):
 path=Path(path);tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2));tmp.replace(path)


def alive(pid):
 try:os.kill(int(pid),0);return True
 except (OSError,TypeError,ValueError):return False


class BusyError(ValueError):pass


class Manager:
 def __init__(self,root=ROOT,check_runtime=True):
  self.root=Path(root).resolve();self.jobs=self.root/'runs/workbench';self.jobs.mkdir(parents=True,exist_ok=True);self.lock=threading.Lock();self.check_runtime=check_runtime
  self.runtime=runtime(self.root) if check_runtime else {'ready':True,'python':engine_python(self.root)}
 def directory(self,id):
  if not re.fullmatch(r'[0-9a-f]{32}',id):raise ValueError('Invalid job id')
  path=self.jobs/id
  if not path.is_dir():raise FileNotFoundError(id)
  return path
 def get(self,id):
  path=self.directory(id);cfg=json.loads((path/'config.json').read_text());state=json.loads((path/'state.json').read_text())
  process=json.loads((path/'process.json').read_text()) if (path/'process.json').exists() else {}
  if state['status'] not in TERMINAL and process and not alive(process.get('pid')):
   state['status']='interrupted';state['message']='Worker is no longer running; resume from the last saved checkpoint.';atomic(path/'state.json',state)
  state.update({'id':id,'config':cfg,'process':process})
  if (path/'export-status.json').exists():state['export']=json.loads((path/'export-status.json').read_text())
  log=path/'job.log'
  if log.exists():
   with log.open('rb') as f:f.seek(max(0,log.stat().st_size-16000));state['log_tail']=f.read().decode(errors='replace')
  return state
 def list(self):
  result=[]
  for path in self.jobs.iterdir():
   if path.is_dir() and (path/'config.json').exists():
    try:result.append(self.get(path.name))
    except (ValueError,OSError):continue
  return sorted(result,key=lambda x:x['config']['created_unix'],reverse=True)
 def start(self,request,resume_from=None):
  with self.lock:
   cfg=make_config(request)
   if not self.runtime['ready']:raise ValueError(self.runtime.get('message','Scientific runtime unavailable'))
   if not availability(cfg['model_id'],self.root)['ready']:raise ValueError('Restore the validated raw dataset before running MD')
   if cfg.get('continuation_from_job') and not resume_from:
    parent=self.get(cfg['continuation_from_job'])
    if parent['status'] not in {'completed','needs_review'} or parent['config']['model_id']!=cfg['model_id']:raise ValueError('Continuation requires a completed, compatible source with saved validation')
    for name in ['validation.json','checkpoint.json','final-state.xml','result.json']:
     if not (self.directory(parent['id'])/name).exists():raise ValueError('Continuation source bundle is incomplete')
   if resume_from:
    old=self.get(resume_from)
    if old['status'] not in {'cancelled','interrupted','failed'}:raise ValueError('Only stopped/interrupted jobs can be resumed into a new job')
    cursor=self.directory(resume_from)/'checkpoint.json'
    if not cursor.exists():raise ValueError('No saved checkpoint is available')
    oldcfg=old['config'];cfg={k:v for k,v in oldcfg.items() if k not in {'id','created_unix','resume_from_job'}};cfg['resume_from_job']=resume_from
   slot=self.jobs/'.active'
   if slot.exists():
    owner=json.loads(slot.read_text())
    if alive(owner.get('pid')):raise BusyError('One calculation already owns the slot; finish or cancel it before starting another')
    slot.unlink()
   fd=os.open(slot,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
   id=uuid.uuid4().hex;os.write(fd,json.dumps({'job_id':id,'pid':os.getpid()}).encode());os.close(fd)
   path=self.jobs/id;path.mkdir();cfg.update({'id':id,'created_unix':time.time()});atomic(path/'config.json',cfg);atomic(path/'state.json',{'status':'queued','message':'Worker starting','updated_unix':time.time()})
   env={**os.environ,'OPENMM_CPU_THREADS':'4','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':str(path/'cache/matplotlib'),'XDG_CACHE_HOME':str(path/'cache')}
   try:
    with (path/'job.log').open('w') as log:p=subprocess.Popen([engine_python(self.root),'-m','workbench.worker','--job',str(path)],cwd=self.root,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    atomic(path/'process.json',{'pid':p.pid,'started_unix':time.time()});atomic(slot,{'job_id':id,'pid':p.pid})
   except BaseException:
    slot.unlink(missing_ok=True);raise
   return self.get(id)
 def control(self,id,action):
  if action not in {'pause','resume','cancel'}:raise ValueError('Invalid action')
  item=self.get(id)
  if action=='resume' and item['status'] in {'cancelled','interrupted','failed'}:
   c=item['config'];request={'model_id':c['model_id'],'protocol_id':c['protocol_id']}
   if c['protocol_id']=='continue-npt':request.update(source_job_id=c['continuation_from_job'],blocks=len(c['plan']),block_ps=c['plan'][0]['schedule'][0]['steps']//500)
   elif c['protocol_id']=='nvt-check':request['steps']=c['plan'][0]['schedule'][0]['steps']
   else:request['max_npt_blocks']=len(c['plan'])-4
   return self.start(request,resume_from=id)
  if item['status'] in TERMINAL:raise ValueError('Job has already stopped')
  atomic(self.directory(id)/'control.json',{'action':action,'updated_unix':time.time()})
  return {'id':id,'requested_action':action}
