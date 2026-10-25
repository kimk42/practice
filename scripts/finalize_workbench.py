"""Export a terminal workbench job, then optionally refresh/push its LFS snapshot.

Runs separately from MD so report/archive failures cannot destroy a checkpoint.
Only the dedicated simulation worktree/branch may be pushed by this command.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from workbench.config import TERMINAL,engine_python
from workbench.manager import Manager,atomic,alive

FILES=['report.html','summary.png','pressure.png','observables.csv','config.json','state.json','provenance.json','validation.json','result.json','density-screen.json','trajectory-index.json','resume-provenance.json','continuation-provenance.json','runtime-conditions.json']

def export(job):
 from workbench.report import publish
 publish(job)
 out=ROOT/'artifacts/workbench/preparation';out.mkdir(parents=True,exist_ok=True)
 for name in FILES:
  if (job/name).exists():shutil.copy2(job/name,out/name)
  else:(out/name).unlink(missing_ok=True)
 # Clear only generated solvent diagnostics before considering this new job.
 for name in ['solvent-audit.json','report.html']:(out/'solvent-audit'/name).unlink(missing_ok=True)
 (out/'solvent-audit-error.json').unlink(missing_ok=True)
 # Analyze only stopped, stable trajectories; never interpret preparation as production.
 for name in ['geometry-analysis.json','geometric-observables.csv','geometry-analysis-error.json']:(out/name).unlink(missing_ok=True)
 state=json.loads((job/'state.json').read_text())
 if state['status'] in TERMINAL and (job/'trajectory-index.json').exists() and (job/'reference.npy').exists():
  try:
   from md.research_metrics import analyze
   geometry=analyze(job,out)
   notice='<h2>추가 기하 분석</h2><p>이 작업의 저장된 궤적 '+str(geometry['frames'])+'프레임을 분석했습니다. 각도·접촉·물 근접성은 준비 구간의 기하 지표이며 흡착 세기나 리텐션이 아닙니다.</p>'
   for name in ['geometry-analysis.json','geometric-observables.csv']:
    notice+='<p><a download="'+name+'" href="data:application/octet-stream;base64,'+__import__('base64').b64encode((out/name).read_bytes()).decode()+'">'+name+'</a></p>'
   report=(out/'report.html').read_text();(out/'report.html').write_text(report.replace('</body>',notice+'</body>'))
  except Exception as e:atomic(out/'geometry-analysis-error.json',{'status':'analysis_failed','error':str(e),'raw_MD_outputs_preserved':True})
 if state['status'] in TERMINAL and (job/'final-state.xml').exists():
  try:
   from md.solvent_audit import audit
   audit(job,out/'solvent-audit')
  except Exception as e:atomic(out/'solvent-audit-error.json',{'status':'analysis_failed','error':str(e)})
 atomic(out/'location.json',{'job_id':job.name,'raw_directory':str(job.relative_to(ROOT)),'report_generated_unix':time.time()})
 shutil.copytree(out,job/'analysis',dirs_exist_ok=True)
 if json.loads((job/'config.json').read_text()).get('protocol_id')=='continue-npt':
  dest=ROOT/'artifacts/pe-continuation'
  for name in FILES+['geometry-analysis.json','geometric-observables.csv','geometry-analysis-error.json','solvent-audit-error.json','solvent-audit/solvent-audit.json','solvent-audit/report.html']:
   if not (out/name).exists():(dest/name).unlink(missing_ok=True)
  shutil.copytree(out,dest,dirs_exist_ok=True)
 return out

def snapshot():
 subprocess.run([engine_python(ROOT),'md/collect_results.py'],cwd=ROOT,check=True)
 src=ROOT/'artifacts/collected-results/results-complete.zip';h=hashlib.sha256()
 with src.open('rb') as f:
  while block:=f.read(8*1024**2):h.update(block)
 with zipfile.ZipFile(src) as z:
  manifest=json.loads(z.read('FILES.json'))
  for entry in manifest['files']:
   data=z.read(entry['path'])
   if len(data)!=entry['bytes'] or hashlib.sha256(data).hexdigest()!=entry['sha256']:raise ValueError('Archive member integrity failed: '+entry['path'])
 from scripts.archive_layers import publish_snapshot,digest
 index=publish_snapshot(src,ROOT/'datasets')
 r={'format':index['format'],'index_path':'datasets/results-index.json','index_sha256':digest(ROOT/'datasets/results-index.json'),'archives':index['archives'],'archives_total_bytes':sum(a['bytes'] for a in index['archives']),'snapshot_files':len(index['files']),'local_full_snapshot_sha256':h.hexdigest(),'local_full_snapshot_bytes':src.stat().st_size,'large_system_atoms':577890,'scope':'Complete latest snapshot selected by the SHA-256 index across an immutable baseline and incremental LFS ZIPs. No equilibrium claim.'}
 atomic(ROOT/'docs/result-archive.json',r);print(json.dumps(r),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('--job',required=True);p.add_argument('--wait',action='store_true');p.add_argument('--archive',action='store_true');p.add_argument('--push',action='store_true');a=p.parse_args();m=Manager(check_runtime=False);job=m.directory(a.job)
 if a.push:
  branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
  if branch!='codex/fibrinogen-pe-md':raise SystemExit('Push is restricted to the simulation branch')
 statusfile=job/'export-status.json';atomic(statusfile,{'status':'waiting' if a.wait else 'exporting'})
 try:
  while True:
   s=m.get(a.job)
   # Worker keeps the slot through its own report generation.
   slot=m.jobs/'.active';owner=json.loads(slot.read_text()) if slot.exists() else {};owned=owner.get('job_id')==a.job and alive(owner.get('pid'))
   if slot.exists() and owner.get('job_id')==a.job and not owned:slot.unlink()
   if s['status'] in TERMINAL and not owned:break
   if not a.wait:raise SystemExit('MD is still running; use --wait to export at termination')
   time.sleep(10)
  # Give the independent campaign follower a bounded chance to record this same terminal MD state.
  campaign=ROOT/'artifacts/research-campaign/status.json'
  for _ in range(12):
   if not campaign.exists():break
   research=json.loads(campaign.read_text())
   if research.get('current_PE_job')!=a.job or research.get('current_PE',{}).get('status')==s['status']:break
   time.sleep(5)
  atomic(statusfile,{'status':'exporting','MD_status':s['status']});export(job)
  if json.loads((job/'config.json').read_text())['model_id']=='pe-core' and (ROOT/'artifacts/pe-continuation/geometry-analysis.json').exists():
   try:
    from md.contact_review import publish as publish_contacts
    publish_contacts(job,ROOT/'artifacts/contact-review01')
    shutil.copytree(ROOT/'artifacts/contact-review01',job/'contact-review',dirs_exist_ok=True)
   except Exception as e:atomic(job/'contact-review-error.json',{'error':str(e)})
  if a.archive:snapshot()
  if a.push:
   paths=['artifacts/contact-review01','artifacts/npt-review01','artifacts/workbench','artifacts/pe-continuation','artifacts/solvent-audit01','artifacts/research-campaign','artifacts/collected-results','docs/result-archive.json']
   if a.archive:
    archive_index=json.loads((ROOT/'datasets/results-index.json').read_text());paths+=['datasets/results-index.json']+['datasets/'+a['name'] for a in archive_index['archives']]
   subprocess.run(['git','add','--',*paths],cwd=ROOT,check=True)
   subprocess.run([sys.executable,'scripts/check_repository.py'],cwd=ROOT,check=True)
   staged=subprocess.run(['git','diff','--cached','--quiet'],cwd=ROOT).returncode
   if staged==1:subprocess.run(['git','commit','-m',f'Record executed PE preparation {a.job[:8]} ({s["status"]})'],cwd=ROOT,check=True)
   elif staged!=0:raise RuntimeError('Cannot inspect staged changes')
   if a.archive:subprocess.run(['git','lfs','push','origin','codex/fibrinogen-pe-md'],cwd=ROOT,check=True)
   subprocess.run(['git','push','origin','codex/fibrinogen-pe-md'],cwd=ROOT,check=True)
  atomic(statusfile,{'status':'finished','MD_status':s['status'],'archive_refreshed':a.archive,'pushed':a.push,'updated_unix':time.time()})
 except BaseException as e:
  atomic(statusfile,{'status':'failed','error':str(e),'updated_unix':time.time()});raise

if __name__=='__main__':main()
