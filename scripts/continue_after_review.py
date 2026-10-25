"""Diagnose the reviewed PE run, execute water controls, conditionally continue NPT."""
import argparse
import base64
import csv
from datetime import datetime,timezone
import html
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from workbench.manager import Manager,atomic
from md.water_control import run


def review(parent):
 rows=list(map(json.loads,(parent/'observations.jsonl').read_text().splitlines()));blocks=[]
 for name in dict.fromkeys(r['stage'] for r in rows):
  samples=[r for r in rows if r['stage']==name];rho=[r['density_g_cm3'] for r in samples];pressure=[r['pressure_bar'] for r in samples if r.get('pressure_bar') is not None];changes=int(np.count_nonzero(np.abs(np.diff([r['volume_nm3'] for r in samples]))>1e-6));blocks.append({'stage':name,'observed_volume_changes':changes,'observation_intervals':len(samples)-1,'density_first_last_g_cm3':[rho[0],rho[-1]],'pressure_sample_count':len(pressure),'pressure_mean_bar':float(np.mean(pressure)),'pressure_range_bar':[min(pressure),max(pressure)]})
 return {'source_job':parent.name,'blocks':blocks,'source_output_validation':json.loads((parent/'validation.json').read_text()),'source_hard_overlap_checks':json.loads((parent/'result.json').read_text())['geometry_after'],'conclusion':'Barostat volume changes are observed and increasing; density continues to rise and sampled pressure is predominantly negative. This is consistent with incomplete volume relaxation, not a demonstrated single root cause. Water-only implementation controls precede further PE NPT.','pressure_caveat':'Sparse correlated instantaneous samples, no converged average or confidence interval.','equilibrium_established':False}


def publish(out,state):
 atomic(out/'status.json',state);body='<h1>PE 압력·밀도 검토와 물 단독 대조</h1><p>상태: '+html.escape(state['status'])+'</p><p>큰 PE 모델의 준비 검사 실패 원인을 검토합니다. 물 단독 계산은 단백질·PE 계산 시간에 더하지 않습니다. 대조 계산 통과를 폴리머 표면이나 전체 계의 평형 검증으로 해석하지 않습니다.</p><pre>'+html.escape(json.dumps(state,ensure_ascii=False,indent=2))+'</pre>'
 child=out/'controls/report.html'
 if child.exists():body+='<hr>'+child.read_text().split('<body>',1)[-1].rsplit('</body>',1)[0]
 (out/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PE NPT 검토와 물 단독 대조</title><style>body{max-width:1100px;margin:auto;padding:24px;font:16px/1.7 system-ui;color:#183247}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style><body>'+body+'</body></html>')


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--source-job',required=True);a=parser.parse_args();manager=Manager();parent=manager.directory(a.source_job);out=ROOT/'artifacts/npt-review01';out.mkdir(exist_ok=True);raw=ROOT/'runs/npt-review01';raw.mkdir(exist_ok=False);state={'status':'running_water_controls','created_utc':datetime.now(timezone.utc).isoformat(),'review':review(parent),'planned_water_controls':{'replicas':2,'water_molecules_per_replica':512,'NVT_ps':5,'NPT_ps':100,'starting_density_g_cm3':.86},'next_PE_plan':'Only if both technical controls pass: checkpoint-continuous 15 ps PE NPT, 310 K, 1 bar, same timestep, no velocity reinitialization.'};publish(out,state)
 try:
  result=run(parent.name,raw/'water-controls')
  from md.plot_water_control import publish as publish_controls
  publish_controls(raw/'water-controls',out/'controls')
  state['water_controls']=result;atomic(out/'water-controls.json',result);publish(out,state)
  if result['status']!='technical_controls_passed':state['status']='needs_review_no_PE_launched';publish(out,state);return
  job=manager.start({'model_id':'pe-core','protocol_id':'continue-npt','source_job_id':parent.name,'blocks':3});state.update(status='water_controls_passed_PE_continuation_started',PE_job_id=job['id']);atomic(ROOT/'runs/pe-continuation-active.json',{'job_id':job['id'],'source_job_id':parent.name});publish(out,state)
  jobdir=manager.directory(job['id']);env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':str(jobdir/'cache/matplotlib'),'XDG_CACHE_HOME':str(jobdir/'cache')}
  for label,argv in [('export',[str(ROOT/'.venv/bin/python'),'scripts/finalize_workbench.py','--job',job['id'],'--wait','--archive','--push']),('research-follow',[str(ROOT/'.venv/bin/python'),'-m','workbench.research','--pe-job',job['id'],'--amberhome','/workspace/tools/ambertools','--follow-existing'])]:
   with (jobdir/(label+'.log')).open('w') as log:p=subprocess.Popen(argv,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
   atomic(jobdir/(label+'-process.json'),{'pid':p.pid,'started_unix':time.time(),'argv':argv})
 except BaseException as e:
  state.update(status='needs_review',error=str(e));publish(out,state);raise

if __name__=='__main__':main()
