"""Reviewable job artifacts; completed execution never implies physical equilibrium."""
import base64
import csv
import html
import json
from pathlib import Path
import struct
from datetime import datetime,timezone
import numpy as np
from openmm import XmlSerializer,unit
from md.dcd_diagnostics import frames
from md.protein_pilot import aligned_rmsd
from .manager import atomic


def verify(engine):
 job=engine.job;by_step={r['context_step']:r for r in engine.rows};counts=[]
 for segment in engine.segments:
  path=job/segment['path'];count=0
  if not path.exists() or path.stat().st_size==0:continue
  for step,xyz,box in frames(path):
   count+=1;context=step+segment['context_step_offset']
   if len(xyz)!=engine.model['atoms'] or not np.isfinite(xyz).all():raise ValueError('Trajectory atom count or coordinates differ')
   row=by_step.get(context)
   if row is None:raise ValueError('Trajectory frame lacks a matching observation')
   if abs(float(np.prod(box))-row['volume_nm3'])>1e-5:raise ValueError('Trajectory and observation cells disagree')
   rmsd=aligned_rmsd(xyz[engine.ca],engine.reference[engine.ca])*10
   if abs(rmsd-row['CA_RMSD_from_pipeline_start_A'])>1e-3:raise ValueError('DCD and observation RMSD disagree')
  counts.append({'path':segment['path'],'frames':count,'context_step_offset':segment['context_step_offset'],'context_time_offset_ps':segment['context_time_offset_ps']})
 state=XmlSerializer.deserialize((job/'final-state.xml').read_text())
 if state.getStepCount()!=engine.rows[-1]['context_step']:raise ValueError('Restart state and observations disagree')
 if abs(state.getTime().value_in_unit(unit.picosecond)-engine.rows[-1]['context_time_ps'])>1e-8:raise ValueError('Restart clock differs')
 result={'status':'saved_output_cross_checks_passed','trajectory_segments':counts,'observations':len(engine.rows),'source_model_atoms':engine.model['atoms'],'native_CA_count':len(engine.ca),'equilibrium_established':False}
 atomic(job/'validation.json',result);return result


def publish(job):
 import os
 os.environ.setdefault('MPLCONFIGDIR',str(Path(job)/'cache/matplotlib'))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 job=Path(job);cfg=json.loads((job/'config.json').read_text());state=json.loads((job/'state.json').read_text());rows=[]
 if (job/'observations.jsonl').exists():
  for line in (job/'observations.jsonl').read_text().splitlines():
   if line.strip():
    try:rows.append(json.loads(line))
    except ValueError:break
 if rows:
  with (job/'observables.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  t=[r['additional_time_ps'] for r in rows];fig,axes=plt.subplots(2,2,figsize=(12,7))
  for ax,key,label in [(axes[0,0],'temperature_K','Kinetic temperature (K)'),(axes[0,1],'density_g_cm3','TOTAL mixture density (g/cm³)'),(axes[1,0],'potential_kJ_mol','Potential + support energy (kJ/mol)'),(axes[1,1],'CA_RMSD_from_pipeline_start_A','Observed CA RMSD\nfrom pipeline start (Å)')]:
   ax.plot(t,[r[key] for r in rows]);ax.set(xlabel='Additional pipeline time (ps)',ylabel=label);ax.grid(alpha=.2)
  axes[0,0].plot(t,[r['target_K'] for r in rows],ls='--',color='orange',label='Target');axes[0,0].legend()
  for ax in axes.flat:
   for a,b in zip(rows,rows[1:]):
    if a['stage']!=b['stage']:ax.axvline(b['additional_time_ps'],ls=':',color='gray',lw=.7)
  fig.suptitle(f'Executed {cfg["model_id"]}: {cfg["protocol_id"]}');fig.tight_layout(rect=(0,.03,1,.94));fig.savefig(job/'summary.png',dpi=140);plt.close(fig)
 pressure_rows=[r for r in rows if r.get('pressure_bar') is not None]
 if pressure_rows:
  fig,ax=plt.subplots(figsize=(10,4));ax.plot([r['additional_time_ps'] for r in pressure_rows],[r['pressure_bar'] for r in pressure_rows],marker='o',ms=3);ax.axhline(1,color='orange',ls='--',label='Applied pressure: 1 bar');ax.set(xlabel='Additional time (ps)',ylabel='Instantaneous molecular-virial pressure (bar)',title='Sparse pressure diagnostics; no equilibrium claim');ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(job/'pressure.png',dpi=140);plt.close(fig)
 image=''
 if (job/'summary.png').exists():image='<img alt="실제 계산 중 기록한 온도, 전체 밀도, 에너지와 구조" src="data:image/png;base64,'+base64.b64encode((job/'summary.png').read_bytes()).decode()+'">'
 if pressure_rows:image+='<img alt="실제로 계산한 순간 압력, 수렴 검증 아님" src="data:image/png;base64,'+base64.b64encode((job/'pressure.png').read_bytes()).decode()+'">'
 result=json.loads((job/'result.json').read_text()) if (job/'result.json').exists() else state
 stamp=datetime.fromtimestamp(state.get('updated_unix',cfg['created_unix']),timezone.utc).isoformat()
 status=html.escape(state['status']);body=f'<h1>실제 PE 계산 작업</h1><p>작업 {html.escape(cfg["id"])} · {html.escape(cfg["model_id"])} · {status}</p><p class="notice">기록된 실행과 검증 결과입니다. 실행 완료 또는 기술적 밀도 추세 검사 통과를 계면 평형·생산 계산 준비·결합 자유에너지 검증으로 해석하지 않습니다.</p>'
 if rows:body+=f'<p>실제로 추가 실행한 시간 {rows[-1]["additional_time_ps"]:.4f} ps · 마지막 실제 온도 {rows[-1]["temperature_K"]:.2f} K · 전체 밀도 {rows[-1]["density_g_cm3"]:.4f} g/cm³</p>'
 body+='<p>실행 상태 기준 시각 (UTC): '+stamp+'</p>'
 body+=image+'<h2>조건과 해석</h2><p>NVT의 평평한 밀도는 상자 고정의 결과입니다. NPT에서는 이 유한 PE 패치에 등방성 barostat를 사용하며 위치 억제를 모두 해제합니다. 이 방식을 무한 표면·막 모델에 그대로 적용하지 않습니다. 서로 다른 위치 억제·ensemble 구간의 총 에너지를 결합 세기로 비교하지 않습니다.</p><p>구간별 온도·전체 밀도·에너지 기울기는 상관된 짧은 표본의 기술 지표입니다. 결손 말단·당쇄 생략·중화 이온만 포함한 조건과 유한 올리고머 표면의 한계가 적용됩니다.</p>'
 body+='<h2>실제 결과</h2><pre>'+html.escape(json.dumps(result,ensure_ascii=False,indent=2))+'</pre><details><summary>고정 실행 설정</summary><pre>'+html.escape(json.dumps(cfg,ensure_ascii=False,indent=2))+'</pre></details>'
 downloads=[]
 for name in ['observables.csv','config.json','provenance.json','validation.json','density-screen.json','resume-provenance.json','continuation-provenance.json','trajectory-index.json']:
  p=job/name
  if p.exists():downloads.append('<a download="'+name+'" href="data:application/octet-stream;base64,'+base64.b64encode(p.read_bytes()).decode()+'">'+name+'</a>')
 body+='<h2>주요 데이터 내려받기</h2><p>'+' · '.join(downloads)+'</p>'
 css='body{max-width:1100px;margin:30px auto;padding:20px;font:16px/1.7 system-ui;color:#183247}img{width:100%}.notice{background:#fff0cf;padding:16px}pre{white-space:pre-wrap;overflow-wrap:anywhere}'
 (job/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PE 계산 작업 결과</title><style>'+css+'</style><body>'+body+'</body></html>')
