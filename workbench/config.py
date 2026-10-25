"""Versioned, allowlisted run configurations. No arbitrary browser commands/paths."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
TERMINAL={'completed','needs_review','failed','cancelled','interrupted'}
MODELS={
 'pe-fragment':{'label':'PE + γ절편 + 물·중화이온','atoms':43288,'native_residues':252,'prmtop':'data/pe-protein-model04/complex.prmtop','manifest':'data/pe-protein-model04/manifest.json','source':'runs/pe-staged-preparation01/03-restrained-interface-NVT-check','source_parent':'runs/pe-staged-preparation01','temperature_K':310,'friction_per_ps':1,'dt_fs':2,'seed':42},
 'pe-core':{'label':'PE + 여섯 사슬 관측 코어 + 물·중화이온','atoms':577890,'native_residues':1913,'prmtop':'data/fibrinogen-core-pe-solvated01/solvated.prmtop','manifest':'data/fibrinogen-core-pe-solvated01/manifest.json','source':'runs/core-pe-100K-hold01/output','temperature_K':100,'friction_per_ps':5,'dt_fs':1,'seed':45},
}

def engine_python(root=ROOT):
 explicit=os.environ.get('FIBRINOGEN_PYTHON')
 if explicit:return str(Path(explicit).expanduser().absolute())
 for p in [Path(root)/'.venv/bin/python',Path(root)/'.venv/Scripts/python.exe']:
  if p.exists():return str(p)
 return sys.executable


def runtime(root=ROOT):
 try:
  p=subprocess.run([engine_python(root),'-c','import openmm,numpy,scipy,matplotlib; print(openmm.__version__)'],capture_output=True,text=True,timeout=30)
  return {'ready':p.returncode==0,'python':engine_python(root),'openmm':p.stdout.strip() if p.returncode==0 else None,'message':None if p.returncode==0 else '계산용 Python에 requirements.txt 패키지를 설치하세요.'}
 except (OSError,subprocess.TimeoutExpired):return {'ready':False,'python':engine_python(root),'message':'계산용 Python을 실행할 수 없습니다.'}


def availability(model_id,root=ROOT):
 m=MODELS[model_id];root=Path(root);source=root/m['source']
 names=[m['prmtop'],m['manifest'],m['source']+'/state.xml',m['source']+'/system.xml',m['source']+'/checkpoint.chk']
 if model_id=='pe-core':names+=['data/fibrinogen-core03/manifest.json',m['source']+'/validation.json',m['source']+'/result.json',m['source']+'/integrator-final.xml']
 else:names += [m['source_parent']+'/result.json',m['source_parent']+'/validation.json']
 missing=[name for name in names if not (root/name).is_file()]
 return {'ready':not missing,'missing':missing,'message':'원본 데이터 복원이 필요합니다.' if missing else '검증된 시작 번들 사용 가능'}


def make_config(request):
 if not isinstance(request,dict):raise ValueError('JSON object required')
 allowed={'model_id','protocol_id','steps','max_npt_blocks','source_job_id','blocks','block_ps'}
 if set(request)-allowed:raise ValueError('Unsupported configuration fields: '+', '.join(sorted(set(request)-allowed)))
 model=request.get('model_id');protocol=request.get('protocol_id')
 if model not in MODELS:raise ValueError('원자 모델이 준비된 PE만 실행할 수 있습니다. PTFE/PVDF는 아직 실행할 수 없습니다.')
 cfg={'schema_version':1,'model_id':model,'protocol_id':protocol,'threads':4,'observation_interval_steps':25,'trajectory_interval_steps':250,'checkpoint_interval_steps':100,'seed':MODELS[model]['seed'],'equilibrium_established':False,'production_ready':False}
 if protocol!='continue-npt' and {'source_job_id','blocks','block_ps'}&set(request):raise ValueError('Continuation fields require continue-npt')
 if protocol=='nvt-check' and model=='pe-fragment':
  if 'max_npt_blocks' in request:raise ValueError('NPT options are not part of this protocol')
  steps=request.get('steps',100)
  if type(steps)is not int or not 25<=steps<=1000 or steps%25:raise ValueError('steps must be a multiple of 25 between 25 and 1000')
  cfg['trajectory_interval_steps']=25;cfg['checkpoint_interval_steps']=25
  cfg['plan']=[{'name':'prepared-fragment-NVT','ensemble':'NVT','dt_fs':2,'schedule':[{'target_K':310,'steps':steps,'restraint_scale':1}]}]
 elif protocol=='prepare-310' and model=='pe-core':
  if 'steps' in request:raise ValueError('Preparation stage lengths are defined by the reviewed protocol')
  blocks=request.get('max_npt_blocks',6)
  if type(blocks)is not int or not 3<=blocks<=6:raise ValueError('max_npt_blocks must be 3..6')
  cfg['plan']=[
   {'name':'warm-to-200K','ensemble':'NVT','dt_fs':1,'schedule':[{'target_K':t,'steps':200,'restraint_scale':1} for t in [120,140,160,180,200]]},
   {'name':'warm-to-310K','ensemble':'NVT','dt_fs':1,'schedule':[{'target_K':t,'steps':100,'restraint_scale':1} for t in range(210,311,10)]},
   {'name':'310K-supported-check','ensemble':'NVT','dt_fs':1,'schedule':[{'target_K':310,'steps':1000,'restraint_scale':1}]},
   {'name':'release-position-support','ensemble':'NVT','dt_fs':1,'schedule':[{'target_K':310,'steps':250,'restraint_scale':.1},{'target_K':310,'steps':250,'restraint_scale':.01},{'target_K':310,'steps':500,'restraint_scale':0}]},
  ]
  cfg['plan'] += [{'name':f'NPT-density-block-{i+1}','ensemble':'NPT','dt_fs':2,'pressure_bar':1,'schedule':[{'target_K':310,'steps':1250,'restraint_scale':0}]} for i in range(blocks)]
  cfg['density_screen']={'minimum_blocks':3,'block_time_ps':2.5,'maximum_adjacent_mean_density_change_percent':.5,'maximum_last_block_density_drift_percent':.5,'maximum_temperature_mean_deviation_K':5,'maximum_relative_potential_slope_percent_per_ps':.5,'interpretation':'Descriptive preparation screen only; correlated short samples do not establish bulk or interface equilibrium.'}
 elif protocol=='continue-npt' and model=='pe-core':
  if {'steps','max_npt_blocks'}&set(request):raise ValueError('Use bounded continuation blocks')
  source=request.get('source_job_id','');blocks=request.get('blocks',3);block_ps=request.get('block_ps',5)
  if not isinstance(source,str) or not re.fullmatch(r'[0-9a-f]{32}',source):raise ValueError('A validated source job id is required')
  if type(blocks)is not int or not 1<=blocks<=3:raise ValueError('Continuation blocks must be 1..3')
  if type(block_ps)is not int or block_ps not in (5,10):raise ValueError('block_ps must be 5 or 10')
  cfg.update(continuation_from_job=source,pressure_interval_steps=250)
  cfg['plan']=[{'name':f'continued-NPT-block-{i+1}','ensemble':'NPT','dt_fs':2,'pressure_bar':1,'schedule':[{'target_K':310,'steps':block_ps*500,'restraint_scale':0}]} for i in range(blocks)]
  cfg['density_screen']={'minimum_blocks':3,'block_time_ps':block_ps,'maximum_adjacent_mean_density_change_percent':.5,'maximum_last_block_density_drift_percent':.5,'maximum_temperature_mean_deviation_K':5,'maximum_relative_potential_slope_percent_per_ps':.5,'minimum_last_block_observed_volume_changes':5,'interpretation':'Short descriptive continuation screen with barostat-activity check; not proof of equilibrium or physical density accuracy.'}
 else:raise ValueError('Protocol/model combination is not supported')
 cfg['maximum_planned_steps']=sum(x['steps'] for stage in cfg['plan'] for x in stage['schedule'])
 cfg['maximum_planned_time_ps']=sum(sum(x['steps'] for x in s['schedule'])*s['dt_fs']/1000 for s in cfg['plan'])
 cfg['stage_notes']=['NVT keeps the original box and anchors; support is released before pressure adjustment.','Isotropic NPT is used only for this finite oligomer patch; all Cartesian support force constants must be zero.','No surface affinity, binding free energy, material ranking or interface-equilibrium conclusion.']
 return cfg


def fingerprint(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  while b:=f.read(8*1024**2):h.update(b)
 return h.hexdigest()
