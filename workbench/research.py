"""Evidence-gated research preparation with experimental inputs explicitly deferred.

Creates real candidate molecules and rigid protein placements. A proposal/queue is
never reported as executed MD, validated crystals, PMF data or a grant submission.
"""
import argparse
import base64
from datetime import datetime,timezone
import hashlib
import html
import json
import os
from pathlib import Path
import time
import numpy as np
from openmm import unit
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,PDBxFile,element
from .config import ROOT
from .manager import atomic,Manager
from md.prepare_fluoropolymers import parameterize,molecule,conformer
from md.build_pe_complex import principal_orientation

MATERIALS=['PE','PTFE','PVDF-alpha','PVDF-beta','PVDF-gamma']
ORIENTATIONS=[('side-0',0,0),('side-90',90,0),('side-180',180,0),('side-270',270,0),('end-plus',0,90),('end-minus',0,-90)]
SEEDS=[101,202,303]


def rotation(roll_deg,pitch_deg):
 a,b=np.deg2rad([roll_deg,pitch_deg]);rx=np.array([[1,0,0],[0,np.cos(a),-np.sin(a)],[0,np.sin(a),np.cos(a)]]);ry=np.array([[np.cos(b),0,np.sin(b)],[0,1,0],[-np.sin(b),0,np.cos(b)]]);return ry@rx


def place(xyz,heavy,roll,pitch,gap_nm=.5):
 r=rotation(roll,pitch);p=np.asarray(xyz)@r.T;p[:,2]+=gap_nm-p[heavy,2].min()
 if not np.allclose(r@r.T,np.eye(3),atol=1e-12) or np.linalg.det(r)<.999999:raise ValueError('Invalid proper rigid rotation')
 return p,r


def placements(root,out):
 core=root/'data/fibrinogen-core03';manifest=json.loads((core/'manifest.json').read_text());p=core/'local-relaxation-restrained/core-relaxed.inpcrd';expected=json.loads((root/'data/fibrinogen-core-pe02/manifest.json').read_text())['source_relaxed_coordinates_sha256']
 if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise ValueError('Reviewed core coordinates differ')
 prm=AmberPrmtopFile(str(core/'core.prmtop'));x=AmberInpcrdFile(str(p));xyz=principal_orientation(x.positions.value_in_unit(unit.nanometer));atoms=list(prm.topology.atoms());heavy=np.array([a.index for a in atoms if a.element!=element.hydrogen]);out.mkdir(parents=True,exist_ok=False);result=[];bonds=np.array([(a.index,b.index) for a,b in prm.topology.bonds()]);d0=np.linalg.norm(xyz[bonds[:,0]]-xyz[bonds[:,1]],axis=1)
 for name,roll,pitch in ORIENTATIONS:
  pos,r=place(xyz,heavy,roll,pitch);error=float(np.max(np.abs(np.linalg.norm(pos[bonds[:,0]]-pos[bonds[:,1]],axis=1)-d0)))
  if error>1e-10:raise ValueError('Rigid placement changed a bond')
  file=out/(name+'.cif')
  with file.open('w') as f:PDBxFile.writeFile(prm.topology,pos*unit.nanometer,f)
  result.append({'id':name,'roll_deg':roll,'pitch_deg':pitch,'rotation_matrix':r.tolist(),'atoms':len(pos),'minimum_solute_heavy_z_nm':float(pos[heavy,2].min()),'maximum_bond_length_change_nm':error,'extent_nm':np.ptp(pos,axis=0).tolist(),'sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'coordinates':str(file.relative_to(root))})
 metadata={'status':'dry_rigid_coordinate_candidates_only','source_coordinate_sha256':expected,'native_residues':manifest['observed_residues'],'reference_plane':'mathematical z=0, no polymer, water or force calculation in these files','Ca_ions':'four source Ca ions co-rotated with capped observed core; this is a coordinate transform, not Ca binding validation','placements':result,'solvated':False,'MD_steps':0,'equilibrium_established':False};atomic(out/'manifest.json',metadata);return metadata


def budget(pe):
 last=pe.get('last') or {};t=last.get('additional_time_ps',0)
 if 'current_process_start_context_time_ps' in pe:t=max(0,last.get('context_time_ps',pe['current_process_start_context_time_ps'])-pe['current_process_start_context_time_ps'])
 wall=pe.get('elapsed_wall_seconds',0);rate=t/1000/max(wall,1)*86400
 # These are explicit planning scenarios, not sample lengths that establish convergence.
 return {'observed_current_process_simulated_ps':t,'observed_preparation_ns_per_wall_day':rate,'source':'Observed current preparation elapsed wall time, including reporting/control pauses. Not a production benchmark.','GPU_throughput':'not measured; no GPU assumed','pilot_matrix':{'materials':5,'orientations':6,'velocity_seeds':3,'replicas':90,'suggested_first_check_ps_per_replica':2,'total_first_check_ns':.18},'production_scenario_only':{'ns_per_replica':2,'replicas':90,'total_ns':180,'estimated_current_CPU_days_at_observed_preparation_rate':180/rate if rate else None,'equilibrium_guaranteed':False},'PMF_scenario_only':{'conditional_orientations':1,'independent_replicas':3,'relative_displacement_windows':41,'suggested_sampling_ns_per_window':10,'total_sampling_ns':1230,'estimated_current_CPU_days_at_observed_preparation_rate':1230/rate if rate else None,'includes_window_preparation':False,'equilibrium_guaranteed':False},'interpretation':'Extrapolation only. Ensemble, timestep, barostat and hardware alter throughput. Actual planning must use a production benchmark; massive scenarios are not automatically launched on this CPU.'}


def registry():
 return {m:{'chemical_family':'PVDF' if m.startswith('PVDF') else m,'requested_phase':m.split('-')[1] if m.startswith('PVDF-') else None,'model_level':'existing finite exploratory oligomer patch' if m=='PE' else 'finite starting-chain candidate only','crystallographic_face':None,'unit_cell_source':None,'bulk_density_validated':False,'water_interface_validated':False,'phase_or_domain_validated':False,'fixed_charge_model_limit':'No induced polarization or ferroelectric switching is represented.','surface_MD_enabled':m=='PE','comparative_adsorption_production_allowed':False,'missing_evidence':['documented phase/packing and exposed face','polymer force-field torsional/structural validation','water-interface validation with uncertainty','surface stability at comparison conditions','finite-chain end-group sensitivity versus bulk surface']} for m in MATERIALS}


def pmf_plan():
 return {'status':'design_only_not_executable','MD_steps':0,'free_energy_values':None,'CV':'Protein+cap mass COM separation from a documented substrate reference group along a validated surface normal; use a differentiable, unwrapped periodic-coordinate implementation.','displacement_grid_nm':[round(i*.1,2) for i in range(41)],'grid_note':'Relative displacements from a prepared bound orientation; these are not absolute surface distances or a claim of a bound state. Window positions must be rebuilt from validated geometry.','initial_spring_k_kJ_mol_nm2':1000,'required_before_execution':['verified bound-state/orientation ensemble, not an arbitrary initial close placement','surface model/force field validated and equilibrated','solvent-only protein and protein-free surface controls','enlarged solvent cell: full displacement range cannot reuse the current cell without image/separation checks','defined protein/substrate atom groups and a continuous CV, including PBC unit tests','consistent lateral/orientation restraints with their entropic effects quantified','prepared and equilibrated window starts with safe forces; budget authorization is not a convergence test'],'required_analysis':['overlap matrix of neighboring CV histograms','autocorrelation-aware effective samples; independent runs and block stability','WHAM/MBAR with uncertainty and coordinate/Jacobian conventions stated','desorbed-region stability, no patch-edge escape and no periodic-image artifacts','orientation/lateral restraint and standard-state/area corrections before any adsorption free-energy claim'],'allowed_interpretation':'At first, an orientation-conditional PMF for a defined finite model. Neither a raw restraint energy nor a conditional PMF alone is experimental retention or a standard binding free energy.'}


class Campaign:
 def __init__(self,path,pe_job,amber,root=ROOT):
  self.root=Path(root);self.path=Path(path);self.pe_job=pe_job;self.amber=Path(amber);self.out=self.root/'artifacts/research-campaign';self.out.mkdir(exist_ok=True);self.data=self.root/'data/fluoropolymer-candidates01';self.models=registry();self.manager=Manager(self.root,check_runtime=False)
  self.state={'schema_version':1,'campaign':'research-campaign01','created_utc':datetime.now(timezone.utc).isoformat(),'status':'preparing','experimental_inputs':'deferred_by_user','current_PE_job':pe_job,'defaults':{'temperature_K':310,'pressure_bar':1,'water_model':'TIP3P','protein_forcefield':'ff14SB','salt_condition':'neutralization only, not 150 mM','pH':None,'protonation':'source fixed ff14SB states; no controlled-pH or pKa validation','radioiodination':'not represented','protein_scope':'3GHG observed capped six-chain core; missing tails/glycans omitted'},'molecular_candidates':{},'surface_models':self.models,'stages':{},'execution_policy':'Continue independent preparation. Never fabricate crystals, equilibrated replicas, PMF values, experimental agreement or grant submission. Stop dependent production at failed scientific/resource gates.'}
 def publish(self,write_private=True):
  try:pe=self.manager.get(self.pe_job)
  except (OSError,ValueError):pe={}
  parent=pe.get('config',{}).get('resume_from_job')
  if parent:
   cursor=self.root/'runs/workbench'/parent/'checkpoint.json'
   if cursor.exists():pe['current_process_start_context_time_ps']=json.loads(cursor.read_text())['context_time_ps']
  literature=self.out/'literature-search.json'
  if literature.exists():
   evidence=json.loads(literature.read_text());self.state['reference_lookup']={'status':'unavailable' if all('error' in x for x in evidence) else 'metadata_only','record':'artifacts/research-campaign/literature-search.json','note':'Network requests returned 403; no paper content, lattice or validated force-field parameters were retrieved.' if all('error' in x for x in evidence) else 'Metadata lookup alone does not validate a model.'}
  self.state['independent_molecular_checks']={}
  for family in ['PTFE','PVDF']:
   p=self.data/family/'independent-validation.json'
   if p.exists():self.state['independent_molecular_checks'][family]=json.loads(p.read_text())
  self.state['current_PE']={k:pe.get(k) for k in ['status','last','elapsed_wall_seconds','message']};self.state['resource_budget']=budget(pe);self.state['updated_utc']=datetime.now(timezone.utc).isoformat()
  if write_private:atomic(self.path/'state.json',self.state)
  atomic(self.out/'status.json',self.state)
  rows=''.join('<tr><td>'+str(n)+'</td><td>'+html.escape(x['status'])+'</td><td>'+html.escape(x['detail'])+'</td></tr>' for n,x in self.state['stages'].items());image=''
  if (self.out/'candidates.png').exists():image='<img alt="실제로 만든 유한 사슬 후보의 초기 좌표, 결정 표면 아님" src="data:image/png;base64,'+base64.b64encode((self.out/'candidates.png').read_bytes()).decode()+'">'
  molecule_downloads=[]
  for p in sorted(self.data.glob('*/*.sdf')):
   molecule_downloads.append('<a download="'+p.parent.name+'-'+p.name+'" href="data:application/octet-stream;base64,'+base64.b64encode(p.read_bytes()).decode()+'">'+p.parent.name+'/'+p.name+'</a>')
  body='<h1>실험 입력을 미룬 후속 연구 준비</h1><p>실제 생성·실행과 아직 계획인 작업을 구분한 기록입니다. 3단계의 실험 데이터 연결은 사용자 요청으로 미뤘습니다.</p><p class="notice">PTFE/PVDF 유한 사슬 후보는 검증된 결정 표면이 아닙니다. 단백질의 dry 배향 파일은 용매화·MD 결과가 아닙니다. 본 계산·PMF의 숫자를 아직 얻지 않았습니다.</p><table><tr><th>단계</th><th>상태</th><th>근거·다음 조건</th></tr>'+rows+'</table>'+image+'<h2>실제 기록과 예산</h2><pre>'+html.escape(json.dumps(self.state,ensure_ascii=False,indent=2))+'</pre><h2>자유에너지 계산 설계</h2><pre>'+html.escape(json.dumps(pmf_plan(),ensure_ascii=False,indent=2))+'</pre>'
  body+='<h2>유한 후보 사슬 파일 (결정 표면 아님)</h2><p>'+' · '.join(molecule_downloads)+'</p>'
  css='body{max-width:1100px;margin:auto;padding:28px;font:16px/1.7 system-ui;color:#183247}table{border-collapse:collapse;width:100%}td,th{border:1px solid #bbb;padding:10px;text-align:left}.notice{background:#fff0cf;padding:16px}pre{white-space:pre-wrap;overflow-wrap:anywhere}img{width:100%}'
  (self.out/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>후속 연구 준비와 실행 상태</title><style>'+css+'</style><body>'+body+'</body></html>')
  # Expose only small, clearly labeled candidate coordinates/status; raw topology/logs stay in the LFS archive.
  import shutil
  for family in ['PTFE','PVDF']:
   source=self.data/family
   if source.exists():
    dest=self.out/'candidate-molecules'/family;dest.mkdir(parents=True,exist_ok=True)
    for p in [*source.glob('*.sdf'),source/'status.json',source/'independent-validation.json']:
     if p.is_file():shutil.copy2(p,dest/p.name)
  self.brief()
 def brief(self):
  text='''# 과제 논의를 위한 예비 자료 초안 — 제출·발송하지 않음

질문: 폴리머의 표면 화학과 PVDF 결정상·배향·분극 조건이 관측 피브리노겐 코어의 배치와 물을 통한 상호작용에 어떤 차이를 만드는가?

현재 증거: PE의 실제 용매화·중화 입력, 짧은 준비 MD, 체크포인트 재시작이 가능한 실행 도구가 있다. 새 PTFE/PVDF 파일은 분자 후보 또는 힘장 기술 검사이며 결정 표면·흡착 비교 데이터가 아니다. 큰 PE 계산의 최신 종료 상태는 함께 제공한 status.json과 준비 보고서가 기준이다.

후속 목표: (1) 표면 구조와 물 계면을 먼저 검증한다. (2) 여러 배향·시드로 접촉·배향·구조·계면 수화를 비교한다. (3) 선택한 모델에서 방향 조건부 탈착 PMF를 계산하고 수렴·불확실성을 제시한다. 나중에 사용자의 I-125 흡착·세척 후 리텐션 자료와 비교할 가설을 정한다.

필요 자원: 현재 577,890원자 CPU 준비 계산의 관측 속도로 긴 반복·PMF를 바로 끝낼 수 없다. resource_budget의 시나리오는 실행 결과가 아닌 자원 산정이다. GPU 실측 후 복제 수·시간·표면 수를 조정한다. GPU 확보 전에는 입력 준비·분석 코드·기술 검사를 진행한다.

해석 한계: 관측 코어는 긴 αC·기타 결손 말단과 당쇄를 생략한다. 유한 PE 패치는 검증된 벌크 결정이 아니다. 중화 이온만 포함하며 실험 조건은 미입력이다. 방사성 표지를 모델링하지 않았다. 접촉 수·총 에너지·짧은 체류시간을 흡착 자유에너지나 실험 리텐션으로 환산하지 않는다. 현재 재료 순위를 주장하지 않는다.

제안: 현재 결과를 ‘도구와 실행 가능성’의 예비 근거로 보여주고, 표면 모델 검증과 반복 계산을 위한 작은 과제를 논의한다. 이 초안은 PI에게 보내거나 제출한 문서가 아니다.
'''
  (self.out/'proposal-brief.md').write_text(text)
 def plot(self):
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig,axes=plt.subplots(2,2,figsize=(12,7),subplot_kw={'projection':'3d'})
  for ax,(family,name) in zip(axes.flat,[('PTFE','ptfe'),('PVDF','pvdf-alpha'),('PVDF','pvdf-beta'),('PVDF','pvdf-gamma')]):
   mol=conformer(molecule(family),name);xyz=mol.GetConformer().GetPositions();xyz-=xyz.mean(axis=0)
   for bond in mol.GetBonds():
    ij=[bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()];ax.plot(*xyz[ij].T,color='#7d8b97',lw=.8)
   for symbol,color in [('C','#315466'),('F','#49aa85'),('H','#c2c8cd')]:
    ids=[a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol()==symbol]
    if ids:ax.scatter(*xyz[ids].T,s=10 if symbol=='H' else 22,c=color)
   ax.set_title(name+' finite starting chain');ax.set(xlabel='Å',ylabel='Å',zlabel='Å');ax.set_box_aspect(np.maximum(np.ptp(xyz,axis=0),2))
  fig.suptitle('Generated coordinates — NO crystal/surface or adsorption validation');fig.tight_layout();fig.savefig(self.out/'candidates.png',dpi=130);plt.close(fig)
 def run(self):
  self.state['stages']={'1':{'status':'running_PE_preparation','detail':'Read the actual existing job; no new equilibrium claim.'},'2':{'status':'implemented_and_engine_tested','detail':'Workbench controls and saved-output checks; further research modules added here.'},'3':{'status':'deferred_by_user','detail':'Experimental data will be connected later.'},'4':{'status':'preparing_molecular_candidates','detail':'PTFE/PVDF finite molecules first; phase-specific crystal surfaces still need construction and evidence.'},'5':{'status':'blocked_physical_validation','detail':'Finite force evaluation cannot validate crystal packing, polymer force fields or water interfaces.'},'6':{'status':'preparing_dry_orientations','detail':'Six proper rotations of the capped observed core; no artificial reflection.'},'7':{'status':'blocked_model_and_preparation_gates','detail':'90-entry seed/orientation/material matrix is a queue, not executed replicas.'},'8':{'status':'blocked_validation_and_compute_budget','detail':'No production MD claimed; requires validated starts, pilots and practical GPU throughput.'},'9':{'status':'design_only_blocked','detail':'No PMF samples/free-energy values; continuous CV, cell size, restraint and convergence checks required.'},'10':{'status':'preliminary_evidence_report_and_draft','detail':'Report what exists and what is missing; no fabricated comparative evidence or external submission.'}}
  self.publish();orient=placements(self.root,self.path/'orientations');atomic(self.out/'orientations.json',orient);self.state['stages']['6']['status']='dry_coordinate_candidates_generated';self.plot();self.publish()
  matrix=[{'material':m,'orientation':o[0],'seed':s,'temperature_K':310,'status':'not_executed','eligible_now':False,'required_gates':['model_surface_and_water_validation','prepared_replica_start','finite_budget_and_benchmark'],'protein_source':'same observed capped core; these initial rotations are not independently equilibrated structures'} for m in MATERIALS for o in ORIENTATIONS for s in SEEDS];atomic(self.out/'replica-plan.json',{'replicas':matrix,'executed_replicas':0,'experimental_input_required_to_prepare':False});atomic(self.out/'pmf-plan.json',pmf_plan())
  for family in ['PTFE','PVDF']:
   out=self.data/family
   if out.exists():raise FileExistsError('Preserve existing candidates; create a new campaign path for another attempt: '+str(out))
   result=parameterize(family,out,self.amber);self.state['molecular_candidates'][family]=result;self.publish()
  self.state['stages']['4']['status']='finite_candidates_prepared_surfaces_pending';self.state['stages']['4']['detail']='Candidate coordinates/parameter attempts preserved; crystal surfaces are not yet constructed or validated.';self.state['status']='independent_preparation_done_waiting_for_PE';self.publish()
  # Waiting never launches dependent MD through an unpassed gate.
  while True:
   pe=self.manager.get(self.pe_job)
   if pe['status'] in {'completed','needs_review','failed','cancelled','interrupted'}:
    self.state['stages']['1']={'status':pe['status'],'detail':pe.get('message','')};break
   time.sleep(10)
  self.state['status']='blocked_at_scientific_and_resource_gates';self.state['stages']['7']['detail']='Existing PE preparation evaluated; comparative surface validation and practical production resources remain prerequisites. No replica MD was launched.';self.publish()


 def follow_existing(self):
  # A previous producer may still be preparing molecules. Wait for its final
  # molecular-preparation state before taking over; never rerun or overwrite it.
  while True:
   state=json.loads((self.path/'state.json').read_text())
   if state['status'] in {'blocked_at_scientific_and_resource_gates','needs_review'}:break
   self.state=state;self.state['current_PE_job']=self.pe_job;self.state['stages']['1']={'status':self.manager.get(self.pe_job)['status'],'detail':'Checkpoint continuation is active; original worker interruption is preserved separately.'};self.publish(write_private=False)
   time.sleep(10)
  if state['status']=='needs_review':return
  from md.validate_polymer_candidates import validate
  for family in ['PTFE','PVDF']:
   path=self.data/family
   if (path/'status.json').exists():validate(path)
  self.state=state;self.state['stages']['8']['detail']='Common geometric analysis implemented and checked against real saved small-model MD; production comparison remains blocked by model validation and compute budget.';self.state['current_PE_job']=self.pe_job;self.state['status']='independent_preparation_done_waiting_for_PE'
  while True:
   pe=self.manager.get(self.pe_job);self.state['stages']['1']={'status':pe['status'],'detail':pe.get('message','')};self.publish()
   if pe['status'] in {'completed','needs_review','failed','cancelled','interrupted'}:break
   time.sleep(50)
  self.state['status']='blocked_at_scientific_and_resource_gates';self.publish()

def main():
 p=argparse.ArgumentParser();p.add_argument('--pe-job',required=True);p.add_argument('--amberhome',required=True);p.add_argument('--campaign-dir',default='runs/research-campaign01');p.add_argument('--follow-existing',action='store_true');a=p.parse_args();path=(ROOT/a.campaign_dir).resolve()
 if not path.is_relative_to(ROOT/'runs'):raise SystemExit('Campaign must be inside runs/')
 if not a.follow_existing:path.mkdir(parents=True,exist_ok=False)
 c=Campaign(path,a.pe_job,a.amberhome)
 try:
  if a.follow_existing:c.follow_existing()
  else:c.run()
 except Exception as e:c.state.update(status='needs_review',error=f'{type(e).__name__}: {e}');c.publish();raise

if __name__=='__main__':main()
