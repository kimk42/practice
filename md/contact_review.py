"""Describe saved-frame geometric encounters without identifying binding sites."""
import argparse
import base64
import html
import json
from pathlib import Path
from workbench.config import ROOT
from workbench.manager import atomic


def summarize(analysis,mapping,observations):
 rows=analysis['geometry_rows'];labels=sorted({x for r in rows for x in r['contact_residues']});records=[]
 for label in labels:
  chain,seq=label.split(':');res=next(m for m in mapping if m['source_chain']==chain and str(m['source_label_seq'])==seq);caps=[m['output_residue'] for m in mapping if m['source_chain']==chain and m.get('artificial_cap')];times=[r['context_time_ps'] for r in rows if label in r['contact_residues']]
  records.append({'source_chain':chain,'source_label_seq':int(seq),'source_auth_seq':res['source_auth_seq'],'residue_name':res['name'],'saved_frames_with_contact':len(times),'total_saved_frames':len(rows),'saved_frame_fraction':len(times)/len(rows),'observed_context_times_ps':times,'model_residues_from_nearest_cap':min(abs(res['output_residue']-c) for c in caps) if caps else None})
 any_contact=[r['context_time_ps'] for r in observations if r['contact_native_residue_count']>0]
 return {'status':'saved_preparation_contacts_described','source_job':analysis['job_id'],'heavy_atom_distance_cutoff_nm':.45,'saved_frame_count':len(rows),'saved_frame_time_range_ps':analysis['saved_frame_context_time_range_ps'],'first_any_contact_count_observation_ps':min(any_contact) if any_contact else None,'identity_sampling_note':'Fine-step observations contain counts only; residue identities below come only from saved DCD frames. No interpolation of identity, contact lifetimes or off-rates.','contacts':records,'equilibrium_established':False,'binding_sites_identified':False,'limitations':['One prepared starting orientation during changing box volume; not an equilibrium sample or material comparison.','Encounter frequency is a saved-frame fraction, not equilibrium occupancy, binding strength or residence probability.','Initial placement, thermal motion and volume adjustment can contribute to contacts.','Native observed residues only; missing tails/glycans and isotope labeling are not modeled.']}


def publish(job,output):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 job=Path(job);output=Path(output);output.mkdir(parents=True,exist_ok=True);analysis_path=ROOT/'artifacts/pe-continuation/geometry-analysis.json';analysis=json.loads(analysis_path.read_text())
 if analysis['job_id']!=job.name:raise ValueError('Geometry report belongs to a different run')
 mapping=json.loads((ROOT/'data/fibrinogen-core03/manifest.json').read_text())['residue_mapping'];observations=[json.loads(x) for x in (job/'observations.jsonl').read_text().splitlines()];summary=summarize(analysis,mapping,observations);atomic(output/'contacts.json',summary);rows=analysis['geometry_rows'];fig,axes=plt.subplots(2,1,figsize=(11,6),sharex=True);times=[r['context_time_ps'] for r in rows];axes[0].plot(times,[r['minimum_protein_polymer_heavy_nm'] for r in rows],'.-');axes[0].axhline(.45,ls='--',color='orange',label='Geometric cutoff: 0.45 nm');axes[0].legend();axes[0].set_ylabel('Closest heavy-atom distance (nm)');axes[1].step(times,[r['contact_residue_count'] for r in rows],where='mid');axes[1].set(xlabel='Protein/PE context time (ps)',ylabel='Contacting native residues')
 for ax in axes:ax.grid(alpha=.2)
 fig.suptitle('Changing-volume preparation encounters — not identified binding sites');fig.tight_layout();fig.savefig(output/'contacts.png',dpi=140);plt.close(fig);image=base64.b64encode((output/'contacts.png').read_bytes()).decode();table=''.join('<tr><td>'+html.escape(str(r['source_chain'])+':'+str(r['source_label_seq']))+'</td><td>'+html.escape(r['residue_name'])+'</td><td>'+str(r['saved_frames_with_contact'])+'/'+str(r['total_saved_frames'])+'</td><td>'+str(r['source_auth_seq'])+'</td></tr>' for r in summary['contacts']);body='<h1>준비 궤적에서 관찰한 기하학적 접촉</h1><p>상자 부피가 변하는 짧은 준비 계산입니다. 결합 부위·흡착 세기·체류시간을 확정하지 않습니다. 원본 사슬과 label 번호, auth 번호를 구분합니다.</p><img style="width:100%" alt="실제 저장 프레임의 최소 중원자 거리와 접촉 잔기 수" src="data:image/png;base64,'+image+'"><table><tr><th>원본 chain:label</th><th>잔기</th><th>접촉 관찰/저장 프레임</th><th>원본 auth 번호</th></tr>'+table+'</table><pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>';(output/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><title>PE 준비 접촉 기록</title><style>body{max-width:1100px;margin:auto;padding:24px;font:16px/1.7 system-ui}pre{white-space:pre-wrap;overflow-wrap:anywhere}td,th{padding:10px;border:1px solid #ccc}table{border-collapse:collapse}</style><body>'+body+'</body></html>');return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--job',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(json.dumps(publish(ROOT/'runs/workbench'/a.job,a.output),indent=2))
