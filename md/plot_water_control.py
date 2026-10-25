"""Portable plots from completed water-only controls, kept separate from protein MD."""
import argparse
import base64
import csv
import html
import json
from pathlib import Path
import numpy as np
from workbench.manager import atomic


def publish(raw,out):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 raw=Path(raw);out=Path(out);out.mkdir(parents=True,exist_ok=True);result=json.loads((raw/'result.json').read_text());series=[];diagnostics=[]
 for replica in result['replicas']:
  rows=[json.loads(x) for x in (raw/f'seed-{replica["seed"]}'/'observations.jsonl').read_text().splitlines()];series.append(rows);npt=[r for r in rows if r['ensemble']=='NPT'];pressure=[r['pressure_bar'] for r in npt if r['pressure_bar'] is not None];diagnostics.append({'seed':replica['seed'],'pressure_samples':len(pressure),'sample_pressure_mean_bar':float(np.mean(pressure)),'sample_pressure_std_bar':float(np.std(pressure)),'pressure_interpretation':'Sparse correlated instantaneous values, not a converged pressure estimate or an uncertainty on the mean.','density_result':replica})
  with (out/f'seed-{replica["seed"]}.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 fig,axes=plt.subplots(3,1,figsize=(11,9),sharex=True)
 for rows in series:
  npt=[r for r in rows if r['ensemble']=='NPT'];label='seed '+str(rows[0]['seed']);axes[0].plot([r['additional_NPT_ps'] for r in npt],[r['density_g_cm3'] for r in npt],label=label);axes[1].plot([r['additional_NPT_ps'] for r in npt],[r['temperature_K'] for r in npt],label=label);p=[r for r in npt if r['pressure_bar'] is not None];axes[2].plot([r['additional_NPT_ps'] for r in p],[r['pressure_bar'] for r in p],'.-',ms=3,label=label)
 axes[0].set_ylabel('Pure-water density (g/cm³)');axes[1].set_ylabel('Temperature (K)');axes[2].set_ylabel('Instantaneous pressure (bar)');axes[2].set_xlabel('Control NPT time (ps)');axes[1].axhline(310,color='gray',ls='--');axes[2].axhline(1,color='gray',ls='--')
 for ax in axes:ax.grid(alpha=.2);ax.legend()
 fig.suptitle('512-water implementation controls — not fibrinogen/PE MD or equilibrium validation');fig.tight_layout();fig.savefig(out/'controls.png',dpi=140);plt.close(fig)
 summary={'status':result['status'],'replicas':diagnostics,'source_conditions':result['source_conditions'],'protein_MD_time_added_ps':0,'equilibrium_established':False,'scope':result['scope']};atomic(out/'control-summary.json',summary);image=base64.b64encode((out/'controls.png').read_bytes()).decode();body='<h1>물 단독 구현 대조</h1><p>물 512분자씩, 독립 초기 배향·속도 시드 두 개. 각 5 ps NVT와 100 ps NPT. 단백질·PE는 포함하지 않으며 그 누적 시간에 더하지 않습니다. 물의 전하·LJ·질량·cutoff·dispersion 조건은 실제 PE 원본과 대조했습니다.</p><p>액체에 가까운 밀도와 기본 압력 조절 동작의 기술 검사입니다. 폴리머 표면·단백질·열역학적 평형을 검증한 결과가 아닙니다. 압력 표본의 표준편차는 평균의 오차가 아닙니다.</p><img style="width:100%" alt="두 물 단독 대조의 실제 밀도·온도·압력" src="data:image/png;base64,'+image+'"><pre>'+html.escape(json.dumps(summary,ensure_ascii=False,indent=2))+'</pre>'
 downloads=[]
 for file in sorted(out.glob('*.csv')):downloads.append('<a download="'+file.name+'" href="data:application/octet-stream;base64,'+base64.b64encode(file.read_bytes()).decode()+'">'+file.name+'</a>')
 body+='<h2>실제 관측 CSV</h2><p>'+' · '.join(downloads)+'</p>'
 (out/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>물 단독 구현 대조</title><style>body{max-width:1100px;margin:auto;padding:24px;font:16px/1.7 system-ui}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style><body>'+body+'</body></html>');return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--raw',required=True);p.add_argument('--output',required=True);a=p.parse_args();print(json.dumps(publish(a.raw,a.output),indent=2))
