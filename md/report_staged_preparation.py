"""Report executed staged preparation without claiming thermodynamic/interface equilibrium."""
import argparse,base64,csv,html,json,os,shutil
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/pe-staged-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def report(run,output):
    run=Path(run);out=Path(output);r=json.loads((run/'result.json').read_text())
    if r['status']!='staged_preparation_completed_equilibrium_not_established':raise ValueError('Selected staged run not completed')
    with (run/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    if int(rows[-1]['cumulative_step'])!=r['steps']:raise ValueError('Final sample missing')
    t=np.array([float(x['time_ps']) for x in rows]);temperature=np.array([float(x['temperature_K']) for x in rows]);density=np.array([float(x['total_density_g_cm3']) for x in rows]);volume=np.array([float(x['volume_nm3']) for x in rows])
    fig,axes=plt.subplots(2,2,figsize=(12,7))
    axes[0,0].plot(t,temperature,color='#a45a35',label='Recorded');axes[0,0].plot(t,[float(x['target_temperature_K']) for x in rows],ls='--',color='#708090',label='Target');axes[0,0].set_ylabel('Temperature (K)');axes[0,0].legend(fontsize=8)
    axes[0,1].plot(t,density,color='#198991');axes[0,1].set_ylabel('TOTAL system density (g/cm³)')
    axes[1,0].plot(t,volume,color='#6a70b8');axes[1,0].set_ylabel('Periodic box volume (nm³)')
    axes[1,1].plot(t,[float(x['min_protein_PE_distance_nm']) for x in rows],color='#7b55a0');axes[1,1].set_ylabel('Closest protein–PE heavy distance (nm)')
    boundaries=[s['end_time_ps'] for s in r['stages'][:-1]]
    axes[0,1].text(.03,.08,'Final flat segment: NVT fixes volume',transform=axes[0,1].transAxes,fontsize=9,color='#526372')
    left=0
    for stage in r['stages']:
        right=stage['end_time_ps']
        for ax in axes.flat:ax.axvspan(left,right,color='#dbeaf3' if stage['ensemble']=='NVT' else '#f7e5ce',alpha=.25,zorder=-1)
        left=right
    for ax in axes.flat:
        for v in boundaries:ax.axvline(v,color='#aab3bb',ls=':',lw=1)
        ax.set_xlabel('New staged preparation time (ps)');ax.grid(alpha=.2)
    fig.suptitle('Small PE/protein/water system: executed preparation stages',fontsize=16)
    fig.text(.5,.055,'Restrained NVT heating → unrestrained isotropic NPT density adjustment → restrained NVT interface check',ha='center',fontsize=9)
    fig.text(.5,.02,'43,288 atoms | Short preparation only; equilibrium/production readiness and adsorption are NOT established.',ha='center',fontsize=10,color='#a12d32')
    fig.tight_layout(rect=(0,.09,1,.94));out.mkdir(parents=True,exist_ok=True);fig.savefig(out/'stages.png',dpi=160);fig.savefig(out/'stages.svg')
    # Descriptive drift is not a statistical equilibrium test.
    npt=[x for x in rows if 'NPT' in x['stage']];end=float(npt[-1]['time_ps']);tail=[x for x in npt if float(x['time_ps'])>=end-2]
    drift=(float(tail[-1]['total_density_g_cm3'])/float(tail[0]['total_density_g_cm3'])-1)*100
    r['descriptive_last_2ps_NPT_density_change_percent']=drift
    r['density_interpretation']='Total mixture density, not a bulk-water-only density measurement; short-window drift does not establish equilibrium.'
    (out/'result.json').write_text(json.dumps(r,indent=2));shutil.copyfile(run/'observables.csv',out/'observables.csv');shutil.copyfile(run/'config.json',out/'config.json')
    image=base64.b64encode((out/'stages.png').read_bytes()).decode();csv_data=base64.b64encode((out/'observables.csv').read_bytes()).decode()
    labels={'01-restrained-NVT-heating':'NVT 단계 가열','02-unrestrained-NPT-density':'NPT 밀도 조정','03-restrained-interface-NVT-check':'NVT 계면 점검'}
    tr=''.join(f"<tr><td>{html.escape(labels.get(s['name'],s['name']))}</td><td>{s['executed_steps']}</td><td>{s['final_temperature_K']:.1f} K</td><td>{s['final_density_g_cm3']:.4f}</td><td>{s['final_volume_nm3']:.2f}</td></tr>" for s in r['stages'])
    page=f'''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>작은 PE 계 단계적 준비 결과</title><style>body{{max-width:1100px;margin:32px auto;padding:0 20px;font:16px/1.7 system-ui;color:#172d39}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;border-bottom:1px solid #ddd;text-align:left}}img{{width:100%}}.note{{padding:18px;background:#fff0dc;border-left:4px solid #b77729}}pre{{white-space:pre-wrap}}a{{color:#087b83}}</style><h1>작은 PE–단백질–물 계: 1차 단계적 준비</h1><p>43,288원자 · 새로운 준비 구간 {r['time_ps']:g} ps · 2 fs 시간 간격 · 이전 1 ps 결과 좌표에서 새로 최소화·속도 초기화</p><p class="note"><b>단계를 실행한 것과 평형화가 확인된 것은 다릅니다.</b> 이 결과는 평형·생산 계산 준비 완료를 뜻하지 않습니다. 단백질은 인공 말단이 있는 작은 γ절편이며, 새 여섯 사슬 코어 모델의 결과가 아닙니다. 중화 이온만 포함하고 150 mM 염 농도는 설정하지 않았습니다.</p><table><tr><th>단계</th><th>실행 스텝</th><th>끝 순간 온도</th><th>끝 전체 밀도 g/cm³</th><th>끝 부피 nm³</th></tr>{tr}</table><img src="data:image/png;base64,{image}" alt="실제로 실행한 단계의 온도, 전체 밀도, 부피와 단백질–PE 거리"><h2>압력·밀도 조정의 의미</h2><p>NPT에서 유한 PE 패치의 부피 조정을 방해하지 않도록 고정 좌표 위치 제약을 모두 해제했습니다. 1 bar 목표 압력으로 등방성 Monte Carlo barostat를 사용했습니다. 이 방법을 무한 평면 또는 막 모델에 그대로 적용하면 안 됩니다. NVT 단계에서는 부피가 고정됩니다.</p><p>마지막 NPT 2 ps 구간의 전체 밀도 변화: {drift:+.2f}%. 이 짧은 구간만으로 안정성이나 통계적 평형을 확정하지 않았습니다. 전체 혼합계 밀도는 벌크 물만의 밀도와 다릅니다.</p><h2>남은 검증</h2><p>충분한 길이의 NPT와 계면 완화를 이어가며 밀도·에너지 추세, 단백질 구조, PE 포장 및 계면/벌크 물 분포를 확인해야 합니다. 단일 접촉이나 거리 감소로 흡착·결합 세기를 판단하지 않습니다.</p><p><a href="data:text/csv;base64,{csv_data}" download="staged-preparation.csv">기록 데이터 CSV 다운로드</a></p><details><summary>상세 실행 메타데이터</summary><pre>{html.escape(json.dumps(r,ensure_ascii=False,indent=2))}</pre></details></html>'''
    (out/'report.html').write_text(page);print('Published',out);return r
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run');p.add_argument('output');a=p.parse_args();report(a.run,a.output)
