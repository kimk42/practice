"""Report an executed 100 K hold with descriptive windows, no equilibrium claim."""
import argparse,base64,csv,html,json,os,shutil
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/large-hold-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def publish(run,output):
    run,out=map(Path,(run,output));out.mkdir(parents=True,exist_ok=True)
    r=json.loads((run/'result.json').read_text());cfg=json.loads((run/'config.json').read_text());check=json.loads((run/'validation.json').read_text())
    with (run/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    t=np.array([float(x['time_ps']) for x in rows]);v=lambda key:np.array([float(x[key]) for x in rows])
    windows=[]
    for left in [0,.5,1,1.5]:
        subset=[x for x in rows if left<float(x['time_ps'])<=left+.5];w={'hold_start_ps':left,'hold_end_ps':left+.5,'cumulative_start_ps':left+1,'cumulative_end_ps':left+1.5,'observations':len(subset)}
        for key in ['temperature_K','potential_including_position_restraints_kJ_mol','protein_CA_RMSD_from_hold_start_A']:
            values=np.array([float(x[key]) for x in subset]);times=np.array([float(x['time_ps']) for x in subset]);w[key]={'mean':float(values.mean()),'standard_deviation':float(values.std()),'slope_per_ps':float(np.polyfit(times-times[0],values,1)[0])}
        windows.append(w)
    (out/'windows.json').write_text(json.dumps({'windows':windows,'interpretation':'Descriptive correlated samples; no standard error, significance or equilibrium test.'},indent=2))
    with (out/'windows.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['hold_start_ps','hold_end_ps','temperature_mean_K','temperature_std_K','potential_mean_kJ_mol','potential_slope_kJ_mol_ps'])
        for w in windows:writer.writerow([w['hold_start_ps'],w['hold_end_ps'],w['temperature_K']['mean'],w['temperature_K']['standard_deviation'],w['potential_including_position_restraints_kJ_mol']['mean'],w['potential_including_position_restraints_kJ_mol']['slope_per_ps']])
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    axes[0,0].plot(t,v('temperature_K'));axes[0,0].axhline(100,ls='--',color='orange',label='Target 100 K');axes[0,0].set_ylabel('Global kinetic\ntemperature (K)');axes[0,0].legend(fontsize=8)
    for ax,key,label in [(axes[0,1],'potential_including_position_restraints_kJ_mol','Potential energy\n(kJ/mol; support included)'),(axes[1,0],'max_relative_constraint_error','Max relative\nbond-constraint error'),(axes[1,1],'protein_CA_RMSD_from_hold_start_A','Observed-core CA RMSD\nfrom hold start (Å)'),(axes[2,1],'minimum_native_protein_PE_distance_nm','Observed protein–PE\nminimum heavy distance (nm)')]:ax.plot(t,v(key));ax.set_ylabel(label)
    axes[1,0].set_yscale('log');axes[1,0].axhline(1e-4,ls=':',color='red')
    for key,label in [('max_protein_heavy_shift_from_hold_start_nm','Protein/caps'),('max_PE_carbon_shift_from_hold_start_nm','PE carbon')]:axes[2,0].plot(t,v(key),label=label)
    axes[2,0].set_ylabel('Maximum restrained atom shift\nfrom hold start (nm)');axes[2,0].legend(fontsize=8)
    for ax in axes.flat:ax.set_xlabel('Additional 100 K hold time (ps)');ax.grid(alpha=.2)
    fig.suptitle('Actual 577,890-atom 100 K hold — 2 additional ps, original restraints retained')
    fig.text(.5,.015,'Checkpoint continuation: cumulative context time 1 → 3 ps. Fixed underdense box; no 310 K, pressure or adsorption result.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.045,1,.95));fig.savefig(out/'hold.png',dpi=150);fig.savefig(out/'hold.svg');plt.close(fig)
    source=resolve_project_path(cfg['source_run'])
    with (source/'observables.csv').open() as f:previous=list(csv.DictReader(f))
    combined=[{'cumulative_time_ps':float(x['time_ps']),'target_K':float(x['target_K']),'temperature_K':float(x['temperature_K']),'potential_kJ_mol':float(x['potential_including_position_restraints_kJ_mol'])} for x in previous]
    combined.extend({'cumulative_time_ps':1+float(x['time_ps']),'target_K':100,'temperature_K':float(x['temperature_K']),'potential_kJ_mol':float(x['potential_including_position_restraints_kJ_mol'])} for x in rows[1:])
    with (out/'combined-temperature-energy.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(combined[0]));writer.writeheader();writer.writerows(combined)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,key,label in [(axes[0],'temperature_K','Global kinetic temperature (K)'),(axes[1],'potential_kJ_mol','Potential + support energy (kJ/mol)')]:
        ax.plot([x['cumulative_time_ps'] for x in combined],[x[key] for x in combined]);ax.axvline(1,ls=':',color='gray',label='Checkpoint boundary');ax.set(xlabel='Cumulative large-model MD time (ps)',ylabel=label);ax.grid(alpha=.2)
    axes[0].step([0,.2,.4,3],[10,50,100,100],where='post',ls='--',color='orange',label='Target');axes[0].legend(fontsize=8)
    fig.suptitle('Large model: original low-T test + verified checkpoint hold');fig.tight_layout();fig.savefig(out/'combined.png',dpi=150);fig.savefig(out/'combined.svg');plt.close(fig)
    embed=lambda name:'data:image/png;base64,'+base64.b64encode((out/name).read_bytes()).decode()
    last=r['final'];table=''.join(f"<tr><td>{w['hold_start_ps']:g}–{w['hold_end_ps']:g}</td><td>{w['temperature_K']['mean']:.3f} ± {w['temperature_K']['standard_deviation']:.3f}</td><td>{w['potential_including_position_restraints_kJ_mol']['slope_per_ps']:+.1f}</td></tr>" for w in windows)
    body=f'''<h1>큰 모델: 100 K 유지 구간 2 ps</h1><p><b>577,890원자 · 추가 MD 2,000스텝 / 2 ps · 큰 모델 누적 MD 3 ps</b></p><div class="notice"><b>수치 검사와 저장 파일 교차 검사를 통과했습니다.</b> 짧은 위치 억제·고정 상자 시험이며, 310 K 가열·액체 밀도·압력·계면 평형 또는 생산 계산 준비가 완료됐다는 뜻은 아닙니다.</div><h2>이전 실행에서 이어받은 조건</h2><p>이전 1 ps 저온 시험의 checkpoint를 같은 OpenMM/CPU 4스레드, 같은 최종 적분기로 불러왔습니다. 좌표·속도가 이전 XML 상태와 일치함을 확인했고, 속도 재초기화·추가 최소화·제약 재투영·난수열 재설정을 하지 않았습니다. 상자와 원래 억제 기준을 유지했습니다. 환경 간 비트 단위 동일성을 보장한다는 의미는 아닙니다.</p><p>100 K, LangevinMiddle 마찰 5 ps⁻¹, 1 fs, HBonds와 강체 TIP3P, 결합 제약 허용 오차 10⁻⁶. 단백질/캡 중원자 및 PE 탄소 31,974개에 1,000 kJ/mol/nm² 위치 억제를 유지했고 barostat는 없습니다. context 스텝은 1,000→3,000, 시간은 1→3 ps로 이어졌습니다. 아래 유지 구간 그래프와 CSV의 시간은 추가 구간 시작을 0으로 둡니다. OpenMM DCDReporter는 DCD 헤더에도 구간 상대 스텝(100→2000)을 기록합니다. DCD 상대 시간에 1 ps를 더하면 누적 context 시간에 대응합니다. 원본 DCD 헤더나 좌표를 수정하지 않았습니다.</p><img src="{embed('combined.png')}" alt="큰 모델의 첫 저온 시험과 추가 유지 구간 전체 3 ps 온도와 에너지"><h2>추가 유지 구간의 실제 기록</h2><img src="{embed('hold.png')}" alt="추가 유지 구간의 온도·에너지·제약 오차·구조·억제 원자 이동·PE 거리"><table><tr><th>지표</th><th>결과</th></tr><tr><td>시작 → 마지막 실제 온도</td><td>{r['initial']['temperature_K']:.3f} → {last['temperature_K']:.3f} K</td></tr><tr><td>마지막 Cα RMSD (유지 구간 시작 기준)</td><td>{last['protein_CA_RMSD_from_hold_start_A']:.4f} Å</td></tr><tr><td>단백질/캡 중원자 최대 이동</td><td>{last['max_protein_heavy_shift_from_hold_start_nm']:.4f} nm</td></tr><tr><td>PE 탄소 최대 이동</td><td>{last['max_PE_carbon_shift_from_hold_start_nm']:.4f} nm</td></tr><tr><td>관측 단백질–PE 최소 중원자 거리 / 0.45 nm 이내 잔기 수</td><td>{last['minimum_native_protein_PE_distance_nm']:.3f} nm / {last['contact_native_residue_count']}</td></tr><tr><td>기록된 최대 상대 결합 제약 오차</td><td>{r['maximum_recorded_constraint_relative_error']:.3g}</td></tr><tr><td>전체 혼합계 밀도 (고정 상자)</td><td>{last['total_density_g_cm3']:.3f} g/cm³</td></tr><tr><td>실행 시간 / 최대 프로세스 RSS</td><td>{r['elapsed_seconds']/60:.1f}분 / {r['peak_process_RSS_GiB']:.2f} GiB</td></tr></table><h2>0.5 ps 구간별 기술 통계</h2><table><tr><th>추가 구간 시간 ps</th><th>평균 온도 ± 표준편차 K</th><th>퍼텐셜 에너지 선형 기울기 kJ/mol/ps</th></tr>{table}</table><p>표본은 시간적으로 상관되어 있습니다. 평균·표준편차·선형 기울기는 기술 통계이며 표준오차·유의성·열적 평형 판정이 아닙니다. 에너지는 위치 억제를 포함하고, thermostat가 있어 총 에너지가 보존되는 시험이 아닙니다. 밀도가 일정한 것은 상자를 고정했기 때문이며, 초기 저밀도가 해결되었다는 뜻도 아닙니다.</p><h2>해석 범위</h2><p>이 억제 조건의 짧은 저온 구간에서 수치 한계와 저장 파일 검사를 통과한 결과입니다. 장기 구조 안정성·정상 기능·PE 결정성·결합 위치·흡착 세기를 확정하지 않았습니다. 원래의 결손 말단·당쇄 생략·유한 PE 모델·중화 이온만 있는 조건의 한계가 적용됩니다. 온도를 더 올리기 전에 이 구간의 잔여 에너지 추세와 국소 구조를 검토하고, 이후 압력·밀도 조정을 별도로 설계해야 합니다.</p><details><summary>독립 교차 검사</summary><pre>{html.escape(json.dumps(check,ensure_ascii=False,indent=2))}</pre></details><details><summary>전체 실행 결과</summary><pre>{html.escape(json.dumps(r,ensure_ascii=False,indent=2))}</pre></details>'''
    first,last_window=windows[0],windows[-1]
    energy_delta=r['final']['potential_including_position_restraints_kJ_mol']-r['initial']['potential_including_position_restraints_kJ_mol']
    measured=f'<h2>추세를 읽는 방법</h2><p>첫 0.5 ps의 평균 온도는 {first["temperature_K"]["mean"]:.3f} K, 마지막 0.5 ps는 {last_window["temperature_K"]["mean"]:.3f} K입니다. 추가 2 ps 전체의 퍼텐셜 변화는 {energy_delta:+,.1f} kJ/mol이고, 마지막 0.5 ps의 선형 기울기는 {last_window["potential_including_position_restraints_kJ_mol"]["slope_per_ps"]:+,.1f} kJ/mol/ps입니다. 이는 큰 혼합계 전체의 값이며 흡착 또는 결합 자유에너지가 아닙니다. 짧고 상관된 표본의 기울기만으로 물리적 수렴이나 평형을 판정하지 않았습니다.</p>'
    body=body.replace('<h2>해석 범위</h2>',measured+'<h2>해석 범위</h2>',1)

    links=[]
    for name,label in [('observables.csv','추가 구간 관측 CSV'),('windows.csv','0.5 ps 구간 통계 CSV'),('combined-temperature-energy.csv','큰 모델 전체 3 ps 온도·에너지 CSV')]:
        src=run/name if name=='observables.csv' else out/name;encoded=base64.b64encode(src.read_bytes()).decode();links.append(f'<a download="{name}" href="data:text/csv;base64,{encoded}">{label}</a>')
    body+='<p>'+' · '.join(links)+'</p>';(out/'body.html').write_text(body)
    css='body{max-width:1120px;margin:30px auto;padding:20px;font:16px/1.75 system-ui;color:#183247}img{width:100%;height:auto}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd6df;padding:10px;text-align:left}.notice{background:#fff0cf;padding:18px;border-left:5px solid #bf8214}pre{white-space:pre-wrap;overflow-wrap:anywhere}h2{margin-top:36px}'
    (out/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>큰 모델 100 K 유지 MD</title><style>'+css+'</style><body>'+body+'</body></html>')
    for name in ['result.json','config.json','observables.csv','validation.json']:shutil.copyfile(run/name,out/name)
    print('Published',out,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');a=p.parse_args();publish(a.run,a.output)
