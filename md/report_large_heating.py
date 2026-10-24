"""Portable record of the executed 1 ps low-temperature large-model MD test."""
import argparse,base64,csv,html,json,os,shutil
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/large-heating-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def publish(run,output):
    run,out=map(Path,(run,output));out.mkdir(parents=True,exist_ok=True)
    r=json.loads((run/'result.json').read_text());cfg=json.loads((run/'config.json').read_text())
    with (run/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    t=[float(x['time_ps']) for x in rows]
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    axes[0,0].plot(t,[float(x['temperature_K']) for x in rows],label='Recorded');axes[0,0].step([0,.2,.4,1],[10,50,100,100],where='post',ls='--',label='Target');axes[0,0].set_ylabel('Global kinetic\ntemperature (K)');axes[0,0].legend()
    for ax,key,label in [(axes[0,1],'potential_including_position_restraints_kJ_mol','Potential energy\n(kJ/mol; support included)'),(axes[1,0],'max_relative_constraint_error','Max relative\nbond-constraint error'),(axes[1,1],'protein_CA_RMSD_from_heating_start_A','Observed-core CA RMSD\nfrom start (Å)'),(axes[2,1],'minimum_native_protein_PE_distance_nm','Observed protein–PE\nmin heavy distance (nm)')]:
        ax.plot(t,[float(x[key]) for x in rows]);ax.set_ylabel(label)
    axes[1,0].set_yscale('log');axes[1,0].axhline(1e-4,ls=':',color='red',label='Safety check limit');axes[1,0].legend(fontsize=8)
    for key,label in [('max_protein_heavy_shift_from_heating_start_nm','Protein/caps'),('max_PE_carbon_shift_from_heating_start_nm','PE carbons')]:axes[2,0].plot(t,[float(x[key]) for x in rows],label=label)
    axes[2,0].set_ylabel('Max restrained atom shift\nfrom start (nm)');axes[2,0].legend(fontsize=8)
    for ax in axes.flat:
        ax.set_xlabel('Executed large-system MD time (ps)');ax.grid(alpha=.2)
        for v in [.2,.4]:ax.axvline(v,ls=':',color='gray',lw=.8)
    fig.suptitle('Actual 577,890-atom low-temperature NVT test — 1 ps, restraints retained')
    fig.text(.5,.015,'Numerical diagnostic only: 10 → 50 → 100 K; not 310 K heating, density adjustment, equilibrium or adsorption.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.045,1,.95));fig.savefig(out/'heating.png',dpi=150);fig.savefig(out/'heating.svg');plt.close(fig)
    image=base64.b64encode((out/'heating.png').read_bytes()).decode();data=base64.b64encode((run/'observables.csv').read_bytes()).decode();f=r['final'];initial=r['initial'];g=r['geometry_after']
    body=f'''<h1>큰 모델: 저온 가열 수치 검사</h1><p><b>{r['atoms']:,}원자 · 실제 MD {r['steps']:,}스텝 · {r['time_ps']:g} ps · 1 fs</b></p><div class="notice"><b>10 → 50 → 100 K의 짧은 위치 억제 시험을 완료했습니다.</b> 수치적 검사는 통과했지만, 310 K 가열 완료·밀도 조정·계면 평형화·생산 계산 준비 완료를 뜻하지 않습니다. 관측 코어의 결손 말단·당쇄 생략과 유한 PE 패치의 한계는 그대로입니다.</div><h2>실행 방법</h2><p>이전 최소화에서 검증한 좌표·주기 상자·원래 위치 억제 기준을 유지했습니다. 결합 길이를 허용 오차 10⁻⁶으로 맞추고, seed 45로 새로운 10 K 속도를 초기화했습니다. LangevinMiddle 적분기, 마찰 5 ps⁻¹, CPU 4스레드, HBonds와 강체 TIP3P를 사용했습니다. 단백질/캡 중원자와 PE 탄소 31,974개에 1,000 kJ/mol/nm² 위치 억제를 유지했고, 물과 Ca/Na는 위치 억제를 받지 않습니다. barostat는 없습니다.</p><table><tr><th>실행 구간</th><th>목표 온도</th><th>스텝</th></tr><tr><td>0–0.2 ps</td><td>10 K</td><td>200</td></tr><tr><td>0.2–0.4 ps</td><td>50 K</td><td>200</td></tr><tr><td>0.4–1.0 ps</td><td>100 K</td><td>600</td></tr></table><p>최소화에는 물리적 MD 시간이 없었습니다. 이번이 이 큰 모델의 첫 실제 MD이며, 이전 작은 γ절편 계의 100 ps와 구분합니다. 최소화 checkpoint의 난수열을 이어가는 방식은 아닙니다.</p><img src="data:image/png;base64,{image}" alt="큰 모델의 실제 저온 시험에서 기록한 온도·에너지·제약 오차·단백질 구조·억제 원자 이동·PE 거리"><h2>실제 측정값</h2><table><tr><th>지표</th><th>결과</th></tr><tr><td>최종 실제 온도 / 최종 목표</td><td>{f['temperature_K']:.2f} / 100 K</td></tr><tr><td>기록된 최고 온도</td><td>{r['maximum_recorded_temperature_K']:.2f} K</td></tr><tr><td>기록된 최대 상대 결합 제약 오차</td><td>{r['maximum_recorded_constraint_relative_error']:.3g}</td></tr><tr><td>최종 관측 코어 Cα RMSD</td><td>{f['protein_CA_RMSD_from_heating_start_A']:.4f} Å</td></tr><tr><td>단백질/캡 중원자 최대 이동</td><td>{f['max_protein_heavy_shift_from_heating_start_nm']:.4f} nm</td></tr><tr><td>PE 탄소 최대 이동</td><td>{f['max_PE_carbon_shift_from_heating_start_nm']:.4f} nm</td></tr><tr><td>물 산소–용질 중원자 최소 거리</td><td>{g['minimum_water_O_solute_heavy_distance_nm']:.3f} nm</td></tr><tr><td>심한 비예외 중원자 / 모든 원자 겹침 쌍</td><td>{g['nonexcluded_heavy_pairs_under_0_12nm']} / {g['nonexcluded_allatom_pairs_under_0_08nm']}</td></tr><tr><td>전체 혼합계 밀도</td><td>{f['total_density_g_cm3']:.3f} g/cm³ (상자 고정으로 일정)</td></tr><tr><td>실행 시간 / 최대 프로세스 RSS</td><td>{r['elapsed_seconds']/60:.1f}분 / {r['peak_process_RSS_GiB']:.2f} GiB</td></tr></table><p>RMSD 기준은 이 저온 MD 시작 좌표이며, 원본 결정 구조가 아닙니다. 밀도가 일정한 것은 부피를 고정했기 때문이고 밀도 평형의 증거가 아닙니다. 온도 조절과 목표 온도 변경 때문에 총 에너지 보존을 기대하는 시험도 아닙니다.</p><h2>해석과 다음 순서</h2><p>이 억제 조건의 1 ps 저온 실행에서 원자·에너지의 비정상 값, 결합 길이 제약 실패, 심한 겹침이나 큰 억제 원자 이동이 나타나지 않았습니다. 이 짧은 시간과 강한 위치 억제만으로 장기 안정성·정상 단백질 기능·PE 결정성·흡착을 검증할 수는 없습니다. 다음은 단계적으로 온도를 더 올리는 검사를 설계하고, 이후 압력·밀도 조정으로 넘어가는 것입니다. 310 K에 도달하거나 물의 밀도·계면이 평형화됐다는 결론은 내리지 않았습니다.</p><p>관측 잔기 접촉 검사는 인공 캡과 Ca를 제외하고 원본 사슬/label sequence로 매핑했습니다. 가까운 위치에서 시작한 유한 표면에 위치 억제가 걸린 모델이므로 접촉을 결합 위치나 결합 세기로 해석하지 않습니다.</p><p><a href="data:text/csv;base64,{data}" download="large-low-temperature.csv">전체 저온 관측 CSV 저장</a></p><details><summary>전체 실행 결과</summary><pre>{html.escape(json.dumps(r,ensure_ascii=False,indent=2))}</pre></details>'''
    extra=f'<tr><td>관측 단백질–PE 최소 중원자 거리</td><td>{f["minimum_native_protein_PE_distance_nm"]:.3f} nm</td></tr><tr><td>PE에서 0.45 nm 이내의 관측 잔기 수</td><td>{f["contact_native_residue_count"]}</td></tr>'
    # Insert in the measured-results table, not in the preceding temperature schedule.
    marker='<tr><td>전체 혼합계 밀도</td>'
    body=body.replace(marker,extra+marker,1)
    early=[float(x['temperature_K']) for x in rows if float(x['time_ps'])<=.2]
    tail=[float(x['temperature_K']) for x in rows if float(x['time_ps'])>=.8]
    note=f'<h2>목표 온도와 실제 온도의 차이</h2><p>목표 10 K 구간에서 실제 온도는 {min(early):.2f}–{max(early):.2f} K였습니다. 초기 퍼텐셜 에너지 변화와 함께 볼 때 잔여 힘에 의한 구조 완화가 초기 온도 상승에 기여했을 가능성이 있습니다. 이 기록만으로 원인을 단독 확정하지는 않습니다. 마지막 0.2 ps의 실제 온도 평균은 {sum(tail)/len(tail):.2f} K입니다. 상관된 짧은 표본의 기술 평균이며 열적 평형 판정이 아닙니다.</p>'
    body=body.replace('<h2>해석과 다음 순서</h2>',note+'<h2>해석과 다음 순서</h2>',1)
    body=body.replace('다음은 단계적으로 온도를 더 올리는 검사를 설계하고, 이후 압력·밀도 조정으로 넘어가는 것입니다.','다음은 100 K 유지 구간에서 실제 온도·에너지·구조 추세를 더 확인한 뒤, 온도를 올리는 단계와 이후 압력·밀도 조정을 설계하는 것입니다.')
    if (Path(__file__).resolve().parents[1]/'runs/core-pe-100K-hold01/output/result.json').exists():
        body=body.replace('다음은 100 K 유지 구간에서 실제 온도·에너지·구조 추세를 더 확인한 뒤, 온도를 올리는 단계와 이후 압력·밀도 조정을 설계하는 것입니다.','이 절은 첫 1 ps 시험의 기록입니다. 이후 추가 2 ps 100 K 유지 시험을 수행했으며, 통합 보고서 최신 결과를 참조하세요. 310 K 가열과 압력·밀도 조정은 아직 완료하지 않았습니다.')
    (out/'body.html').write_text(body);css='body{max-width:1120px;margin:30px auto;padding:20px;font:16px/1.75 system-ui;color:#183247}img{width:100%;height:auto}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd6df;padding:10px;text-align:left}.notice{background:#fff0cf;padding:18px;border-left:5px solid #bf8214}pre{white-space:pre-wrap;overflow-wrap:anywhere}h2{margin-top:36px}'
    (out/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>큰 모델 저온 MD 검사</title><style>'+css+'</style><body>'+body+'</body></html>')
    for name in ['result.json','config.json','observables.csv','validation.json']:shutil.copyfile(run/name,out/name)
    print('Published',out,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');a=p.parse_args();publish(a.run,a.output)
