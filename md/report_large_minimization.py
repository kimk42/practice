"""Portable report of executed large-system restrained minimization, no MD claim."""
import argparse,base64,csv,html,json,os,shutil
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/large-min-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def publish(run,output):
    run,out=map(Path,(run,output));out.mkdir(parents=True,exist_ok=True)
    r=json.loads((run/'result.json').read_text());cfg=json.loads((run/'config.json').read_text())
    with (run/'minimization.csv').open() as f:rows=list(csv.DictReader(f))
    fig,axes=plt.subplots(2,2,figsize=(12,7));x=[int(v['callback']) for v in rows]
    for ax,key,label in [(axes[0,0],'physical_plus_position_restraint_kJ_mol','Force-field + position restraint energy (kJ/mol)'),(axes[0,1],'constraint_penalty_kJ_mol','Internal constraint penalty (kJ/mol)'),(axes[1,0],'max_relative_constraint_error','Max relative constraint error'),(axes[1,1],'objective_gradient_RMS_kJ_mol_nm','Objective gradient RMS (kJ/mol/nm)')]:
        ax.plot(x,[float(v[key]) for v in rows]);ax.set(xlabel='Recorded minimizer callback (not MD steps)',ylabel=label);ax.grid(alpha=.2)
    for ax in [axes[0,1],axes[1,0],axes[1,1]]:ax.set_yscale('symlog',linthresh=1e-7)
    fig.suptitle(f'Actual {r["atoms"]:,}-atom restrained periodic minimization — no dynamics')
    fig.text(.5,.02,'Optimizer iteration may reset during constraint repair. Internal penalty is separate from positional support and binding energy.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.95));fig.savefig(out/'minimization.png',dpi=150);fig.savefig(out/'minimization.svg');plt.close(fig)
    picture=base64.b64encode((out/'minimization.png').read_bytes()).decode();csv_data=base64.b64encode((run/'minimization.csv').read_bytes()).decode()
    b,a=r['before'],r['after'];status='최소화 실행 및 예비 검사 통과' if r['allow_heating_preparation'] else '최소화 실행 후 추가 검토 필요'
    body=f'''<h1>큰 모델: 원자 겹침 검사와 위치 억제 에너지 최소화</h1><p><b>{r['atoms']:,}원자</b> · 여섯 사슬 관측 코어 + PE + TIP3P 물 + 중화이온 · {html.escape(status)}</p><div class="notice"><b>MD는 0스텝입니다.</b> 이번 결과는 고정된 주기 상자에서 배치를 정리한 에너지 최소화입니다. 가열·밀도 조정·계면 평형화·흡착이나 결합 자유에너지 계산이 아닙니다. 최소화의 완전 수렴도 인증하지 않았습니다.</div><h2>실행 조건</h2><p>Amber ff14SB 단백질, GAFF2/AM1-BCC PE, TIP3P 물. PME 1 nm cutoff, HBonds 및 강체 물, 원래 원자 질량을 사용했습니다. 단백질 중원자와 모든 PE 탄소에 1,000 kJ/mol/nm² 위치 억제를 걸고, 수소·물·Ca²⁺·Na⁺는 결합 제약 외에는 자유롭게 두었습니다. barostat 없이 고정 상자에서 CPU 4스레드로 실행했습니다.</p><p>제약 보정 시도마다 최적화 반복 상한 {cfg['max_iterations_per_constraint_repair_attempt']}회, 기록된 callback {r['recorded_minimizer_callbacks']}개. 반복 번호는 제약 보정 때 다시 0으로 돌아갈 수 있어 callback 수나 반복을 MD 시간으로 환산하면 안 됩니다.</p><h2>실제 결과</h2><table><tr><th>지표</th><th>결과</th></tr><tr><td>결합 등 예외 쌍을 제외한 중원자 거리 &lt; 0.12 nm</td><td>{b['nonexcluded_heavy_pairs_under_0_12nm']} → {a['nonexcluded_heavy_pairs_under_0_12nm']}쌍</td></tr><tr><td>예외 쌍을 제외한 모든 원자 거리 &lt; 0.08 nm</td><td>{b['nonexcluded_allatom_pairs_under_0_08nm']} → {a['nonexcluded_allatom_pairs_under_0_08nm']}쌍</td></tr><tr><td>물 산소–용질 중원자의 최소 거리</td><td>{b['minimum_water_O_solute_heavy_distance_nm']:.3f} → {a['minimum_water_O_solute_heavy_distance_nm']:.3f} nm</td></tr><tr><td>제약 적용 후 시작 → 최종 총 에너지</td><td>{r['initial_constrained_total_energy_kJ_mol']:,.1f} → {r['final_total_energy_kJ_mol']:,.1f} kJ/mol</td></tr><tr><td>최종 위치 억제 에너지</td><td>{r['final_position_restraint_energy_kJ_mol']:,.1f} kJ/mol</td></tr><tr><td>단백질 중원자 최대 이동</td><td>{r['max_protein_heavy_displacement_nm']:.4f} nm</td></tr><tr><td>PE 탄소 최대 이동</td><td>{r['max_PE_carbon_displacement_nm']:.4f} nm</td></tr><tr><td>최종 최대 상대 결합 길이 제약 오차</td><td>{r['max_relative_constraint_error']:.3g}</td></tr><tr><td>실행 시간 / 프로세스 최대 RSS</td><td>{r['elapsed_seconds']/60:.1f}분 / {r['peak_process_RSS_GiB']:.2f} GiB</td></tr></table><p>거리 검사는 주기 경계의 최소 이미지에서 시행하며 NonbondedForce 예외 쌍(결합·1–3·1–4 포함)을 제외합니다. 이 임계값은 심한 겹침을 찾는 기준이며 모든 입체 충돌이나 단백질 구조의 정상성을 검증한 것은 아닙니다.</p><img src="data:image/png;base64,{picture}" alt="실제로 기록한 최소화 에너지·제약 페널티·제약 오차·최적화 기울기"><p>여기서 내부 제약 페널티는 최적화기가 결합 길이를 유지하려고 사용하는 항입니다. 단백질·PE 위치 억제 에너지와 별개입니다. 총 에너지 감소도 결합 세기나 흡착 자유에너지를 뜻하지 않습니다.</p><h2>현재 해석과 다음 작업</h2><p>고정 상자이므로 전체 밀도는 초기 {r['initial_total_density_g_cm3']:.3f} g/cm³에서 변하지 않습니다. 물의 밀도·온도·압력 평형을 확인한 결과가 아닙니다. 다음에는 낮은 온도에서 짧은 위치 억제 가열을 실행하고, 힘·제약·단백질/PE 구조를 확인한 뒤 압력·밀도 조정을 설계해야 합니다. 불완전한 관측 코어·당쇄 생략·유한 PE 패치의 한계는 계속 적용됩니다.</p><p><a href="data:text/csv;base64,{csv_data}" download="large-minimization.csv">최소화 기록 CSV 저장</a></p><details><summary>전체 결과 메타데이터</summary><pre>{html.escape(json.dumps(r,ensure_ascii=False,indent=2))}</pre></details>'''
    if cfg.get('restart_state'):
        note='<div class="notice"><b>재실행 경위:</b> 첫 설정(보정 시도마다 최대 200회)은 30분 안에 결합 길이 제약 보정을 끝내지 못해 중단됐습니다. 첫 시도를 성공 또는 수렴으로 처리하지 않았습니다. 보존한 중간 좌표를 불러와 결합 길이를 재투영한 뒤, 시도마다 최대 20회로 짧게 최적화하면서 제약 보정을 더 자주 적용한 두 번째 시도입니다. 첫 시도의 전체 CSV·실패 로그·중간 상태도 원본 ZIP에 보존했습니다. 아래 에너지 시작값은 재실행 직전 제약을 적용한 좌표의 값입니다.</div>'
        body=body.replace('<h2>실행 조건</h2>',note+'<h2>실행 조건</h2>',1)
    if (Path(__file__).resolve().parents[1]/'runs/core-pe-low-temperature01/output/result.json').exists():
        body=body.replace('다음에는 낮은 온도에서 짧은 위치 억제 가열을 실행하고, 힘·제약·단백질/PE 구조를 확인한 뒤 압력·밀도 조정을 설계해야 합니다.','이 절은 최소화 단계의 기록입니다. 이후 수행한 저온 MD 검사는 통합 보고서의 최신 결과에 기록합니다. 310 K 가열과 압력·밀도 조정은 아직 완료하지 않았습니다.')
    (out/'body.html').write_text(body)
    css='body{max-width:1120px;margin:30px auto;padding:20px;font:16px/1.75 system-ui;color:#183247;background:#f7f9fc}img{width:100%;height:auto}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd6df;padding:10px;text-align:left}.notice{background:#fff0cf;padding:18px;border-left:5px solid #bf8214}pre{white-space:pre-wrap;overflow-wrap:anywhere}h2{margin-top:36px}'
    (out/'report.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>큰 모델 위치 억제 최소화 결과</title><style>'+css+'</style><body>'+body+'</body></html>')
    for name in ['result.json','config.json','minimization.csv','geometry-before.json']:shutil.copyfile(run/name,out/name)
    print('Published',out,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');a=p.parse_args();publish(a.run,a.output)
