"""Publish real progress and executed follow-up results; never label queued work completed."""
import argparse,base64,csv,html,json,os,shutil,time
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/pe-followup-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load(p):return json.loads(p.read_text()) if p.exists() else None

def rows(p):
    if not p.exists():return []
    with p.open() as f:return list(csv.DictReader(f))


def publish(run,output,large):
    display_time=datetime.now(ZoneInfo('America/Los_Angeles')).strftime('%Y-%m-%d %H:%M:%S %Z')
    run,out,large=map(Path,(run,output,large));out.mkdir(parents=True,exist_ok=True)
    status=load(run/'pipeline.json') or {'status':'NPT_running','current_stage':'01-npt-100ps'}
    npt=run/'01-npt-100ps';p=load(npt/'progress.json') or {};nr=load(npt/'result.json');analysis=load(run/'02-structure-check/analysis.json');interface=load(run/'03-interface-release/result.json');big=load(large/'manifest.json');obs=rows(npt/'observables.csv');ir=rows(run/'03-interface-release/observables.csv')
    sections=['<h1>추가 PE 계산과 후속 세 단계</h1>']
    labels={'NPT_running':'100 ps NPT 계산 중','structure_check_running':'구조 검사 중','interface_running':'계면 억제 조정·해제 계산 중','solvation_running':'큰 모델 용매화·중화 중','completed':'요청한 계산 및 준비 단계 완료','failed':'오류 발생 — 후속 단계 중지'}
    label=labels.get(status['status'],status['status'])
    sections.append(f'<div class="notice"><b>현재 상태: {html.escape(label)}</b><p>마지막 보고서 갱신: {display_time}. 계산이 완료됐다는 것과 평형이 확인됐다는 것은 다릅니다.</p></div>')
    items=[('100 ps NPT',bool(nr)),('구조 검사',bool(analysis)),('10 ps 계면 준비·단백질 억제 해제',bool(interface)),('여섯 사슬 큰 모델 용매화·중화',bool(big))]
    sections.append('<table><tr><th>단계</th><th>상태</th></tr>'+''.join(f'<tr><td>{x}</td><td>{"완료" if done else "미완료 / 대기 또는 실행 중"}</td></tr>' for x,done in items)+'</table>')
    if p:
        elapsed=p.get('extension_time_ps',0);remaining=p.get('estimated_remaining_seconds')
        sections.append(f'<h2>1. NPT 연장</h2><p>실제 실행 기록: <b>{elapsed:g} / 100 ps</b>. 43,288원자 · 2 fs · 310 K · 1 bar · CPU 4스레드 · 위치 억제 없음. 이전 마지막 NVT의 좌표·속도·상자를 유지하고 thermostat/barostat 난수열을 새로 시작했습니다. 이전 궤적과 비트 단위로 동일한 연속 계산은 아닙니다.</p>')
        if remaining and not nr:sections.append(f'<p>최근 속도로 추정한 남은 시간: 약 {remaining/60:.0f}분. 실행 속도와 후속 단계 시간에 따라 달라집니다.</p>')
    if len(obs)>1:
        fig,axes=plt.subplots(3,2,figsize=(12,10));t=[float(x['extension_time_ps']) for x in obs]
        for ax,(key,label) in zip(axes.flat,[('temperature_K','Temperature (K)'),('total_density_g_cm3','TOTAL mixture density (g/cm³)'),('volume_nm3','Box volume (nm³)'),('potential_kJ_mol','Potential energy (kJ/mol)'),('min_protein_PE_distance_nm','Minimum protein–PE heavy distance (nm)'),('CA_RMSD_from_extension_start_A','CA RMSD from extension start (Å)')]):
            ax.plot(t,[float(x[key]) for x in obs]);ax.set_xlabel('Additional NPT time (ps)');ax.set_ylabel(label);ax.grid(alpha=.2)
        fig.suptitle('Executed NPT extension — finite PE patch / gamma fragment / explicit water')
        fig.tight_layout(rect=(0,.02,1,.96));fig.savefig(out/'npt.png',dpi=140);fig.savefig(out/'npt.svg');plt.close(fig)
        encoded=base64.b64encode((out/'npt.png').read_bytes()).decode();sections.append(f'<img src="data:image/png;base64,{encoded}" alt="실제로 기록한 NPT 결과 그래프 6개">')
        last=obs[-1];sections.append(f'<p>마지막 기록: 온도 {float(last["temperature_K"]):.1f} K, 전체 밀도 {float(last["total_density_g_cm3"]):.4f} g/cm³, 최소 단백질–PE 거리 {float(last["min_protein_PE_distance_nm"]):.3f} nm. 총 퍼텐셜 에너지는 결합 자유에너지가 아닙니다.</p>')
    if analysis:
        sections.append('<h2>2. 구조 검사</h2><p>검사 상태: '+html.escape(analysis['status'])+'. 결합 길이 제약, 물–용질 겹침, 단백질 Cα 크기, PE 사슬 길이·배향 및 물 분포를 점검했습니다. 통과는 후속 준비 계산을 허용하는 기하학적 검사이며 평형·결정성·단백질 기능 검증은 아닙니다.</p>')
        a,b=analysis['initial_structure'],analysis['final_structure']
        sections.append(f'<p>Cα 회전반경 {a["CA_radius_of_gyration_nm"]:.3f} → {b["CA_radius_of_gyration_nm"]:.3f} nm. PE 사슬 평균 끝간 거리 {a["PE_chain_end_to_end_mean_nm"]:.3f} → {b["PE_chain_end_to_end_mean_nm"]:.3f} nm. 최대 상대 결합 제약 오차 {b["max_relative_constraint_error"]:.2g}.</p>')
        fig,ax=plt.subplots(figsize=(10,4))
        for name,snap in [('Initial',a),('After 100 ps NPT',b)]:
            profile=snap['water_O_z_profile'];ax.plot(profile['bin_center_nm'],profile['oxygen_number_density_nm_minus3'],label=name)
        ax.set(xlabel='Z relative to PE carbon circular mean (nm)',ylabel='Water O number density (nm⁻³)',title='Finite patch: water profile along initial slab Z axis');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'water-profile.png',dpi=140);plt.close(fig)
        sections.append('<img src="data:image/png;base64,'+base64.b64encode((out/'water-profile.png').read_bytes()).decode()+'" alt="물 산소의 Z방향 분포"><p>이 분포에는 유한 패치의 가장자리와 용질이 차지하는 공간이 섞여 있습니다. 순수한 벌크 물 밀도나 정밀 계면 구조의 검증 결과로 해석하지 않습니다.</p>')
        structural_rows=rows(run/'02-structure-check/trajectory-structure.csv')
        if structural_rows:
            fig,axes=plt.subplots(1,2,figsize=(12,4));t=[float(x['extension_time_ps']) for x in structural_rows]
            for ax,key,label in [(axes[0],'CA_radius_of_gyration_nm','Protein CA radius of gyration (nm)'),(axes[1],'PE_nematic_order_max_eigenvalue','PE chain orientational order')]:
                ax.plot(t,[float(x[key]) for x in structural_rows]);ax.set(xlabel='Additional NPT time (ps)',ylabel=label);ax.grid(alpha=.2)
            fig.tight_layout();fig.savefig(out/'trajectory-structure.png',dpi=140);plt.close(fig)
            sections.append('<img src="data:image/png;base64,'+base64.b64encode((out/'trajectory-structure.png').read_bytes()).decode()+'" alt="실제 궤적의 단백질 크기와 PE 배향 변화"><p>실제 DCD 궤적 100프레임을 읽어 구조 변화를 확인했습니다. 회전반경과 배향 지표만으로 단백질의 정상 기능이나 PE 결정성을 확정할 수는 없습니다.</p>')
        sections.append(f'<p>PE 사슬 중심의 평균 최근접 거리 {a["PE_nearest_chain_COM_distance_mean_nm"]:.3f} → {b["PE_nearest_chain_COM_distance_mean_nm"]:.3f} nm. 최종 다른 사슬 간 최소 탄소 거리는 {b["PE_min_interchain_carbon_distance_nm"]:.3f} nm입니다. 국소 포장 변화의 진단 수치이며 벌크 PE 결정의 검증을 대신하지 않습니다.</p>')
        sections.append('<details><summary>10 ps 구간별 밀도·에너지 기록</summary><pre>'+html.escape(json.dumps(analysis['descriptive_10ps_windows'],ensure_ascii=False,indent=2))+'</pre></details><p>구간 평균·표준편차·기울기는 상관된 표본의 기술 통계이며, 평형 판정이나 통계적 유의성을 의미하지 않습니다.</p>')
    if interface:
        sections.append('<h2>3. 계면 억제 조정·해제</h2><p>밀도 조정 후 상자를 고정하고 NVT 10 ps를 실행했습니다. 새 끝 좌표에 임시 단백질 위치 억제를 걸어 100 → 10 → 0 kJ/mol/nm²로 2 → 2 → 6 ps 동안 조정했습니다. 아래 PE 탄소의 지지 억제 1,000은 유지했고 위 PE는 자유롭게 두었습니다. 앞 NPT에서는 이미 모든 위치 억제가 풀려 있었으므로 별도의 지지 표면 진단 단계입니다.</p>')
        if ir:
            fig,axes=plt.subplots(1,2,figsize=(12,4));t=[float(x['interface_time_ps']) for x in ir]
            axes[0].plot(t,[float(x['min_protein_PE_distance_nm']) for x in ir]);axes[0].set_ylabel('Minimum protein–PE distance (nm)');axes[1].plot(t,[float(x['CA_RMSD_from_NPT_endpoint_A']) for x in ir]);axes[1].set_ylabel('CA RMSD from NPT endpoint (Å)')
            for ax in axes:ax.set_xlabel('Interface preparation time (ps)');ax.grid(alpha=.2);ax.axvline(2,ls=':',color='gray');ax.axvline(4,ls=':',color='gray')
            fig.tight_layout();fig.savefig(out/'interface.png',dpi=140);plt.close(fig);sections.append('<img src="data:image/png;base64,'+base64.b64encode((out/'interface.png').read_bytes()).decode()+'" alt="계면 준비의 거리 및 구조 변화">')
        sections.append('<p>마지막 접촉 잔기 원본 label 번호: '+html.escape(str(interface['final_contact_residue_source_labels']))+'. 짧은 접촉 관측이며 전체 피브리노겐의 결합 부위 또는 결합 세기로 확정하지 않습니다.</p>')
    if big:
        sections.append(f'<h2>4. 여섯 사슬 큰 모델 용매화·중화</h2><p>실제로 생성한 크기: <b>{big["atoms"]:,}원자</b>, 물 <b>{big["water_molecules"]:,}분자</b>. 순전하 {big["net_charge_e"]:.3g}e. 원본 이황화결합 25개와 건조 모델 원자 순서를 보존했습니다. 이 수치는 이전의 약 69만 원자 추정과 구분되는 실제 생성 결과입니다.</p><p>물과 중화이온을 넣은 입력 준비까지입니다. 큰 모델의 에너지 최소화·평형화·흡착 MD는 실행하지 않았습니다. 긴 결손 αC와 당쇄는 생략된 모델이며 중화 이온만 넣었습니다.</p>')
    if big:
        sections.append(f'<p>새 큰 모델의 초기 전체 밀도는 {big["initial_total_density_g_cm3"]:.3f} g/cm³입니다. 초기 물 배치에는 용질 주변 배제 공간이 있으며, 아직 용액 속 최소화·압력/밀도 조정을 거치지 않았습니다. 이전 약 69만 원자 추정은 용질 부피 등을 제외하기 전의 거친 계산이었고, 실제 물 배치·용질 주변 배제 및 중화 후에는 {big["atoms"]:,}원자가 생성됐습니다.</p>')
    if big and big.get('visualization'):
        image_path=Path(__file__).resolve().parents[1]/big['visualization']
        sections.append('<img src="data:image/png;base64,'+base64.b64encode(image_path.read_bytes()).decode()+'" alt="실제로 용매화한 여섯 사슬 코어와 PE 및 물·나트륨 좌표"><p>그림의 물 산소 점은 표시용으로 일부만 골랐습니다. 전체 좌표는 solvated.cif와 Amber 입력 파일에 보존했습니다. 그림은 용매화 입력이며 큰 모델 MD 결과가 아닙니다.</p>')
    sections.append('<h2>해석 범위와 전체 데이터</h2><p>PE의 유한 올리고머 표면과 인공 절단 말단이 있는 γ절편의 준비 계산입니다. 결합 자유에너지·자발 흡착·평형 흡착·PTFE/PVDF 비교는 계산하지 않았습니다. 큰 모델 준비와 작은 계 MD는 별개입니다.</p><p>최종 전체 자료는 artifacts/collected-results/report.html 및 results-complete.zip에 병합합니다. 원본·관측 CSV·구조·궤적·입력·검증·실행 코드·설명을 포함하며, 이전 결과를 덮어쓰지 않습니다.</p>')
    data_files=[(npt/'observables.csv','추가 NPT 관측 CSV'),(run/'02-structure-check/trajectory-structure.csv','궤적 구조 지표 CSV'),(run/'02-structure-check/water-profile.csv','물 분포 CSV'),(run/'03-interface-release/observables.csv','계면 준비 관측 CSV')]
    links=[]
    for path,label in data_files:
        if path.exists():
            encoded=base64.b64encode(path.read_bytes()).decode()
            links.append(f'<a download="{path.parent.name}-{path.name}" href="data:text/csv;base64,{encoded}">{label}</a>')
    sections.append('<p>'+' · '.join(links)+'</p>')
    if status.get('error'):sections.append('<pre>'+html.escape(status['error'])+'</pre>')
    if any((Path(__file__).resolve().parents[1]/f'runs/core-pe-minimization{i:02d}/output/result.json').exists() for i in (1,2)):
        sections=[s.replace('큰 모델의 에너지 최소화·평형화·흡착 MD는 실행하지 않았습니다.','이 절은 용매화 단계의 기록입니다. 이후 최소화는 최신 결과를 참조하세요. 큰 모델의 가열·평형화·흡착 MD는 아직 실행하지 않았습니다.') for s in sections]
    if (Path(__file__).resolve().parents[1]/'runs/core-pe-low-temperature01/output/result.json').exists():
        sections=[s.replace('큰 모델의 가열·평형화·흡착 MD는 아직 실행하지 않았습니다.','이후 큰 모델의 1 ps 저온 시험을 수행했습니다. 최신 결과를 참조하세요. 310 K 가열·평형화·흡착 계산은 아직 완료하지 않았습니다.') for s in sections]
    body=''.join(sections);(out/'body.html').write_text(body)
    css='body{max-width:1120px;margin:30px auto;padding:20px;font:16px/1.75 system-ui;color:#183247;background:#f7f9fc}img{width:100%;height:auto}table{border-collapse:collapse;width:100%;background:white}td,th{border:1px solid #ccd6df;padding:10px;text-align:left}.notice{background:#fff0cf;padding:18px;border-left:5px solid #bf8214}pre{white-space:pre-wrap;overflow-wrap:anywhere}h2{margin-top:40px}'
    temp=out/'report.tmp.html';temp.write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>추가 계산 진행·결과 보고서</title><style>'+css+'</style><body>'+body+'</body></html>');temp.replace(out/'report.html')
    for src in [npt/'observables.csv',npt/'config.json',npt/'progress.json']:
        if src.exists():shutil.copyfile(src,out/src.name)
    return status

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('output');p.add_argument('large');a=p.parse_args();publish(a.run,a.output,a.large)
