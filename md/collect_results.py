"""Build a portable results report and a source-preserving results archive."""
from pathlib import Path
import base64, hashlib, html, json, zipfile
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/collected-results'
OUT.mkdir(exist_ok=True)
for body_path,raw_path in [
    ('artifacts/large-hold/body.html','runs/core-pe-100K-hold01/output/result.json'),
    ('artifacts/large-heating/body.html','runs/core-pe-low-temperature01/output/result.json'),
    ('artifacts/large-minimization/body.html','runs/core-pe-minimization02/output/result.json')]:
    if (ROOT/body_path).exists() and not (ROOT/raw_path).exists():
        raise SystemExit('Raw results are required to rebuild this snapshot. Run: python scripts/restore_raw_results.py /path/to/results-complete.zip')

# Re-render completed follow-up results from current analysis code before archiving.
pipeline=ROOT/'runs/pe-followup01/pipeline.json'
if pipeline.exists() and json.loads(pipeline.read_text()).get('status')=='completed':
    from report_pe_followup import publish
    publish(ROOT/'runs/pe-followup01',ROOT/'artifacts/pe-followup',ROOT/'data/fibrinogen-core-pe-solvated01')

min_run=ROOT/'runs/core-pe-minimization02/output'
if not (min_run/'result.json').exists():min_run=ROOT/'runs/core-pe-minimization01/output'
if (min_run/'result.json').exists():
    from report_large_minimization import publish as publish_minimization
    publish_minimization(min_run,ROOT/'artifacts/large-minimization')
heat_run=ROOT/'runs/core-pe-low-temperature01/output'
if (heat_run/'result.json').exists():
    from report_large_heating import publish as publish_heating
    publish_heating(heat_run,ROOT/'artifacts/large-heating')
hold_run=ROOT/'runs/core-pe-100K-hold01/output'
if (hold_run/'result.json').exists():
    from report_large_hold import publish as publish_hold
    publish_hold(hold_run,ROOT/'artifacts/large-hold')
def image(path,caption):
    payload=base64.b64encode((ROOT/path).read_bytes()).decode()
    return f'<figure><img src="data:image/png;base64,{payload}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'
def download(path,label):
    p=ROOT/path
    payload=base64.b64encode(p.read_bytes()).decode()
    return f'<a download="{p.parent.name}-{p.name}" href="data:application/octet-stream;base64,{payload}">{html.escape(label)}</a>'
st=json.loads((ROOT/'artifacts/staged-preparation/result.json').read_text())
sections=[]
sections.append('''<h1>피브리노겐–폴리머 시뮬레이션: 지금까지의 결과</h1><p>이 보고서에는 구조 확인, 실제 실행한 계산, 새 입력 모델 준비 결과를 모았습니다. 그림과 주요 데이터가 파일 안에 들어 있어 인터넷 없이 열 수 있습니다.</p><div class="notice"><b>현재 도달한 지점:</b> PE와 단백질 조각을 포함한 43,288원자 계에서 짧은 준비 MD를 실행했습니다. 평형 흡착, 결합 세기, 재료별 우열을 판단할 단계에는 아직 도달하지 않았습니다.</div>
<nav><a href="#overview">전체 현황</a> · <a href="#structure">입력 구조</a> · <a href="#small">실제 PE 계산</a> · <a href="#stages">이전 12 ps</a> · <a href="#core">여섯 사슬 모델</a> · <a href="#scope">해석 범위</a> · <a href="#files">원본 파일</a></nav>
<h2 id="overview">1. 무엇을 어디까지 했나</h2><table><tr><th>작업</th><th>규모</th><th>실행 범위</th></tr>
<tr><td>계산 엔진 확인</td><td>2입자</td><td>1 ps. 실제 단백질·폴리머 결과가 아님</td></tr>
<tr><td>3GHG γ사슬 절편 단독</td><td>252잔기 / 3,897원자</td><td>암묵적 용매에서 0.1 ps MD</td></tr>
<tr><td>γ절편 + PE + 물 + Na⁺</td><td>43,288원자</td><td>첫 1 ps MD</td></tr>
<tr><td>같은 작은 혼합계의 단계적 준비</td><td>43,288원자</td><td>재최소화·속도 재설정 후 별도 12 ps. 평형 미확립</td></tr>
<tr><td>여섯 사슬 관측 코어</td><td>1,913잔기 / 30,215원자</td><td>말단 캡·수소·Ca 및 힘장 준비, 국소 구조 정리</td></tr>
<tr><td>여섯 사슬 코어 + PE</td><td>물 없는 81,275원자</td><td>당시 물 없는 입력 준비 결과. 후속 용매화는 위 새 결과 참조; 큰 모델 MD 미실행</td></tr>
<tr><td>큰 모델의 물 포함 규모</td><td>약 69만 원자 추정</td><td>용질 부피 제외 전 추정. 실제 생성한 원자 수가 아님</td></tr>
<tr><td>PTFE / PVDF α·β·γ상</td><td>미구축</td><td>설정 구상만 있음. 실제 MD·비교 결과 없음</td></tr></table>''')
sections.append('''<h2 id="structure">2. 사용한 피브리노겐 구조</h2><p>사용자가 제공한 <b>3GHG.cif</b>에서 출발했습니다. 사람 피브리노겐의 2.90 Å X선 구조입니다. Assembly 1의 A–F 여섯 사슬에는 관측된 아미노산 1,913개가 있습니다. 긴 αC 및 기타 말단의 결손이 있어 완전한 전체 길이 모델은 아닙니다.</p><p>작은 계산에는 C 사슬의 원본 label 번호 143–394 절편만 사용했습니다. 새 큰 모델에는 여섯 사슬을 사용하되 당쇄와 knob 펩타이드를 제외하고, Ca²⁺ 4개와 원본 이황화결합 25개(사슬 간 15개)를 보존했습니다.</p>'''+image('artifacts/fibrinogen-3ghg.png','3GHG 구조 확인 그림. 구조 관측 범위와 결손을 함께 고려해야 합니다.')+'''<p>구조 출처: <a href="https://www.rcsb.org/structure/3GHG">RCSB 3GHG</a> · <a href="https://doi.org/10.1021/bi802205g">원 논문</a>. 사람 α사슬과 연결된 PDB 40건을 검토한 기록과 1LT9 비교 자료도 압축 파일에 포함했습니다.</p>''')
sections.append('''<h2>3. 단백질 단독 계산</h2><p>γ절편 252잔기, 3,897원자에 Amber14/GBn2 암묵적 용매를 적용해 100스텝(0.1 ps)을 실행했습니다. 물·폴리머·중화이온은 없습니다. 최소화 구조 대비 최종 Cα RMSD는 약 <b>0.209 Å</b>였습니다. 초기 구조를 계산할 수 있음을 확인한 짧은 준비 실행입니다.</p>'''+image('artifacts/protein-pilot/results.png','단백질 단독 계산: 단계별 퍼텐셜 에너지와 잔기별 변위. 총 에너지는 결합 에너지가 아닙니다.')+download('artifacts/protein-pilot/result.json','계산 요약 JSON')+' · '+download('artifacts/protein-pilot/observables.csv','관측 데이터 CSV'))
sections.append('''<h2 id="small">4. PE를 포함한 첫 실제 혼합계 계산</h2><p>γ절편 + PE 올리고머 48개 + 물 11,945분자 + Na⁺ 4개로 <b>43,288원자</b>를 만들었습니다. PE는 유한 길이 C₂₄H₅₀ 사슬이며 검증된 벌크 결정 표면은 아닙니다. 단백질 ff14SB, PE GAFF2/AM1-BCC, 물 TIP3P를 사용했습니다. PE 원시 전하의 반올림 잔차(+0.006e)는 원자별 균등 보정했으며 원본과 보정 기록을 보존했습니다.</p><p>명시적 물에서 1 ps MD를 실행했습니다. 단백질–PE 최소 중원자 거리는 약 <b>0.505 → 0.427 nm</b>로 변했습니다. 마지막 접촉은 원본 label Thr393 근처였지만 인공 절단 말단에 가까우므로 피브리노겐의 고유 결합 위치로 해석할 수 없습니다. 접촉 기준은 PE 탄소에서 0.45 nm 이내의 단백질 중원자입니다.</p>'''+image('artifacts/pe-pilot/results.png','첫 PE 혼합계 계산: 실제 배치와 거리·접촉 변화.')+download('artifacts/pe-pilot/observables.csv','첫 PE 관측 CSV')+' · '+download('artifacts/pe-pilot/result.json','첫 PE 요약 JSON'))
sections.append(f'''<h2 id="stages">5. 이전 실행: 물·밀도·계면의 단계적 준비 12 ps</h2><p>앞의 작은 계 좌표를 가져와 다시 최소화하고 100 K에서 속도를 재설정했습니다. 앞선 1 ps와 동일 조건으로 연속 실행한 13 ps 궤적은 아닙니다. 2 fs 간격으로 총 6,000스텝을 실행했습니다.</p><table><tr><th>단계</th><th>길이</th><th>조건</th><th>마지막 온도</th></tr><tr><td>가열 NVT</td><td>2 ps</td><td>100 → 200 → 310 K, 단백질·PE 중원자 위치 억제</td><td>245.1 K</td></tr><tr><td>밀도 조정 NPT</td><td>8 ps</td><td>310 K / 1 bar, 위치 억제 해제</td><td>307.4 K</td></tr><tr><td>계면 점검 NVT</td><td>2 ps</td><td>단백질 약한 억제, PE 아래층 억제</td><td>310.3 K</td></tr></table><p>전체 혼합계 밀도 <b>0.785 → 0.928 g/cm³</b>, 부피 <b>550.25 → 465.41 nm³</b>. 최종 최소 거리는 <b>{st['final_min_distance_nm']:.3f} nm</b>, 접촉 잔기의 원본 label 번호는 280, 393입니다. 이는 짧은 준비 구간의 접촉 관측이며 안정 결합 부위를 뜻하지 않습니다.</p><div class="notice">NPT 마지막 2 ps에서도 밀도가 약 1.63% 변했습니다. 마지막 NVT에서 밀도가 일정한 것은 부피를 고정했기 때문입니다. 평형화 완료·생산 계산 준비 완료로 판단하지 않았습니다. 이온은 중화용 Na⁺만 있으며 150 mM 염 조건은 아닙니다.</div>'''+image('artifacts/staged-preparation/stages.png','단계별 온도·전체 혼합계 밀도·부피·단백질–PE 거리. NVT의 평평한 밀도는 평형의 증거가 아닙니다.')+download('artifacts/staged-preparation/observables.csv','이전 12 ps 관측 CSV')+' · '+download('artifacts/staged-preparation/result.json','이전 12 ps 결과 JSON')+' · '+download('runs/pe-staged-preparation01/validation.json','실행 검증 JSON')+'''<p>입자 수, 단계별 barostat·억제 조건, 결합 길이 제약 오차와 60개 궤적 프레임을 검증했습니다. 최신 실행은 CPU 4스레드에서 준비·파일 출력을 포함해 약 11분 소요됐습니다.</p>''')
sections.append('''<h2 id="core">6. 다음 규모를 위한 여섯 사슬 코어 + PE</h2><p>관측된 여섯 사슬에 ACE/NME 캡 12개를 붙여 인공 말단 전하를 줄였습니다. 결손 사슬을 복원한 것은 아닙니다. 국소 정리 후 1.2 Å 미만의 비결합 원자 쌍은 18개에서 0개로 줄었지만, 용액 평형화 결과는 아닙니다.</p><p>여기에 PE 올리고머 690개를 배치해 물 없는 <b>81,275원자</b> 입력 모델을 만들었습니다. 이 건조 입력의 순전하는 약 −10e입니다. 후속 용매화·중화 결과는 추가 작업 절을 참조하세요. 예상 상자는 약 52.8 × 13.0 × 8.9 nm이고, 물을 채우면 약 69만 원자 규모가 될 것으로 추정합니다. 용질 부피를 빼기 전 추정치이므로 실제 생성 원자 수는 달라집니다.</p>'''+image('artifacts/core-model/model.png','새 여섯 사슬 관측 코어와 PE 표면 입력 모델. 물 없는 모델 준비 결과이며 흡착 MD 결과가 아닙니다.')+download('artifacts/core-model/summary.json','코어·PE 모델 검증 요약 JSON'))
sections.append('''<h2 id="scope">7. 지금 답할 수 있는 것과 남은 것</h2><table><tr><th>관심 질문</th><th>현재 결과가 말해 주는 범위</th></tr><tr><td>결합 위치</td><td>짧은 γ절편 계산의 접촉 잔기만 관측. 전체 단백질 결합 위치 확정 불가</td></tr><tr><td>결합 세기</td><td>미계산. 총 퍼텐셜 에너지는 흡착 자유에너지가 아님</td></tr><tr><td>표면 흡착</td><td>가까이 배치한 작은 계의 준비 거동 확인. 자발 흡착·평형 흡착 검증 없음</td></tr><tr><td>PE vs PTFE vs PVDF 결정상</td><td>PE만 실제 모델·MD 있음. 비교 결론 없음</td></tr></table><p>다음 순서는 작은 PE 계에서 더 긴 압력·밀도 조정과 물/PE/단백질 구조 검증 → 억제 조건의 단계적 해제 → 여러 배향·시드 반복입니다. 큰 코어 계는 용매화·중화부터 별도로 진행해야 합니다. 그 뒤 PTFE와 PVDF α·β·γ의 표면 모델을 검증하고 같은 조건에서 비교할 수 있습니다.</p>
<h2>8. 초기 시험과 실패 기록</h2><p>2입자 계산은 계산 엔진 확인용입니다. 이전 데이터에 PVDF 등의 설정 문구가 있어도 실제 폴리머 계산이 아닙니다. 단백질 첫 실행 중단, PE 전하·잔기명·중복 결합 수정, 코어 캡 이름 및 국소 정리 실패 시도는 성공 결과에 포함하지 않았습니다. 압축 파일의 실패 기록은 추적용으로 별도 경로에 보존했습니다.</p>'''+image('artifacts/engine-result.png','초기 2입자 엔진 시험 — 과학적 흡착 결과로 사용하지 않음.'))
sections.append('''<h2 id="files">9. 파일을 보는 방법</h2><p><b>이 HTML:</b> 브라우저로 열면 모든 그림을 볼 수 있습니다. 주요 CSV·JSON은 위 링크에서 저장할 수 있습니다. CSV는 엑셀에서 열 수 있고 PDB/CIF 구조는 ChimeraX·PyMOL 등에서 볼 수 있습니다.</p><p><b>results-complete.zip:</b> 이 보고서, 기존 개별 보고서, 그래프, 관측 데이터, 업로드 원본, 올바른 최종 입력 모델, MD 궤적·상태·체크포인트, 실행 스크립트와 설정, 실패 기록을 포함합니다. 압축을 풀면 원래의 상대 경로가 유지됩니다. MD 프로그램과 Python 가상환경 자체는 제외했습니다. 재현에는 별도 환경 설치가 필요하며 체크포인트는 실행 환경에 의존할 수 있습니다.</p><p>압축 파일의 FILES.json에는 포함 파일의 크기와 SHA-256이 기록되어 있습니다. 경로의 로컬 환경 이름은 출처 추적용이며 보고서 열기에 필요하지 않습니다.</p>''')
# Include the new follow-up report as a portable section when it exists.
followup=ROOT/'artifacts/pe-followup/body.html'
if followup.exists():
    sections.insert(1,'<hr><section id="followup">'+followup.read_text()+'</section><hr>')
    sections[0]=sections[0].replace('<nav>','<nav><a href="#followup">추가 100 ps와 후속 작업</a> · ',1)
    additions=[]
    for label,path,key in [('추가 NPT','runs/pe-followup01/01-npt-100ps/result.json','time_ps'),('계면 억제 조정·해제','runs/pe-followup01/03-interface-release/result.json','time_ps')]:
        p=ROOT/path
        if p.exists():
            r=json.loads(p.read_text());additions.append(f'<tr><td>{label}</td><td>{r["atoms"]:,}원자</td><td>{r[key]:g} ps 실행. 평형 미확립</td></tr>')
    p=ROOT/'data/fibrinogen-core-pe-solvated01/manifest.json'
    if p.exists():
        r=json.loads(p.read_text());additions.append(f'<tr><td>큰 코어 + PE + 실제 물·중화이온</td><td>{r["atoms"]:,}원자</td><td>용매화·중화 입력 준비 기록. 후속 MD는 최신 결과 참조</td></tr>')
    sections[0]=sections[0].replace('</table>',''.join(additions)+'</table>',1)
    latest=ROOT/'runs/pe-followup01/pipeline.json'
    if latest.exists() and json.loads(latest.read_text()).get('status')=='completed':
        sections[0]=sections[0].replace('43,288원자 계에서 짧은 준비 MD를 실행했습니다.','43,288원자 계의 100 ps 추가 NPT와 10 ps 계면 준비를 완료하고, 여섯 사슬 큰 모델에 물과 중화이온을 넣었습니다.')
min_body=ROOT/'artifacts/large-minimization/body.html'
if min_body.exists():
    sections.insert(1,'<hr><section id="large-minimization">'+min_body.read_text()+'</section><hr>')
    sections[0]=sections[0].replace('<nav>','<nav><a href="#large-minimization">큰 모델 최소화 최신 결과</a> · ',1)
    m=json.loads((min_run/'result.json').read_text())
    label='위치 억제 최소화 실행 및 예비 검사 통과' if m['allow_heating_preparation'] else '최소화 실행 후 추가 검토 필요'
    row=f'<tr><td>큰 모델 위치 억제 최소화</td><td>{m["atoms"]:,}원자</td><td>{label}. 이 최소화 단계는 MD 0스텝</td></tr>'
    sections[0]=sections[0].replace('</table>',row+'</table>',1)
    sections[0]=sections[0].replace('평형 흡착, 결합 세기, 재료별 우열을 판단할 단계에는 아직 도달하지 않았습니다.','큰 모델의 위치 억제 최소화도 실행했습니다. 평형 흡착, 결합 세기, 재료별 우열을 판단할 단계에는 아직 도달하지 않았습니다.',1)
heat_body=ROOT/'artifacts/large-heating/body.html'
if heat_body.exists():
    sections.insert(1,'<hr><section id="large-heating">'+heat_body.read_text()+'</section><hr>')
    sections[0]=sections[0].replace('<nav>','<nav><a href="#large-heating">큰 모델 저온 MD 최신 결과</a> · ',1)
    row='<tr><td>큰 모델 위치 억제 저온 MD</td><td>577,890원자</td><td>10→50→100 K 목표, 1 ps 실제 MD. 310 K·평형 미확립</td></tr>'
    sections[0]=sections[0].replace('</table>',row+'</table>',1)
    sections[0]=sections[0].replace('큰 모델의 위치 억제 최소화도 실행했습니다.','큰 모델의 위치 억제 최소화와 1 ps 저온 MD 검사도 실행했습니다.',1)
hold_body=ROOT/'artifacts/large-hold/body.html'
if hold_body.exists():
    sections.insert(1,'<hr><section id="large-hold">'+hold_body.read_text()+'</section><hr>')
    sections[0]=sections[0].replace('<nav>','<nav><a href="#large-hold">100 K 추가 유지 최신 결과</a> · ',1)
    row='<tr><td>큰 모델 100 K 유지</td><td>577,890원자</td><td>추가 2 ps, 큰 모델 MD 누적 3 ps. 압력·밀도 조정과 평형 미확립</td></tr>'
    sections[0]=sections[0].replace('</table>',row+'</table>',1)
    sections[0]=sections[0].replace('큰 모델의 위치 억제 최소화와 1 ps 저온 MD 검사도 실행했습니다.','큰 모델의 위치 억제 최소화, 1 ps 저온 시험과 추가 2 ps 100 K 유지(큰 모델 MD 누적 3 ps)도 실행했습니다.',1)
# Embed the latest executed workbench report and downloadable CSV in the full snapshot.
workbench_report=ROOT/'artifacts/workbench/preparation/report.html'
if workbench_report.exists():
    body=workbench_report.read_text().split('<body>',1)[-1].rsplit('</body>',1)[0]
    sections.insert(1,'<hr><section id="workbench-preparation">'+body+'</section><hr>')
    sections[0]=sections[0].replace('<nav>','<nav><a href="#workbench-preparation">실행기 기반 준비 계산 최신 상태</a> · ',1)
guide=ROOT/'artifacts/data-guide.md'
if guide.exists():
    sections.append('<h2>데이터·단위·파일 해설</h2><details open><summary>각 데이터와 지표를 읽는 방법</summary><pre style="white-space:pre-wrap;font-family:inherit">'+html.escape(guide.read_text())+'</pre></details>')
css='body{font-family:system-ui,sans-serif;line-height:1.75;max-width:1120px;margin:auto;padding:28px;color:#182d42;background:#f6f8fb}h1,h2{color:#123b62}h2{margin-top:48px}table{width:100%;border-collapse:collapse;background:white}td,th{border:1px solid #cbd6e0;padding:10px;text-align:left}th{background:#e3edf7}img{width:100%;height:auto;background:white}figure{margin:24px 0}figcaption{color:#526477;font-size:14px}.notice{background:#fff0cf;border-left:5px solid #d49014;padding:18px}a{color:#075cba}nav{margin:24px 0}p{overflow-wrap:anywhere}'
(OUT/'report.tmp.html').write_text('<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>피브리노겐–폴리머 전체 결과</title><style>'+css+'</style><body>'+''.join(sections)+'</body></html>',encoding='utf-8')
(OUT/'report.tmp.html').replace(OUT/'report.html')
# Preserve successful final bundles and all run/attempt records; omit environments and caches.
paths=set()
for dirname in ['artifacts','runs','data','md','environment','tests','schemas','docs','scripts','workbench','.github']:
    d=ROOT/dirname
    if d.exists():
        for p in d.rglob('*'):
            if p.is_file() and OUT not in p.parents and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name not in {'archive-validation.json','result-archive.json'} and 'cache' not in p.parts and not p.is_symlink() and p.suffix!='.tmp' and p.name!='.active': paths.add(p)
for name in ['README.md','requirements.txt','package.json','package-lock.json','index.html','app.mjs','simulation.mjs','workbench.mjs','model.mjs','styles.css','.gitignore','.gitattributes']:
    if (ROOT/name).is_file(): paths.add(ROOT/name)
paths.add(OUT/'report.html')
entries=[]
archive=OUT/'results-complete.zip'
archive_temp=OUT/'results-complete.tmp.zip'
with zipfile.ZipFile(archive_temp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(paths):
        rel=str(p.relative_to(ROOT)); b=p.read_bytes()
        z.writestr(rel,b)
        entries.append({'path':rel,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
    z.writestr('OPEN-FIRST.txt','먼저 artifacts/collected-results/report.html을 브라우저로 여세요.\n실패 시도는 추적용이며 성공 계산과 구분하세요.\n')
    manifest={'generated_utc':datetime.now(timezone.utc).isoformat(),'files':entries,'excluded':['Python virtual environment','AmberTools executables','cache files','Git internals'],'note':'Includes unsuccessful preparation attempts for provenance; use the successful paths identified in report and README.'}
    z.writestr('FILES.json',json.dumps(manifest,ensure_ascii=False,indent=2))
archive_temp.replace(archive)
(OUT/'FILES.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print(json.dumps({'report_bytes':(OUT/'report.html').stat().st_size,'archive_bytes':archive.stat().st_size,'file_count':len(entries)},indent=2))
