# Protein · Surface Explorer

피브리노겐과 PTFE/PVDF/PE 표면 상호작용 연구를 위한 기초 구조 탐색기입니다. 브라우저 탐색기는 외부 패키지 없이 실행하고, MD 실행기는 OpenMM을 사용합니다.

## 실행

```sh
cd /workspace/practice
npm start
```

기본 포트는 8000입니다. 테스트: `npm test` (Node.js 18 이상). 서버: Python 3.

브라우저에서 PDB 파일을 열면 첫 모델의 Cα 좌표를 표시합니다. 파일은 업로드되지 않습니다. 캔버스를 드래그하거나 회전 슬라이더를 사용하세요. 외부 구조 다운로드는 필요하지 않습니다. RCSB PDB에서 직접 받은 구조를 사용할 수 있습니다.

## 해석 범위

- 기본 모식도는 실제 피브리노겐 좌표가 아닙니다.
- 원본 PDB의 종, 조립체, 사슬, 누락 영역을 확인하세요. 구조 복원과 생물학적 조립체 생성은 지원하지 않습니다.
- Cα–평면 거리만 계산합니다. 재료 선택은 설명과 색상만 바꿉니다.
- 브라우저 구조 화면은 MD를 실행하지 않습니다. 별도 MD 실행기는 아래 절차를 사용하세요. 도킹, 결합 에너지, 흡착 자유에너지 계산은 미구현이며 현재 기능으로 재료별 흡착 순위를 도출할 수 없습니다.
- 다음 단계: 검증된 구조 및 표면 모델 선택, 힘장 검토, 용매 포함 MD와 반복 계산.

## MD 실행기 (신규)

PVDF α/β/γ상, 분극 조건, 노출 결정면을 UI에서 지정하고 `simulation.json`을 저장합니다. PTFE/PE의 결정 구조는 모델 자체에서 명시합니다. 결정면의 지수는 사용한 결정 셀 기준입니다. α상은 이 버전에서 비극성 배열만 지원합니다. β/γ상의 비극성·미분극 조건은 도메인 구성에 대한 별도 모델링이 필요합니다. 상 선택만으로 원자 배열이나 전하가 생성되지는 않습니다.

```sh
cd /workspace/practice
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m unittest discover -s tests -p 'test_*.py'
.venv/bin/python md/run.py simulation.json prepared-model runs/first --platform CPU
```

`prepared-model/`에는 다음 세 파일이 필요합니다.

- `coordinates.pdb`: 동일한 원자 순서의 단백질/고분자/물/이온 좌표와 토폴로지.
- `system.xml`: OpenMM으로 직렬화한 완전히 파라미터화된 System. 주기 경계·힘장·전하·구속·표면 유지 조건은 이 파일에 포함해야 합니다. 실행기가 자동으로 추가하지 않습니다.
- `manifest.json`: 아래 형식의 모델 출처와 파일 SHA-256. `surface`는 UI에서 내보낸 조건과 정확히 일치해야 합니다.

```json
{
  "surface": {"polymer":"PVDF","phase":"beta","face":"(100)","polarization":"plus-normal"},
  "protein_source":"실제로 사용한 구조 출처와 조립체",
  "polymer_forcefield":"실제로 사용한 파라미터와 문헌/버전",
  "protein_forcefield":"실제로 사용한 힘장/버전",
  "water_model":"실제로 사용한 물 모델",
  "preparation_notes":"용매화, 이온 농도, 결함/누락 처리, 결정 셀, 표면 배향 및 도메인 구성",
  "sha256":{"system.xml":"실제 SHA-256","coordinates.pdb":"실제 SHA-256"}
}
```

새 출력 디렉터리에 `trajectory.dcd`, `energy.csv`, `final.pdb`, `result.json`을 저장합니다. 기존 출력은 덮어쓰지 않습니다. CPU는 2스레드를 사용하고 CUDA/OpenCL도 선택할 수 있습니다. 100회 이하의 에너지 최소화 후 지정 온도로 속도를 초기화하는 단기 NVT 실행기로, 평형화와 생산 계산을 대체하지 않습니다. 시스템 총 퍼텐셜 에너지는 결합 에너지가 아닙니다.

**현재 미구현:** 검증된 PVDF/PTFE/PE 원자 모델 생성 및 힘장 통합, 피브리노겐 준비, 반복 흡착 계산, 궤적 기반 접촉/구조 변화 분석과 자유에너지 계산. 모델 메타데이터와 체크섬 검사는 물리적 타당성을 인증하지 않습니다. 엔진 통합 테스트는 합성 2입자계이며 단백질–고분자 연구 검증이 아닙니다.

## 실제 피브리노겐 구조 준비

```sh
.venv/bin/python md/prepare_protein.py fetch 3GHG data/fibrinogen-candidate-3GHG-new
.venv/bin/python md/prepare_protein.py audit data/fibrinogen-candidate-3GHG-new/source.pdb data/fibrinogen-candidate-3GHG-new/audit.json
```

3GHG는 확인할 후보 ID이며 구조 다운로드·메타데이터 검토 전에는 최종 피브리노겐 모델로 확정하지 않습니다. 다운로드는 원본 PDB 및 RCSB 엔트리 메타데이터를 저장합니다. 연결 실패 시 `download.json`에 실패 원인을 남기고 종료 코드 2로 실패합니다. 출력 폴더를 재사용하지 않아 기존 데이터를 덮어쓰지 않습니다.

점검은 사슬·누락 잔기·누락 중원자·비표준 잔기·원본 헤더와 SHA-256을 기록합니다. 원본에 없는 긴 구간을 임의 생성하거나 원자를 삭제하지 않습니다. 점검 보고서는 MD 준비 완료를 의미하지 않습니다. 종, biological assembly, 당화, 이황화 결합, 누락 영역의 취급을 먼저 검토해야 합니다.

현재 실측 차단 기록: `data/fibrinogen-candidate-3GHG/download.json`. `files.rcsb.org`의 HTTPS CONNECT 요청이 중계 서버에서 403으로 거부되어 실제 구조 확보 및 흡착 계산은 진행되지 않았습니다. 두 입자 시험계는 이 작업을 대체하지 않습니다.

## 업로드된 실제 3GHG

원본 `data/3ghg-upload/source.cif`는 수정하지 않고 보존했습니다. `audit.json`은 PDBFixer의 원본 전체 점검이며 두 조립체를 모두 포함합니다. `assembly1-summary.json`은 첫 조립체의 단백질 6사슬만을 대상으로 한 관측 Cα/전체 서열 비교입니다. 둘의 잔기 수 범위는 다릅니다.

`.venv/bin/python md/import_3ghg.py`는 mmCIF에 명시된 첫 조립체가 항등 변환인지 검사하고, 관측된 단백질 Cα 좌표만 보기용 PDB 및 이미지로 내보냅니다. `assembly1-view.pdb`는 힘장 입력이 아닙니다. 원본의 당쇄·칼슘·knob는 원본에 보존하고 보기에서만 생략합니다. 브라우저는 이 실제 구조를 기본으로 불러옵니다. 원본 mmCIF의 출처와 누락 영역은 반드시 유지해야 합니다. 사슬별 누락 목록은 summary에 기록되어 있습니다.

## 실제 γ사슬 조각 CPU 계산

```sh
OPENMM_CPU_THREADS=2 .venv/bin/python md/protein_pilot.py data/3ghg-upload/source.cif runs/3ghg-gamma-pilot-new --steps 1000
```

첫 조립체의 γ사슬 C, label 서열 143–394의 연속된 252잔기를 대상으로 합니다. 절단한 조각의 인공적인 말단, 다른 사슬과 당쇄·칼슘의 부재라는 근사 조건이 있습니다. 원본 3GHG는 변경하지 않습니다. 없는 잔기 구간을 자동 복원하지 않습니다. 중원자 누락이 발견되면 중단하고 말단 원자와 수소만 추가합니다. 출력 PDB의 사슬/잔기 번호는 원본과 다를 수 있으므로 원본 범위와 대조해야 합니다.

Amber14 + GBn2 implicit solvent, 310 K, 1 fs, CPU 최대 4스레드, 결합 제약 없이 전체 결합항을 유지하고 25회 한도의 최소화 후 짧은 MD를 실행합니다. 명시적 물·염이나 폴리머가 포함된 계산이 아닙니다. `observables.csv`에는 총 퍼텐셜 에너지, 순간 온도, 최소화 구조 대비 정렬한 Cα RMSD가 기록됩니다. RMSD는 강체 이동과 회전을 제거합니다. 짧은 계산으로 평형화·흡착·결합 세기를 판단하지 않습니다. 실제 폴리머 흡착 단계는 별도의 검증된 표면 구조와 힘장이 필요합니다.

## PE를 포함한 실제 혼합계 파일과 실행

이 단계는 실제 3GHG γ사슬 조각 + 유한 PE 올리고머 표면 + 명시적 TIP3P 물 + 중화 이온의 짧은 준비 계산입니다. 무한 결정 표면이나 평형 흡착을 재현한 모델로 해석하지 않습니다.

```sh
.venv/bin/python md/install_ambertools.py
.venv/bin/python md/prepare_pe.py data/pe-gaff2/new-oligomer
.venv/bin/python md/build_pe_complex.py runs/3ghg-gamma-pilot-02/final.pdb data/pe-gaff2/new-oligomer data/new-pe-protein-bundle
OPENMM_CPU_THREADS=4 .venv/bin/python md/run_pe_complex.py data/new-pe-protein-bundle runs/new-pe-protein-pilot --steps 1000
.venv/bin/python md/plot_pe_complex.py data/new-pe-protein-bundle runs/new-pe-protein-pilot artifacts/pe-pilot
```

AmberTools 24.8의 Linux x86-64 설치는 `environment/ambertools-linux64.explicit`의 SHA-256으로 고정된 패키지를 사용합니다. 설치 도구 자체도 conda-forge 아카이브 SHA-256을 검사합니다. 사용자 홈에는 환경을 등록하지 않으며 `/workspace/tools/ambertools`에 설치합니다. 설치 재실행 시 기존 버전을 확인하고 재설치하지 않습니다.

모델 준비:

- C24H50(12 반복 단위) PE 올리고머: RDKit로 초기 all-trans 구조 생성, UFF는 초기 구조 완화에만 사용. 실제 MD에는 AmberTools의 GAFF2(데이터 파일 version 2.2.20)와 AM1-BCC 전하를 적용.
- SQM 출력의 전하 반올림 잔차를 검사하고, 0.01 e 이하일 때만 알려진 중성 전하를 만족하도록 원자 전체에 균등 보정. 원본 `PE.raw.mol2`, 보정 크기, 보정 후 파일과 SHA-256을 보존. 실제 준비에서는 합계 +0.006 e, 원자당 −0.000081081 e 보정. 큰 잔차는 실패 처리.
- 2열 × 12행 × 2층, 총 48사슬의 **유한 올리고머 패치**. 대략적인 평행 배열이지 실험으로 검증한 PE 결정상은 아님. 메틸 말단과 패치 가장자리가 있음.
- 실제 단백질 조각은 ff14SB 사용. 히스티딘 수소 상태, 이황화 결합 및 말단 이름을 Amber 템플릿과 대응. PDB CONECT와 명시적 이황화 결합 명령이 중복되지 않게 처리.
- TIP3P 물과 전하 중화에 필요한 이온만 추가. 150 mM 염 농도나 생리적 용액 조건을 재현하지 않음.
- 첫 실제 준비 파일: `data/pe-gaff2/model02/`, `data/pe-protein-model04/`. 앞선 번호의 실패 산출물을 완료 모델로 사용하지 않음.

실행과 해석:

- PME, 1 nm 비결합 cutoff, HBonds 제약, 310 K Langevin, 1 fs, NVT. PE 아래층 탄소만 1000 kJ/mol/nm² 위치 제약. 위층과 단백질은 자유롭게 움직임.
- 최소화의 `maxIterations=50`은 OpenMM에서 결합 제약 재조정 시 여러 번 수행될 수 있음. 최소화 수렴이나 압력·밀도 평형화를 보장하지 않음.
- `observables.csv`의 최소 거리는 주기 경계를 고려한 단백질 중원자–PE 탄소 최근접 거리. 접촉 잔기는 0.45 nm 이내 중원자를 하나 이상 가진 단백질 잔기이며 원본 γ사슬 label 서열 번호로 기록. 이는 수소결합이나 결합 에너지의 정의가 아님.
- 단백질은 미리 표면 가까이 배치됨. 가까워짐 또는 접촉 형성만으로 자발적 흡착을 증명하지 않음. 전체 피브리노겐이나 PTFE/PVDF와의 비교로 일반화하지 않음.
- 결과는 원자 수·중성 전하·체크섬·유한 좌표를 검증하고, 실제 실행된 스텝과 결과 파일을 보존. 궤적, 체크포인트, 최종 상태도 남김.
- 다음 과학적 단계는 표면 모델 검증, 압력/밀도 및 계면 평형화, 여러 초기 배향·시드·시간 길이 반복, PTFE와 각 PVDF 결정상의 별도 파라미터 검토, 자유에너지 계산임.

첫 완료 혼합계 실행은 `runs/pe-protein-pilot01`입니다. 1,000스텝(1 ps), 43,288원자, 물 11,945분자, Na⁺ 4개로 완료되었습니다. 단백질–PE 최소 중원자 거리는 0.50465 → 0.42674 nm, 0.45 nm 기준의 최종 근접 잔기는 γ Thr393 한 개였습니다. 이 잔기는 인공 절단 말단 394 바로 옆이므로 전체 피브리노겐의 선호 흡착 부위로 해석할 수 없습니다. 단백질 Cα RMSD는 최소화 구조 대비 0.931 Å, 최종 순간 온도는 301.3 K입니다. 목표 온도 도달과 평형화는 별개이며 이 계산으로 평형화를 확립하지 않았습니다.

현재 인스턴스의 기존 8000 포트 서버는 응답에 실패하여 결과 검증용 서버는 8001에서 실행했습니다. `PORT=8001 npm start`로 다른 포트를 선택할 수 있습니다. 원래 서버를 임의로 종료하지 않고 다른 포트에서 결과 JSON·CSV·PNG의 응답과 내용을 검증했습니다. 환경 설정의 설치·시작 지침 초안도 갱신했으며 게시 여부와 새 작업에서의 복원은 별도로 확인해야 합니다.

## 여섯 사슬 관측 코어 + 중성 캡 + 크기를 맞춘 PE (새 모델)

```sh
.venv/bin/python md/prepare_fibrinogen_core.py data/3ghg-upload/source.cif data/new-six-chain-core
.venv/bin/python md/relax_core_hydrogens.py data/new-six-chain-core --iterations 100
.venv/bin/python md/build_core_pe.py data/new-six-chain-core data/pe-gaff2/model02 data/new-six-chain-PE
.venv/bin/python md/report_core_model.py data/new-six-chain-core data/new-six-chain-PE artifacts/core-model
```

실제 성공 산출물은 `data/fibrinogen-core03/`, 그 아래 `local-relaxation-restrained/`, `data/fibrinogen-core-pe02/`입니다. 다른 실패·진단용 디렉터리의 좌표를 사용하지 않습니다. 원본 CIF는 보존합니다.

- 조립체 1의 α·β·γ 두 세트(관측 1,913잔기)를 함께 유지합니다. 원본 `_struct_conn`의 이황화 결합 25개(사슬 간 15개)를 정확히 반영하고 생성된 토폴로지와 대조합니다. 원본 칼슘 이온 4개를 유지합니다.
- 각 사슬의 관측된 연속 구간을 유지하고, 알려진 좌표가 없는 긴 αC 및 기타 말단 구간은 생성하지 않습니다. 부재 좌표는 단일한 임의 구조로 메울 수 없으므로 전체 길이 복원으로 주장하지 않습니다. 이 작업 모델에서는 당쇄와 knob 펩타이드를 생략합니다.
- 잘린 경계에는 ACE/NME 캡 각 6개를 붙입니다. 캡별 전하가 0인지 검사합니다. 캡이 중성이라는 뜻이며 단백질 전체가 중성이라는 뜻은 아닙니다. 코어 자체는 약 −10 e이므로 이후 용매화에서 중화가 필요합니다.
- 캡 연결은 이상적인 아미드 결합 길이·각도로 시작하고, 원본 중원자는 원래 파일과 처음에 정확히 일치하는지 검사합니다. 수소/캡을 추가한 코어는 30,215원자입니다.
- 수소·캡의 국소 사전 완화에는 HBonds 제약과 관측 중원자에 대한 100,000 kJ/mol/nm² 위치 제약을 사용합니다. 실제 최대 관측 중원자 이동은 0.03707 nm(0.3707 Å), X–H 결합 상대 오차 최대는 약 1.82e−6이었습니다. 1.2 Å 미만 비결합 원자 접촉은 18→0개로 감소했습니다. 이 단일 거리 기준은 완전한 구조 검증이 아닙니다.
- 사전 완화는 용매 없는 비주기 1 nm reaction-field 설정의 진단적 국소 준비이며, 용매화·평형화 또는 흡착 MD로 해석하지 않습니다. 완전 고정 중원자/무제약 H 시도는 비정상 짧은 접촉이 남아 실패 진단으로 표시했고 사용하지 않습니다.
- 코어는 길이 약 45 nm라 기존 γ조각용 작은 PE 패치가 부적합합니다. 새 패치는 C24H50 올리고머 15×23×2층=690사슬이며, 물 없는 혼합계는 81,275원자입니다. 기존 GAFF2/AM1-BCC 파라미터를 사용하며 이황화 결합 25개를 재검증했습니다.
- 물을 넣기 전 입력 파일: `dry.prmtop`, `dry.inpcrd`, `dry-complex.pdb`. 물 상자 추정은 약 52.8×13.0×8.9 nm, 용질 배제 전 총 원자 추정 약 69만 개입니다. 이 숫자는 실제 용매화 완료나 실행된 MD 원자 수가 아닙니다. 새 모델의 용매화·중화·압력/밀도 및 계면 평형화·흡착 MD는 아직 수행하지 않았습니다.
- 앱은 이 새 준비 모델과 이전 작은 γ절편의 실제 MD 결과를 구분해 표시합니다. 보고서: `artifacts/core-model/report.html`.

## 작은 계: 단계적 온도·밀도 준비 절차

전체 코어로 확대하기 전에 기존 43,288원자 γ절편–PE–물–Na⁺ 계에서 준비 절차를 실행합니다. 새로운 단계적 준비 구간이며 이전 1 ps 실행의 같은 앙상블 연속 계산으로 표시하지 않습니다.

```sh
OPENMM_CPU_THREADS=4 .venv/bin/python md/equilibrate_small_pe.py data/pe-protein-model04 runs/pe-protein-pilot01 runs/new-staged-preparation --npt-ps 8
.venv/bin/python md/validate_staged_preparation.py runs/new-staged-preparation
.venv/bin/python md/report_staged_preparation.py runs/new-staged-preparation artifacts/staged-preparation
```

기본 단계는 다음과 같습니다. 모든 힘장·원자 수·입력 체크섬 및 중성 전하를 먼저 확인합니다.

1. 이전 좌표에서 단백질/PE 중원자에 1000 kJ/mol/nm² 위치 제약을 걸고 최소화. 이전 속도는 이어 쓰지 않고 100 K에서 초기화.
2. NVT 100 K 0.5 ps → 200 K 0.5 ps → 310 K 1 ps. 2 fs, HBonds와 rigid TIP3P 사용. 이 구간은 부피가 고정되므로 밀도가 일정한 것이 평형의 증거가 아님.
3. 310 K, 1 bar 등방성 NPT 8 ps. 고정 좌표 제약을 모두 끄고 유한 PE 패치와 물 상자의 부피 조정을 관찰. 무한 평면 또는 막에 이 압력 조절 방식을 그대로 적용하지 않음.
4. 현재 NPT 좌표·상자를 기준으로 단백질에는 약한 위치 제약 100, PE 아래층에는 1000, 위층에는 0을 적용한 NVT 2 ps 계면 점검. 이는 준비 절차 실행 점검이지 밀도/계면 평형 확립이 아님. 앞선 NPT의 미완료 평형화를 대신하지 않음.

기본 총 구간은 12 ps(6000스텝)입니다. 길이만으로 평형화를 선언하지 않습니다. 결과에 `equilibration_established=false`, `production_ready=false`를 유지하며, 추가 NPT와 충분히 긴 추세/독립 반복 검증이 필요합니다. 단백질은 아직 인공 말단이 있는 작은 γ절편입니다. 새 여섯 사슬 모델의 결과로 일반화하지 않습니다.

각 단계의 System XML, State XML, 체크포인트, DCD, 최종 좌표와 별도 단계 보고서를 저장합니다. 검증기는 실제 스텝, 유한 수치, NVT 부피 고정, NPT 부피 반응과 고정 좌표 제약 해제, X–H 제약 오차 및 DCD 프레임 수를 검사합니다. 최종 PDB의 CRYST1은 권위 있는 저장 State의 상자에서 생성하며, 기존 실행 프로세스의 작성기 수정 전 출력도 원자 좌표를 변경하지 않고 바로잡습니다.

전체 밀도는 혼합계 전체 질량/상자 부피입니다. 벌크 물만의 밀도 측정이 아니고, 목표 압력과 순간 압력도 구분합니다. 이 구현은 1 bar를 barostat 목표로 설정하지만 순간 압력을 직접 기록하지는 않습니다. 중화용 Na⁺ 4개만 사용하고 염을 추가하지 않습니다.

첫 실제 단계적 준비 실행: `runs/pe-staged-preparation01/`. 총 6000스텝, 새로운 준비 구간 12 ps가 실행됐습니다. 끝 온도 310.29 K, 전체 밀도 0.92794 g/cm³, 부피 465.405 nm³였습니다. 초기 밀도 0.78486 g/cm³, 부피 550.253 nm³에서 변했습니다. 마지막 NPT 2 ps의 밀도 변화는 +1.63%로, 아직 충분한 평형화가 확인되지 않았습니다. 뒤의 NVT에서 밀도/부피가 평평한 것은 설정으로 고정한 결과입니다.

실행 무결성 검증에서 단계별 DCD 10/40/10프레임(총 60), 실제 기록 61개, 총 6000스텝을 확인했습니다. X–H 제약 상대 오차 최대는 각 단계 약 1.3e−6으로 검사 기준 안에 있었습니다. NVT 부피 불변, NPT 부피 반응(기록 사이 27번 변화), NPT 고정 좌표 제약 0개를 검증했습니다. 이것이 과학적 평형 검증을 대체하지 않습니다. 보고서: `artifacts/staged-preparation/report.html`.

### 추가 실행: 100 ps NPT → 구조 검사 → 계면 준비 → 큰 모델 용매화

사용자 요청으로 `runs/pe-followup01`에서 후속 작업을 실행한다. 진행 중인 작업을 완료 결과로 해석하지 말고 `pipeline.json`, 각 단계의 `progress.json`과 `result.json`을 확인한다.

- `01-npt-100ps`: 이전 `pe-staged-preparation01/03-restrained-interface-NVT-check`의 좌표·속도·상자를 유지한 100 ps NPT. 310 K, 1 bar, 2 fs, 50,000스텝, 위치 억제 없음. thermostat/barostat 난수열은 새로 시작하므로 비트 단위 체크포인트 연속 실행이 아니다. 관측 간격 0.2 ps, DCD·재시작 저장 간격 1 ps.
- `02-structure-check`: 10 ps 구간별 밀도·에너지 기술 통계, 단백질 Cα 회전반경, PE 사슬 길이·배향, 최소 이미지 물–용질 거리 및 물 산소 분포. 검사는 상관된 짧은 표본에서 평형이나 흡착을 확정하지 않는다. 심각한 기하학적 문제 시 후속 계산을 중지한다.
- `03-interface-release`: NPT 끝 상자를 고정한 10 ps NVT. 새로운 끝 좌표에 단백질 억제를 임시 적용한 뒤 100→10→0 kJ/mol/nm²로 2→2→6 ps 조정한다. 아래층 PE 탄소 지지 억제 1,000은 유지하며 위층은 자유롭다. 이미 억제가 없는 NPT와 별개의 지지 표면 준비 진단이다.
- `data/fibrinogen-core-pe-solvated01`: 검토된 큰 건조 모델의 원래 LEaP 레시피로 25개 이황화결합을 재구성하고 TIP3P 물(10 Å 패딩)과 중화 Na⁺를 추가한다. 150 mM 염 조건이 아니며 큰 모델의 최소화·MD는 이 작업에 포함하지 않는다. 실제 생성 원자 수와 원자 순서·좌표·전하·이황화결합 검사를 manifest에 기록한다.

`md/run_pe_followup_pipeline.py`가 이미 실행 중인 NPT를 감시하고 위 단계들을 순서대로 실행한다. 중복 실행 방지 잠금이 있으므로 현재 작업이 살아 있는지 확인하기 전에 다시 실행하지 않는다. 프로세스 중단 시 체크포인트·상태는 보존되지만 자동 재시작을 구현했다고 가정하지 않는다. 오류나 10분 이상 진행 중단 시 후속 단계를 멈추고 실패 상태를 기록한다.

진행·최종 후속 보고서: `artifacts/pe-followup/report.html`. 완료 후 `md/collect_results.py`로 기존 결과까지 `artifacts/collected-results/report.html`과 `results-complete.zip`에 병합하고 전체 SHA-256 검사를 수행한다. 보고서와 ZIP은 파일 링크에서 새로 다운로드해야 최신 내용을 받을 수 있다.

### 큰 모델 위치 억제 최소화

`md/minimize_large_pe.py data/fibrinogen-core-pe-solvated01 runs/core-pe-minimization01/output --iterations 200`으로 577,890원자 주기 PME 모델의 심한 겹침을 검사하고 에너지 최소화를 실행한다. 단백질 중원자와 PE 탄소 전체에 1,000 kJ/mol/nm² 위치 억제를 적용하며, Ca/Na와 물은 위치 억제에서 제외한다. 모든 원자에 실제 질량을 유지하고 HBonds 및 강체 TIP3P 제약을 사용한다. MD·가열·압력 조정은 실행하지 않는다.

`minimization.csv`의 callback은 MD 스텝이 아니다. L-BFGS 반복 번호는 내부 제약 페널티 강도를 조정할 때 다시 0으로 돌아갈 수 있다. `constraint_penalty`는 최적화기의 결합 길이 보정 항이며 단백질/PE 위치 억제와 다르다. 반복 상한은 제약 보정 시도마다 적용되며 완전한 수렴을 인증하지 않는다. 최대 30분 벽시계 예산을 초과하면 실패 상태와 점검용 중간 파일을 남기고 후속 계산을 멈춘다.

출력 `minimized.inpcrd`는 원본 `solvated.prmtop`과 짝을 이루며, 전체 구조는 원자 번호 한계를 피한 `minimized.cif`로 저장한다. 고정 상자이므로 전체 밀도는 변하지 않으며 물의 밀도·온도·압력 평형이 성립했다는 뜻이 아니다. 단독 보고서는 `artifacts/large-minimization/report.html`, 통합 보고서와 원본 묶음은 기존 `artifacts/collected-results` 경로를 갱신한다.

큰 모델 첫 최소화 `core-pe-minimization01`은 30분 벽시계 예산 안에 제약 보정이 끝나지 않아 실패로 보존했다. `core-pe-minimization02`는 그 시도의 보존된 `inspection-state.xml`을 복사해 원본 좌표를 덮어쓰지 않고 다시 시작한다. 결합 길이를 재투영하고 최적화 상한을 보정 시도마다 20회로 줄여 제약 강도 조정을 더 자주 수행한다. 두 번째 시도의 결과·수렴 한계는 별도로 판단하며 첫 시도를 성공으로 합치지 않는다. `minimize_large_pe.py`의 `--start-state` 옵션과 config의 재시작 SHA-256으로 좌표 출처를 기록한다.

### 큰 모델 1 ps 저온 MD 수치 시험

`md/heat_large_pe.py runs/core-pe-minimization02/output data/fibrinogen-core-pe-solvated01 runs/core-pe-low-temperature01/output`는 실제 577,890원자 계의 첫 MD다. 기존 위치 억제(단백질/캡 중원자·PE 탄소, 1000 kJ/mol/nm²)와 원래 억제 기준, 상자를 유지한다. 1 fs, LangevinMiddle(마찰 5/ps), seed45, HBonds 및 강체 TIP3P, 10 K 초기 속도로 10 K 200스텝 → 50 K 200스텝 → 100 K 600스텝을 실행한다. 25스텝마다 온도·에너지·제약 오차·억제 원자 이동·1913개 Cα RMSD·네이티브 단백질–PE 최소 거리를 기록하고 100스텝마다 전체 DCD 프레임과 재시작을 저장한다. 310 K 가열·밀도 조정·계면 평형화를 실행했다고 해석하지 않는다.

최소화가 완전 수렴한 것은 아니므로 초기 퍼텐셜 에너지 감소와 목표보다 높은 순간 온도가 나타날 수 있다. 실제 온도와 에너지 추세를 함께 검토하며 단지 비정상 값이 없다는 것으로 열적 평형을 선언하지 않는다. 온도>350 K, 상대 제약 오차>1e-4, 큰 억제 원자 이동>0.3 nm 또는 심한 겹침은 중단 조건이다. 실패 시 최근 검증 가능한 재시작 및 실패 기록을 검토하고 이후 단계로 자동 진행하지 않는다. 이 안전 한계는 물리적 안정성·평형의 통계 검정이 아니다.

### 큰 모델 100 K 유지 2 ps (checkpoint 계속 실행)

`md/hold_large_pe.py runs/core-pe-low-temperature01/output data/fibrinogen-core-pe-solvated01 runs/core-pe-100K-hold01/output`는 검증된 이전 checkpoint를 같은 System·최종 LangevinMiddle 적분기·CPU 4스레드로 불러온다. 좌표·속도를 이전 XML과 비교해 일치 여부를 검사하고, 새 속도 초기화·최소화·제약 재투영·난수열 재설정을 하지 않는다. 100 K, 1 fs, 마찰5/ps, 기존 위치 억제 및 상자를 그대로 유지하여 추가 2000스텝(2 ps)을 실행한다. Context는 1000→3000스텝, 1→3 ps로 계속된다. CSV의 `step`, `time_ps`는 추가 구간의 상대 값(0→2000, 0→2 ps)이며 State의 시계는 누적 값이다. OpenMM DCDReporter가 저장한 DCD 헤더는 구간 상대 스텝(100→2000)이므로 DCD 시간에 1 ps를 더하면 누적 context 시간이 된다. 원본 헤더를 수정하지 않고 교차 검사에 이 오프셋을 명시한다.

25스텝 관측 및 100스텝 저장을 유지한다. 추가 구간 시작 대비 억제 중원자 이동>0.3 nm 등 보수적인 중단 한계를 적용하며, 이 이동 한계에 걸릴 경우 거부된 상태와 원자 이동 정보를 보존해 국소 구조를 검토한다. 한계 통과와 실제 평형을 혼동하지 않는다. `validate_large_hold.py`는 XML/CSV/DCD, context시계, 속도 연속성, 동일 System 및 적분기 설정을 교차 검사한다. 0.5 ps별 평균·표준편차·기울기는 상관된 짧은 표본의 기술 통계이며 평형 검정이 아니다. 큰 모델의 총 MD는 3 ps이며 작은 γ절편 계의 100 ps와 합산하지 않는다.

## 실험 입력을 미룬 후속 준비와 환경 중단 복구

- 이전 준비 작업 `4bdbdd36b0984e948232ae503e35b10b`는 마지막 관측 6,175스텝 이후 갱신이 없었고, 현 환경에서 기록된 계산·취합 PID가 없었습니다. 수치 실패로 단정하지 않고 `interrupted`로 기록했습니다.
- 검증된 체크포인트 6,100스텝, Context 시간 6.1 ps에서 새 작업 `68bdacf68ea34dbc8f4b4c5dcc69a4bc`로 재개했습니다. pipeline 시작 3 ps 대비 3.1 ps의 관측 접두부만 이어받았습니다. 이전 6,125–6,175스텝의 관측은 원본에 보존하되 새 계보에 중복 포함하지 않습니다. 위치·속도·제약을 재설정한 새 시작이 아니라 바이너리 체크포인트 연속입니다.
- PTFE/PVDF 유한 분자 생성·후보 파라미터화, 여섯 dry 코어 배향, 미실행 반복 조건표·PMF 설계, 실제 작은 PE 궤적의 공통 기하 분석을 추가했습니다. 결정 표면·생산 MD·자유에너지 결과로 표시하지 않습니다. [후속 연구 보고서](../artifacts/research-campaign/report.html)의 현재 상태가 기준입니다.
- 큰 MD의 재개와 취합, 후속 연구 상태 추적은 별도 프로세스이며 클라우드 환경 종료를 넘어 계속된다고 보장하지 않습니다.

## 추가 NPT와 물 공간 진단

완료 작업 `5654118d6c7f48129395d353d5cf6387`의 세 NPT 블록에서 관측된 체적 변화는 각 49개 관측 간격 중 3, 3, 8번이었습니다. 전체 혼합계 밀도는 마지막에 0.8520 g/cm³였고 계속 증가했습니다. 한 프레임의 용질에서 0.8 nm 이상 떨어진 공간을 적분한 물 밀도는 약 0.8422 g/cm³(공간 적분 표준오차 약 0.0018)였습니다. 이는 MD 불확실성이나 평형 밀도 판정이 아닙니다.

새 작업 `6ce66315455240338927597da4bc6757`에서 같은 checkpoint를 이어 추가 15 ps의 NPT와 순간 압력·체적 변화 횟수를 기록합니다. 누적 14.6 ps에서 시작하며, 새 RMSD 기준은 그때의 좌표입니다. 초기 순간 분자 압력은 약 -189 bar였습니다. 이 한 값만으로 수렴이나 평균 압력을 판단하지 않습니다. 이전의 실행 완료·짧은 밀도 추세 검사 통과가 물/계면 평형을 뜻하지 않는다는 한계를 유지합니다.
