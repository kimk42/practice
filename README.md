# 피브리노겐–폴리머 시뮬레이션

3GHG 관측 구조에서 출발해 PE 표면 입력을 만들고, 준비 MD와 그 한계를 기록하는 연구용 프로토타입입니다. 현재 실제 MD는 **PE**만 수행했습니다. PTFE 및 PVDF α·β·γ는 설정 항목이며 실제 비교 계산 결과는 없습니다.

## 지금까지 확인한 결과

- 작은 γ절편 혼합계: 43,288원자, 100 ps 추가 NPT와 10 ps 계면 준비.
- 여섯 사슬 관측 코어 혼합계: **577,890원자**, 물 165,535분자·Na⁺ 10개·Ca²⁺ 4개.
- 큰 모델: 초기 최소화, 1 ps 저온 시험, checkpoint에서 이어받은 **추가 2 ps 100 K 유지**까지 실행. 큰 모델 MD 누적은 3 ps입니다.
- 100 K 유지 끝 실제 온도 약 **101.02 K**, 유지 시작 대비 Cα RMSD 약 **0.170 Å**. 에너지는 계속 완화 중입니다.

이것은 위치 억제를 건 짧은 수치 시험입니다. 310 K 가열·압력/밀도 조정·계면 평형·결합/흡착 자유에너지를 검증하지 않았습니다. 긴 결손 말단·당쇄가 생략된 관측 코어이며 전체 길이의 피브리노겐 모델이 아닙니다.

## 결과부터 보기

- [최신 100 K 유지 보고서](artifacts/large-hold/report.html)
- [전체 통합 보고서](artifacts/collected-results/report.html)
- [데이터·단위·파일 설명](artifacts/data-guide.md)
- [실행 이력과 이전 명령](docs/EXPERIMENT_LOG.md)
- [재현·대용량 데이터 복원](docs/REPRODUCING.md)

HTML은 그림을 포함한 독립 파일입니다. GitHub에서 소스가 보이면 **원본 다운로드 후 브라우저로 열면** 됩니다. 보고서와 작은 입력은 일반 Git으로, 전체 좌표·궤적·재시작 파일이 든 `datasets/results-complete.zip`은 **이 브랜치의 Git LFS**로 관리합니다. `artifacts/collected-results/FILES.json`은 그 원본 ZIP의 체크섬 색인이며 현재 Git 코드의 체크섬 목록이 아닙니다. 코드 정리 전 실제 계산에 쓴 스크립트도 원본 ZIP에 보존되어 있습니다.

## 화면 실행

Python 3와 Node.js가 있는 환경에서 프로젝트 디렉터리 안에서 실행합니다. Node 패키지 설치는 필요 없습니다.

```bash
npm start
```

자신의 PC에서 실행했다면 브라우저에 `http://localhost:8000`을 입력합니다. `PORT` 환경변수로 포트를 바꿀 수 있습니다. 화면은 완료 결과를 보여주고 구조 배치를 살펴보거나 조건 JSON을 내보냅니다. 화면의 표면 선택이 실제 힘장을 바꾸거나 MD를 자동 실행하지는 않습니다.

## 계산 환경과 검사

Python 3.12 및 Linux x86-64에서 사용한 패키지를 기록했습니다.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m unittest discover -s tests -p 'test_*.py'
npm test
```

실제 모델 구축에는 AmberTools 24.8이 추가로 필요합니다. 설치기는 홈 디렉터리에 환경을 등록하지 않고 프로젝트의 `.tools/`에 설치합니다. 이미 설치했다면 `AMBERHOME`을 그 경로로 지정하거나 각 명령의 `--amberhome`을 사용합니다.

```bash
.venv/bin/python md/install_ambertools.py
```

CI는 작은 단위·엔진 검사만 실행합니다. 큰 MD, AmberTools 설치나 수 시간짜리 재계산은 자동 실행하지 않습니다.

## 대용량 결과 보존

`runs/`와 준비된 큰 `data/` 모델, 풀린 checkpoint·DCD는 일반 Git에 넣지 않고, 전체 원본 ZIP만 Git LFS로 보관합니다. 전체 원본 ZIP은 이 브랜치의 Git LFS에 포함했습니다. 먼저 LFS 데이터를 내려받은 뒤 작업용 디렉터리로 복원합니다. 채팅의 ZIP 링크도 같은 원본입니다.

```bash
git lfs install --local
git lfs pull
python scripts/restore_raw_results.py
```

복원기는 ZIP 내부 `FILES.json`의 SHA-256과 크기를 확인하며 **`data/`, `runs/`만** 복원합니다. README나 코드·설정을 덮어쓰지 않고, 내용이 다른 기존 데이터도 덮어쓰지 않습니다. 과거 환경의 절대 경로는 기록으로 유지하고, 코드에서는 현재 체크아웃의 대응 데이터 경로를 먼저 찾습니다.

새 실행은 새 출력 디렉터리에 저장하세요. 기존 실행을 재시도할 때에는 실패 기록·원본 입력·재시작 출처를 보존합니다.

## 디렉터리

| 위치 | 내용 |
|---|---|
| `md/` | 구조 준비, PE 파라미터화, MD, 검증 및 보고서 생성 |
| `data/3ghg-upload/` | 제공받은 원본 CIF·FASTA와 구조 점검 |
| `data/pe-gaff2/model02/` | 보정 기록을 보존한 작은 PE 파라미터 번들 |
| `artifacts/` | Git에 보관하는 작은 그래프·CSV·JSON·HTML |
| `datasets/` | Git LFS로 관리하는 검증된 전체 원본 ZIP |
| `runs/` | 원본 ZIP에서 복원한 작업용 계산 데이터; Git 제외 |
| `docs/` | 실행 이력, 재현 방법, 결과 단계 설명 |
| `scripts/` | 검증된 원본 복원과 Git 용량 검사 |
| `environment/` | AmberTools Linux 패키지 고정 목록 |
| `tests/` | 작은 물리·형상·경로·복원 검사 |

두 입자 시험은 계산 엔진 검사이며 폴리머 연구 결과가 아닙니다. 총 퍼텐셜 에너지 감소나 짧은 거리·접촉을 결합 세기 또는 평형 흡착으로 해석하지 않습니다.
