# 설계 검토 후 다섯 단계 실행 결과

이 기록은 `DESIGN_REVIEW_2026-10-28.md` 이후 사용자가 승인한 다섯 단계의 실제 결과다. 기존 계산의 조건을 바꾸거나 결과를 새 모델 결과로 재표시하지 않았다.

| 단계 | 수행 결과 | 아직 주장하지 않는 것 |
|---|---|---|
| 1. 기존 표면 진단·재실행 기준 | 10ps/100ps 밀도 CSV를 같은 규칙으로 요약하고, 변경별 재분석/새 계산 규칙과 표면 비교 설계를 작성 | 밀도 프로파일만으로 평형·실제 PC 물성 입증 |
| 2. PC 사슬 길이 후보 | PC8(264원자),PC16(528원자) 연결·중성·분자식·말단·기하 검사 및 MMFF 기하 완화 완료 | 새 PC 힘장/벌크 기판/코팅 평형화 완료 |
| 3. 버퍼 수량 기준 | 상자 전체 부피 기준을 역사적 표로 보존하고 물 분자 수 기준의 잠정 계산표로 변경 | 측정된 접근 가능 부피,완성된 pH7.4 혼합물 |
| 4. 인산염 후보 | 원본 보존,정수 전하·대칭 제약 보정;GAFF2 topology 생성;각각 중성 Na/물 상자에서 EM+2ps와41프레임 기하/에너지 검사 통과 | 수화/다가 음이온 응집 정확성 또는 최종 cPBSzI 승인 |
| 5. 아지드 후속 | 기존2ps checkpoint/동일 TPR 조건에서100ps 연장,총102ps;2001프레임·RDF·Na 거리 검사 완료 | 이온쌍 평형,실험 희석농도 재현,단백질 흡착 검증 |

## 정량 결과

기존 소형 표면100ps의 PC 질량5–95%를 포함하는 z 구간은2.7–5.3nm이며,그 구간의 xy/시간 평균 PC 밀도631.59kg/m³와 물421.67kg/m³다. 이는 공간 평균이며 국소 벌크 밀도나 공극률 측정은 아니다. PC4의 OH말단/BPA단위 비0.5를 PC8에서0.25,PC16에서0.125로 낮춘 구조를 준비했다. 새로운 매개변수는 아직 배정하지 않았다. `models/surface-revision-v2/`와 `config/surface-revision-design.json` 참조.

물 기준은 현재 기준계의 물 개수와 **55.5M라는 명시적 근사 가정**을 사용한다. Coated 기준계112801개 물에서 NaCl/NaI/citrate/azide 수는224/20/20/6, bare113579개 물에서는225/20/20/6이다. 수화물에 따른 phosphate 대안도 따로 유지했다. 이 수량은 최종 용질 삽입·재용매화 후 물 개수가 바뀌면 다시 계산해야 하며,완화 후 지정한 bulk-water 영역 농도를 확인해야 한다. 과거 상자전체 기준258/23/23/7을 그대로 적용하지 않는다. 기존 중화용 Cl6에 새 버퍼를 중복 추가하지 않도록 inventory 정책을 명시했다. `config/cpbszi-system-preparation.json` 참조.

H₂PO₄⁻/HPO₄²⁻의 최대 원자별 전하 보정은 각각0.000142857143e/0.000166666667e다. 좌표·원자형·연결·원래 QM/BCC 출력은 보존했다. 추가 소수 자릿수는 산술 정밀도이며 복원된 QM 결과가 아니다. 단일 음이온의 LEaP 비중성 경고2건은 기대한 정수 전하에 대한 것으로 검토했고,각 MD 상자는 Na1/Na2로 중성화했다. GROMACS 전처리의 경고 우회는 사용하지 않았다. 원본 ACPYPE6자리 전하 출력이 재반올림되므로 이름/순서를 확인한 뒤 별도 GROMACS 후보에 projected mol2 전하를 그대로 기록했다.

두 인산염 EM 최대힘은848.49/836.34kJ/mol/nm,2ps의41개 에너지 프레임은 유한했다. 정의된 O–H/P–O 결합의 전체 저장 프레임 범위는 각각0.09725–0.16034nm/0.09725–0.16343nm였다. `models/buffer-parameters/projected-phosphate-water-audit.json` 참조.

아지드100ps 연장 구간의 평균 온도300.048K,최대 NN거리 오차2.23×10⁻⁷nm다. 앞/뒤50ps의 terminal-N–water-O RDF 최대위치는 지정한0.20–0.45nm 탐색 구간에서 모두0.27nm,0.35nm까지의 물O 배위 적분은약4.01/4.02였다. 두 창은 독립 반복이 아니며,이 유사성으로 수렴을 단정하지 않는다. Terminal-N–Na 최단거리는0.605–2.496nm이고0.35nm 이하 저장 프레임은 없었다. 따라서 이번 궤적은 이온쌍 형성/해리 평형을 평가하지 못한다. 이 작은 상자의 명목상 각 이온 농도는약61.5mM로 실험 azide 약3.08mM와 다르다. `models/buffer-parameters/azide-solution-102ps-audit.json` 및 `models/analysis/azide-*-rdf.csv` 참조.

## 재계산 규칙과 보존

`config/model-revision-policy.json`에 규칙을 저장했다. 분석 정의 변경은 저장 궤적 재분석으로 처리하고,힘장·조성·구조 변경은 영향을 받는 계에 별도 버전·출력 경로를 만든다. 온도/압력 변경은 새 조건의 재평형화를 요구하며,기하/힘이 바뀌지 않으면 EM을 기계적으로 반복할 필요는 없다. 동일 조건의 시간 연장은 checkpoint continuation을 사용한다. 이전 조건의 통계와 새 조건의 통계를 합쳐서는 안 된다.

이번 작업 전38개 checkpoint의 SHA256가 모두 동일함을 확인했다. 원래 phosphate mol2 두 개도 변경되지 않았다. `models/validation/five-step-revision-preservation.json`에 검증했다. 원래 worktree의 기존 삭제 파일5개는 그대로 보존했다.

## 재현 경로

설치 환경 활성화 후 `prepare_surface_revision_candidates.py`, `prepare_phosphate_charge_projection.py`, `prepare_cpbszi_manifest.py`, `audit_buffer_charge_balance.py`가 준비/기록 단계다. 기존 생성 디렉터리를 보존하기 위해 구조/전하 생성기는 존재하는 출력에 덮어쓰지 않는다.

인산염의 매개변수 변환은 각 분리된 작업 폴더에서 projected mol2를 `input.mol2`로 복사한 뒤 아래 명령을 사용했다. `leap.in`은 `source leaprc.gaff2`, `loadamberparams candidate.frcmod`, `M = loadmol2 input.mol2`, `check M`, `saveamberparm M candidate.prmtop candidate.inpcrd`, `quit` 순서다. PH2는 basename을 PH2로 바꾼다. 원래 출력·명령·로그는 `runs/phosphate-projected-parameterization-v1/`에 보관했다. 신규 QM 계산은 하지 않았다.

```bash
parmchk2 -i input.mol2 -f mol2 -o candidate.frcmod -s gaff2
tleap -f leap.in
acpype -p candidate.prmtop -x candidate.inpcrd -b PH1 -a gaff2
```

`run_projected_phosphate_controls.py`와 `audit_phosphate_control_geometry.py`는 분리된 중성 물 제어계 생성 및 수치/기하 검사를 담당한다. `continue_azide_solution.py`는 검증된 source TPR/checkpoint만 연장하며 중복 worker를 파일 lock으로 차단한다. `analyze_azide_solution.py`는 결과를 새 분석 폴더로 읽는다. 완성된 분석이나 계산을 다시 실행해 덮어쓰지 않는다.

다음 과학적 단계는 긴 PC의 일관된 내부/말단 매개변수 검토와 기판 패킹·완화,그리고 citrate/phosphate/azide의 용액 상호작용 및 pH 종 분포를 명시한 후보 조성 평가다. 새로운 전체 FXII/cPBSzI 계산은 아직 시작하지 않았다. 현재 단계에 추가 사용자 입력은 필수적이지 않다.
