# 2026-10-29 재개: 주기 경계 생성기 수정과 citrate 수치 검사

이 작업은 `/workspace/practice-md`, Git 브랜치 `md/cbaa-surface`에서 수행했다. 기존 다른 작업과 체크포인트는 보존했다. 계산 완료·기하 검사·물성 검증을 구분한다.

## PC 전하 계산

재접속 당시 실행 중인 프로세스는 없었다. PC8은 evaluation3938, PC16은 evaluation515의 채택된 좌표와 해시가 남아 있었고, 각각 RMS gradient0.012593/0.117627kcal/mol/Å였다. 이 좌표에서 재개했다. 실제 수렴 목표0.0005는 완화하지 않았다. 최적화 이력은 다시 시작하지만 기존 좌표·로그·이전 forcefield는 지우지 않는다. 진행 중 상태는 `runs/pc-length-sqm-lbfgs-v3/pc*/status.json`과 실제 프로세스를 함께 확인해야 한다.

## 이전 periodic packing 해석 수정

**v2–v4에서 입력의 atom1만 제한하면 나머지 사슬이 상자 경계를 넘을 수 있다고 설명한 것은 잘못됐다.** 설치된 Packmol21.0.2는 PBC를 켜면 모든 자유 원자에 inside-box 제한을 자동으로 추가했다. 로그의 `PBC on: Constraint automatically added to non-fixed atoms`와 동일 동작이 있는21.2.3 소스의 `app/packmol.f90`으로 확인했다. 원래 입력의 atom1 제한은 이 자동 제한을 없애지 못했다.

낮은 밀도 PC8/16 네 조건은 기존 생성기에서 실패했다. 첫 독립 진단의 최소 사슬 간 거리는0.524–1.724Å로 기준 미달이었다. 단순 지역 최적화도 이 제한을 제거하지 못해 실패했다. 해당 구조를 MD 시작점으로 채택하지 않았다. 기존 원본·배치·로그·pose 파일을 보존하고 잘못된 제한으로 실행 중인 **배치 생성기만** 종료했다. QM 계산은 그대로 진행했다.

별도 variant는 공개 Packmol21.2.3 commit `14e50c65fa9120b58e1ba33ad482c7c7260f72b2`에서 자동 box 제한을 분자당 `iatom==1`에만 적용한다. 같은21.2.3 소스로 수정하지 않은 reference executable도 별도로 빌드했다. 소스의 의미 있는 변경은 자동 제약 배정 조건과 진단 문구이며, 주기 거리 계산·반경·최적화 루틴은 **같은21.2.3 reference 기준으로** 그대로다. 설치된21.0.2 binary나 원본 source tree는 교체하지 않았다. 버전 차이까지 단일 패치 효과라고 주장하지 않는다.

실행 검증에서12Å 상자보다 긴30Å 두 점 rigid rod는 원래/동일 버전 reference에서 실패하고 variant에서 길이 보존·경계 통과에 성공했다. 경계 양쪽의1Å 근접 원자 쌍도 주기 거리를 고려해 분리됐다. 이 toy rod는 화학 분자가 아닌 기하 회귀 검사다. 근거는 `packmol-anchor-build.json`, `packmol-anchor-validation.json`이다.

기존 두 PC4 완료 조건은 재실행하지 않았다. 나머지10조건은 새 `anchor-0001` 디렉터리에서 마지막 pose부터 이어가며, 아직 시도하지 않은 조건만 원래 seed로 시작한다. 독립적 사슬 간·자기 주기 이미지 거리 기준2.0Å 및 PDB 반올림 허용0.002Å는 그대로다. 첫 PC16 seed2027 결과는 최소 사슬 간2.00611Å,자기 이미지23.63542Å,최대 결합 길이 편차0.00204Å 미만으로 통과했다. 이는 forcefield나 기판 평형 승인과 다르다. 다른 조건은 실행 상태/독립 검사를 따로 확인한다.

기존 PC4 두 조건에는 모든 원자를 기본 상자에 가둔 초기 배치 편향이 있었다. 완료된 EM/2ps NVT/20ps NPT의 수치 증거는 유지하지만, 새로운 생성기의 PC8/16과 직접 물성 비교를 확정하지 않는다. 이전 PC4 밀도 드리프트도 남아 있다.

추가 확인에서 PC16 seed2026도 최소 사슬 간2.03709Å로 통과했다. PC8 seed2027은 수정된 생성기에서도 지역 최적화가 정체돼 보류했다. 이 경우에는 저장된 pose/seed/box/distance 기준을 유지하면서 이전에 꺼 둔 표준 molecule-relocation 기능을 다시 켜는 후속 시도를 검토·준비했다. 현재 queue 종료 후 supervisor가 이 명시적으로 검토된 작업만 이어받는다. 미검토 실패를 자동 재시도하지 않는다. 살아 있는 queue lock을 확인하여 중복 실행을 막는 검사도 통과했다.

## Citrate 추가 계산과 실제 MD 결과

기존 최종 SQM 좌표에서 citrate³⁻와 한쪽 말단이 양성자화된 hydrogen citrate²⁻를 별도 경로로 정밀 최적화했다. 실제 최종 RMS gradient는 각각0.00049798/0.00045222kcal/mol/Å로 목표를 충족했다. 채택 좌표를 배정밀도로 저장했다. 재계산한 AM1-BCC 원자 전하는 기존 후보와 **저장된 정밀도에서 모두 동일**했으며 전하 재정규화는 하지 않았다. 정수 총전하−3/−2와 카복실레이트 전하 대칭·연결을 확인했다.

상세 SQM 출력의 derivative header와 빈 `QMMM:` 진단 줄을 Antechamber가 원자로 오독하는 문제도 발견했다. 원본 출력은 그대로 보관하고 이 진단 줄만 제거한 별도 호환용 사본을 사용한다. 좌표·전하·에너지 줄은 그대로이며, 수렴 판정은 항상 원본 출력으로 한다. 실제 acetate/SQM→Antechamber 회귀 검사를 통과했고, 후속 PC8/16 변환 경로에도 적용했다.

| 후보 | EM 최대힘(kJ/mol/nm) | MD 길이 | 에너지/기하 프레임 | 평균 온도(K) | 전체 결합 길이 범위(nm) |
|---|---:|---:|---:|---:|---|
| Citrate³⁻ + Na3 | 778.35 | 2ps | 41/41 | 296.12 | 0.09724–0.16709 |
| Hydrogen citrate²⁻ + Na2 | 762.29 | 2ps | 41/41 | 295.00 | 0.09724–0.16790 |

둘 다 GAFF2/AM1-BCC,Amber/TIP3P,JC Na의 중성3nm 물 상자에서 수행했다. LEaP의 두 단일 음이온 비중성 경고는 각각 정수 전하에 대한 것으로 검토하고 MD 상자는 Na로 중성화했다. GROMACS 경고 우회 없이 전처리, EM, MD 및 유한 에너지·결합 거리·LINCS 검사를 통과했다. 이 상자에서는 용질 명목 농도가약61.5mM로 SOP의 총 citrate약10mM와 다르다. 짧은 수치 검사이며 수화/이온쌍 평형이나 실험 농도 재현을 주장하지 않는다.

작은 후보 topology는 `models/buffer-parameters/refined-citrate-candidates/`, 결과 근거는 `citrate-gradient-refinement-audit.json`과 `refined-citrate-water-audit.json`, 원본 계산은 `runs/citrate-*`에 보관한다. 원래 citrate와 phosphate/azide/PC/FXII 결과는 변경하지 않았다.

## 남은 작업과 재현

PC8/16 최종 기울기·전하/힘장 검토와 나머지 배치/기판 수치 검사를 마쳐야 한다. 버퍼는 citrate의 terminal/central 양성자화 미시적 상태, pH 종 분포와 수용액 상호작용 검증이 남는다. 전체 코팅/FXII/cPBSzI 시스템은 아직 시작하지 않았다. 현재 작업에 추가 사용자 입력은 필수가 아니다.

새 환경의 기하 도구는 기본 설치 후 `bash scripts/install_packmol_anchor.sh`로 준비한다. Debian 서명/체크섬으로 확인한 compiler 패키지를 로컬 디렉터리에 추출하며 system install을 하지 않는다. 원본 source revision이 다르면 덮어쓰지 않고 멈춘다. 실행 중인 worker와 lock을 먼저 확인하고, 기존 `complete_pc_screen.py --watch`와 `run_pc_packing_margin.py run`을 재개한다. `prepare_packmol_anchor_recovery.py`는 최초 준비 전용이며 기존 attempt에 덮어쓰지 않는다.
