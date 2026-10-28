# PC 사슬 길이 검사 재개와 실제 완료 결과

작업 위치는 `/workspace/practice-md`, 브랜치는 `md/cbaa-surface`다. 원래 `/workspace/practice` 작업 내용은 수정하지 않았다. 이 기록은 코팅/FXII/cPBSzI 준비 중 **건조 PC 기판 후보 검사**의 결과다.

## 완료된 결과

환경이 다시 연결된 뒤 실제 프로세스를 검사했을 때 계산 프로세스는 없었다. 그러나 연결이 끊어진 동안 PC4 AM1 최적화, AM1-BCC/GAFF2 변환과 두 개의 건조 PC4 수치 검사가 완료된 파일이 남아 있었다. 결과 JSON만 믿지 않고 에너지 파일을 다시 읽고 체크포인트 SHA256을 확인했다.

| 조건 | EM 최대 힘(kJ/mol/nm) | NVT/NPT | NPT 평균 온도(K) | NPT 처음/마지막5ps 평균 밀도(kg/m³) |
|---|---:|---|---:|---|
| PC4,초기0.6g/cm³,seed2026 | 695.65 | 2ps/20ps | 300.75 | 604.76 / 643.96 |
| PC4,초기0.6g/cm³,seed2027 | 542.36 | 2ps/20ps | 300.48 | 604.83 / 638.53 |

각각6336원자·48사슬이며, 두 경우 모두 NPT 에너지401프레임이 유한하고 LINCS 경고가 없었다. 밀도 창 평균은6.48%/5.57% 증가했다. 따라서 **짧은 수치 검사는 완료됐지만 기판 평형화는 입증되지 않았다.** 두 packing/velocity seed는 같은 단일 사슬 conformer를 사용한다. 300K/1bar와 barostat compressibility4.5e-5bar⁻¹는 계산 프로토콜 가정이며 측정된 PC 물성이 아니다. 건조 균질 벌크이므로 이 단계의 분산 장거리 보정은 EnerPres를 사용한다.

PC4 native SQM 최종 출력의 RMS gradient는0.0003kcal/mol/Å이며 출력 반올림 상한0.00035가 목표0.0005보다 작다. BCC 원시 전하 합−0.005e는 출력된 Mulliken 전하 합과 일치한다. ACPYPE의 잔여 전하 분배를 재구성하여 기록했으며, 이 관례는 제한된 수치 검사에만 채택했다. 기존 PC4 힘장과 새 후보는 별도로 보존한다. 물성 정확성을 승인한 것은 아니다.

근거는 `models/surface-revision-v2/pc4-dry-controls-interim-audit.json` 및 `runs/pc-length-forcefield-v2/pc4/parameter-audit.json`이다.

## 진행 중인 계산과 중단 복구

PC8과 PC16은 native SQM의 완료된 첫 segment에 각각50/10회의 좌표 최적화 기록이 남아 있었다. 미완료 두 번째 segment와 원래 v1 계산을 보존하고, 완료된 좌표에서 AM1 계산을 이어갔다. `config/pc-qm-continuation.json`이 현재 경로를 선택한다.

새 경로는 SQM AM1의 실제 에너지·해석적 기울기를 L-BFGS에 전달하고 채택된 좌표를 배정밀도 JSON/NumPy 파일로 매번 저장한다. 물과 메탄올에서 기울기 부호·단위를 유한차분으로 확인했고, 물 최적화는 native SQM 최소 에너지와1e-4kcal/mol 이내로 일치했다. SQM의 단순 완료 문구만으로 수렴을 인정하지 않는다. 최종 실제 RMS gradient≤0.0005kcal/mol/Å를 요구하며 원자 순서·좌표·결합 길이도 검사한다. 재개 시 좌표는 유지하지만 최적화 이력과 SCF 상태는 다시 시작하므로 비트 단위 재시작은 아니다.

PC16의 첫 unscaled L-BFGS 시험점에서는 원자301의 기울기가 SQM 출력 필드 범위를 넘었다. 이 시험점을 채택하지 않았고, 누락된 숫자를 추정하지 않았다. 마지막 정상 좌표에서 optimizer 좌표 길이 척도를0.125로 줄여 재개했다. 에너지는 E(0.125y), optimizer 기울기는0.125·dE/dx로 일관되게 변환한다. 물리적 에너지·전하 방법과 수렴 기준은 그대로다. 작은 계에서 이 변경 후 같은 최소점과 실제 기울기 기준을 재검증했다. 향후 계산/기울기 오류도 자동 우회하지 않는다. 원본과 검토는 `pc16-optimizer-step-review.json`에 연결된다.

Packing v3의 precision0.0001은 Packmol 내부 목적함수 정지 바닥과 충돌했다. v4는 generator 목표2.004Å와 precision0.01을 사용한다. 독립 검사 기준은 기존2.0Å 및 PDB 출력 반올림 허용0.002Å 그대로다. 성공한 두 PC4 조건의 최소 사슬 간 거리는2.00150/2.00166Å다. PC8 중단 작업은 마지막 rigid-body pose를 새 `resume-*` 디렉터리에 복사하여 이어간다. 이전 packed.pdb,입력,로그,pose 파일은 덮어쓰지 않는다. 미착수 조건은 별도 대기하며 기하 검사를 통과해야 MD로 넘어간다.

## 자동 진행 범위와 보존

`complete_pc_screen.py --watch`는 PC4/8/16 전하 계산, 수렴/힘장 검토, 통과한 배치의 EM→2ps NVT→20ps NPT 순서를 감시한다. 총12조건 중 이 기록 시점에는2조건이 완료됐다. QM과 packing은 각각 최대2개, 건조 MD는 순차 실행한다. 공유 lock이 중복 계산을 막고 실패는 보류 상태에 남긴다. 기판 기하·힘장·수치 검사가 완료되기 전 전체 코팅/FXII/cPBSzI로 확장하지 않는다.

현재 상태는 `runs/pc-completion-status.json`, `runs/pc-length-sqm-lbfgs-v3/pc*/status.json`, `runs/pc-dry-packing-v4/*/status.json`과 실제 프로세스를 함께 확인해야 한다. 프로세스가 사라져도 파일의 running 문자열은 남을 수 있다. 환경 재시작 시 README의 두 worker를 기존 디렉터리에 재개하며, 오류가 있는 후보는 원인 검토 없이 초기화하거나 재시도하지 않는다.

재개 전 기존 MD 체크포인트43개 및 새 PC4 체크포인트6개를 확인했다. 전체 로컬 백업 `cbaa-surface-20261028T191757Z`는 새 계산보다 오래된 시점이다. 실행 중인 원본 계산을 멈추고 전체 백업을 만들지는 않는다. 코드·설정·작은 검토 결과는 기존 GitHub 브랜치에 저장하고, 원본 궤적·체크포인트·논문은 계속 로컬에 보관한다.

## 코팅/FXII/cPBSzI로 이어지는 과학적 검증

cPBSzI 조성은 업로드된 SOP에 이미 있다. 현재 차단점은 문서 부재가 아니다. 기판 밀도/사슬 길이/말단 효과와 완화, 이후 코팅 수화·접합 구조를 검토해야 한다. 버퍼에서는 citrate/phosphate/azide 수화·이온쌍 거동과 pH7.4 종 분포 및 정수 이온 수를 일관되게 준비해야 한다. 인산염 수화물 미상은 명시적 농도 민감도로 유지할 수 있다. 이번 PC 준비에 추가 사용자 입력은 필수가 아니다.

최종 전체 시스템은 아직 만들거나 실행하지 않았으며, 기존 코팅100ps·FXII baseline 결과를 새 PC 또는 cPBSzI 결과로 바꾸어 해석하지 않는다.

진행 중 독립적으로 수행한 citrate 표기 검토는 `CITRATE_REPRESENTATION_REVIEW_2026-10-28.md`에 기록했다. 기존 카복실레이트 표기와 전하 대칭은 확인됐지만, hydrogen citrate의 terminal/central 양성자화 분율과 수용액 검증은 아직 미정이다.
