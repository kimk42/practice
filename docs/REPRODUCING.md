# 재현과 데이터 복원

## 코드와 실제 계산 기록

이 브랜치는 시뮬레이션 프로젝트의 코드·환경·작은 원본·보고서를 분리한 스냅샷입니다. 원래 작업 디렉터리는 그대로 유지했습니다. 원본 전체 결과 ZIP에는 계산 당시 코드와 모든 기록 데이터, 실패 시도, 구조·상태·궤적이 들어 있습니다. 이 브랜치는 경로 이식성과 데이터 복원·화면 표시를 정리했으므로 ZIP의 코드 파일 체크섬이 현재 정리된 코드와 반드시 같지는 않습니다. 물리 입력·완료 결과의 원본은 바꾸지 않았습니다.

## 복원

전체 원본은 이 브랜치의 `datasets/results-complete.zip`에 Git LFS로 포함되어 있습니다. 일반 Git에는 작은 포인터가, LFS에는 실제 ZIP이 저장됩니다.

```bash
git lfs install --local
git lfs pull
python scripts/restore_raw_results.py --verify-only
python scripts/restore_raw_results.py
```

기본 복원 위치는 이 프로젝트 루트입니다. `--destination`으로 별도 디렉터리를 선택할 수 있습니다. 데이터 경로와 크기·SHA-256을 검사하며, 소스 코드·Git 메타데이터·README는 복원 대상으로 삼지 않습니다. 이미 있는 파일이 다르면 중단합니다. 중단 전에 복원한 파일은 각각 검증된 것이며 재실행은 같은 파일을 확인하고 건너뜁니다.

입력은 `data/3ghg-upload/source.cif`와 `rcsb_pdb_3GHG.fasta`입니다. PDBFixer 또는 재매개변수화가 원본 좌표를 조용히 바꾸었다고 가정하지 않습니다. 원본 파일과 수정 정책은 manifest로 추적합니다. PE 파라미터는 `model02`를 사용합니다.

## 완료 결과 독립 검사

패키지 설치와 원본 복원 후 프로젝트 루트에서 실행합니다. 이 명령은 긴 MD를 다시 실행하지 않고 저장 출력만 검사합니다.

```bash
.venv/bin/python md/validate_large_hold.py runs/core-pe-100K-hold01/output
```

현재 체크아웃의 `runs/`를 우선 찾도록 경로 조회를 정리했습니다. JSON 안 과거 `/workspace/practice/...` 경로는 출처 기록으로 보존합니다. 체크포인트는 OpenMM 버전·플랫폼·장비에 의존할 수 있습니다. 동일한 설정에서도 다른 환경의 비트 단위 궤적 일치를 보장하지 않습니다.

## 시간과 상태

큰 모델의 첫 저온 실행은 1 ps, 100 K 유지는 추가 2 ps입니다. 후자는 이전 checkpoint의 좌표·속도·상자·난수 상태를 이어받았으며 새 속도나 제약 재투영을 적용하지 않았습니다. Context는 1→3 ps, 1000→3000스텝입니다. 유지 CSV와 기본 DCD 헤더는 구간 상대 값이며, DCD 시간에 1 ps를 더하면 누적 context 시간에 대응합니다.

유지 시작 기준 RMSD와 첫 저온 시작 기준 RMSD를 그대로 이어 붙이지 않습니다. 같은 정의인 온도·퍼텐셜 에너지의 연결 자료만 `combined-temperature-energy.csv`로 제공합니다. 고정 상자에서 밀도가 평평한 것은 압력·밀도 평형의 증거가 아닙니다. 구간 통계는 상관된 표본의 기술 통계입니다.

## 보고서 재생성

원본을 복원한 뒤 실행합니다.

```bash
.venv/bin/python md/report_large_hold.py runs/core-pe-100K-hold01/output artifacts/large-hold
.venv/bin/python md/collect_results.py
```

통합 수집기는 임시 ZIP을 완성한 뒤 교체하므로 기존 정상 ZIP을 생성 도중 지우지 않습니다. 새 ZIP은 코드 정리 후의 별도 내보내기이며 원본 제공 ZIP과 구분합니다. SHA-256 색인과 입력 데이터 출처를 함께 보존합니다.

## 이어서 계산할 때

현재 저온 유지에서 수치 한계는 통과했지만 온도·에너지 추세와 국소 구조를 검토한 후 다음 가열 단계를 설계해야 합니다. 310 K, 압력·밀도 조정, 위치 억제 해제, 반복 배향·시드, 흡착 자유에너지는 별도 작업입니다. 지금 결과로 PTFE/PVDF 우열이나 전체 길이·당쇄를 가진 피브리노겐의 결합을 확정하지 않습니다.
