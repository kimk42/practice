# 원본 보관 색인

현재 `datasets/results-complete.zip`은 누적 44.6 ps까지의 기존 검증 스냅샷이며 이후 덮어쓰지 않습니다. 대용량 ZIP 한 개를 계속 키우는 대신, 새로 생기거나 내용이 바뀐 파일을 `results-delta-0001-<hash>.zip` 형태의 Git LFS 원본으로 추가합니다. 각 증분 ZIP의 비압축 파일 합계는 1.25 GB 이하로 나눕니다. 개별 큰 파일이 이 예산을 넘으면 색인을 변경하지 않고 실패를 기록합니다.

`datasets/results-index.json`은 최신 스냅샷의 각 파일 경로·크기·SHA-256과 어느 ZIP에서 읽을지를 지정합니다. 기존 ZIP과 증분 ZIP 자체의 크기·SHA-256도 기록합니다. 이전 버전은 불변 ZIP들에 그대로 남고 색인은 최신 버전을 선택합니다. 이는 분할 압축 조각을 이어붙이는 방식이 아니라, 각자 열 수 있는 독립 ZIP과 검증된 파일 색인입니다.

## 복원

```bash
git lfs pull
python scripts/restore_raw_results.py
```

기본 복원기는 색인이 있으면 자동으로 사용합니다. 모든 ZIP의 체크섬을 먼저 확인한 뒤, 최신 파일 중 `data/`, `runs/`만 복원합니다. 코드·README를 덮어쓰지 않고, 내용이 다른 기존 데이터도 덮어쓰지 않습니다. 다른 시점의 데이터가 이미 있다면 새 목적지에 복원하세요.

```bash
python scripts/restore_raw_results.py datasets/results-index.json --destination /path/to/fresh-checkout
```

기존 `results-complete.zip`을 명시해서 복원하면 그 불변 기준 스냅샷까지만 복원합니다. 최신 전체 상태에는 색인을 사용해야 합니다. 증분 ZIP 하나만 풀어서는 전체 모델이 되지 않습니다.

## 내보내기와 무결성

통합 수집기는 `artifacts/collected-results/results-complete.zip`에 로컬 전체 스냅샷을 만들고 내부 `FILES.json`을 기록합니다. 이 로컬 전체 ZIP은 Git에 넣지 않습니다. 취합기는 그 전체 내용을 검증한 다음, 이전 색인과 경로·크기·해시가 다른 파일만 새 LFS ZIP에 씁니다. 새 파일들의 해시를 다시 확인하고 색인을 마지막에 교체합니다. 실패하더라도 기존 기준 ZIP이나 기존 색인을 덮어쓰지 않습니다.

`docs/result-archive.json`은 취합 후 색인 경로·색인 해시·원본 ZIP 목록과 로컬 전체 스냅샷 해시를 기록합니다. Git에 있는 코드는 현재 코드이며, 원본 ZIP 안의 생산 코드와 반드시 같은 버전은 아닙니다. 각 계산의 생산 코드 사본과 조건·체크포인트 출처를 함께 확인하세요. 완료된 실행별 보고서도 그 실행의 `analysis/`, `contact-review/`에 보존합니다.
