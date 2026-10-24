# 검증된 전체 원본

`results-complete.zip`은 이 브랜치에서 Git LFS로 관리하는 전체 계산 원본입니다. 작은 포인터 파일만 존재하면 아직 실제 데이터가 내려받아지지 않은 상태입니다.

```bash
git lfs install --local
git lfs pull
python scripts/restore_raw_results.py
```

원본은 759,413,919바이트이며 SHA-256은 아래와 같습니다.

```text
b59c38939187091b4c2211c503567d77b73575cc2235e6e568ef7c24f7a42e24
```

390개 원본 파일, 생성 당시 코드, 성공·실패 기록, 구조·좌표·속도·궤적·재시작·관측 데이터와 보고서를 포함합니다. ZIP 내부 `FILES.json`은 각 파일의 크기와 SHA-256 색인입니다. `docs/result-archive.json`도 같은 전체 ZIP을 식별합니다.

복원기는 `data/`, `runs/`만 검증해 복원하므로 이 브랜치의 정리된 코드를 이전 스냅샷으로 덮어쓰지 않습니다. 풀린 작업용 데이터는 계속 Git에서 제외합니다. 원본 ZIP을 다시 압축하거나 수정하면 새 데이터 버전이 되므로 기존 파일을 조용히 교체하지 마세요.

이 환경에서는 다른 작업의 Git hook을 바꾸지 않도록 LFS filter만 저장소 범위에 설정했습니다. 데이터 업로드는 `git lfs push origin codex/fibrinogen-pe-md`로 별도 수행한 뒤 브랜치를 push합니다. 일반적인 새 클론에서는 `git lfs install --local`이 설정한 pre-push hook을 사용할 수 있습니다.
