"""Fail if raw MD payloads or oversized files are added to Git."""
from pathlib import Path
import subprocess

root=Path(__file__).resolve().parents[1]
tracked=subprocess.check_output(['git','ls-files','-z'],cwd=root).decode().split('\0')
problems=[]
for name in filter(None,tracked):
    path=root/name
    if not path.is_file():continue
    if path.stat().st_size>10*1024**2:problems.append(f'Over 10 MiB: {name}')
    if name.startswith(('runs/','.venv/','.tools/')) or path.suffix in ('.dcd','.chk','.zip'):
        problems.append(f'Raw/local payload in Git: {name}')
if problems:raise SystemExit('\n'.join(problems))
print('Tracked-file size and raw-data checks passed.')
