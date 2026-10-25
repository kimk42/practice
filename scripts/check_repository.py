"""Reject raw Git payloads; allow the verified originals only as an LFS pointer."""
from pathlib import Path
import re
import subprocess

root=Path(__file__).resolve().parents[1]
tracked=subprocess.check_output(['git','ls-files','-z'],cwd=root).decode().split('\0')
problems=[]
for name in filter(None,tracked):
    size=int(subprocess.check_output(['git','cat-file','-s',':'+name],cwd=root))
    if name=='datasets/results-complete.zip' or re.fullmatch(r'datasets/results-delta-[0-9]{4}-[0-9a-f]{12}\.zip',name):
        pointer=subprocess.check_output(['git','show',':'+name],cwd=root) if size<1024 else b''
        attr=subprocess.check_output(['git','check-attr','--cached','filter','--',name],cwd=root).decode().strip()
        if not re.fullmatch(rb'version https://git-lfs.github.com/spec/v1\noid sha256:[0-9a-f]{64}\nsize [0-9]+\n',pointer) or not attr.endswith('filter: lfs'):
            problems.append('Original archive must be an LFS pointer, never a raw Git ZIP.')
        continue
    if size>10*1024**2:problems.append(f'Git object over 10 MiB: {name}')
    if name.startswith(('runs/','.venv/','.tools/')) or Path(name).suffix in ('.dcd','.chk','.zip'):
        problems.append(f'Raw/local payload in Git: {name}')
if problems:raise SystemExit('\n'.join(problems))
print('Git object size, raw-data exclusion and LFS pointer checks passed.')
