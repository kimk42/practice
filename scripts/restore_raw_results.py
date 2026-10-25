"""Restore verified data/ and runs/ entries, preserving this checkout's source code."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tempfile
import zipfile


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk := f.read(8 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def restore(archive, destination, verify_only=False, entries_override=None):
    if Path(archive).suffix==".json":
        return restore_index(archive,destination,verify_only)
    with Path(archive).open('rb') as source:
        if source.read(128).startswith(b'version https://git-lfs.github.com/spec/v1'):
            raise ValueError('The archive is a Git LFS pointer. Run git lfs install --local and git lfs pull first.')
    destination = Path(destination).resolve()
    counts = {'verified': 0, 'restored': 0, 'already_present': 0}
    with zipfile.ZipFile(archive) as z:
        entries = json.loads(z.read('FILES.json'))['files'] if entries_override is None else entries_override
        selected = []
        seen = set()
        for entry in entries:
            path = PurePosixPath(entry['path'])
            if path.is_absolute() or '..' in path.parts or '\\' in str(path) or not path.parts:
                raise ValueError(f'Unsafe archive path: {entry["path"]}')
            if str(path) in seen:
                raise ValueError(f'Duplicate manifest path: {path}')
            seen.add(str(path))
            if path.parts[0] not in ('data', 'runs'):
                continue
            target = destination.joinpath(*path.parts)
            if not target.resolve().is_relative_to(destination):
                raise ValueError(f'Archive target escapes destination: {path}')
            if target.exists() and (not target.is_file() or digest(target) != entry['sha256']):
                raise ValueError(f'Existing file differs; refusing overwrite: {path}')
            selected.append((entry, target))
        for entry, target in selected:
            if target.exists():
                counts['already_present'] += 1
            h = hashlib.sha256()
            size = 0
            temporary = None
            try:
                if not verify_only and not target.exists():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    temporary = tempfile.NamedTemporaryFile(dir=target.parent, suffix='.part', delete=False)
                with z.open(entry['path']) as source:
                    while chunk := source.read(8 * 1024 * 1024):
                        h.update(chunk)
                        size += len(chunk)
                        if size > entry['bytes']:
                            raise ValueError(f'Archive size exceeds manifest: {entry["path"]}')
                        if temporary is not None:
                            temporary.write(chunk)
                if h.hexdigest() != entry['sha256'] or size != entry['bytes']:
                    raise ValueError(f'Archive integrity mismatch: {entry["path"]}')
                counts['verified'] += 1
                if temporary is not None:
                    temporary.close()
                    os.link(temporary.name, target)  # Atomic, exclusive install: never overwrite a concurrently created file.
                    counts['restored'] += 1
            finally:
                if temporary is not None:
                    temporary.close()
                    Path(temporary.name).unlink(missing_ok=True)
    return counts


def restore_index(index_path,destination,verify_only=False):
    try:
        from .archive_layers import load_index
    except ImportError:
        from archive_layers import load_index
    index_path=Path(index_path);index=load_index(index_path);destination=Path(destination).resolve();groups={a['name']:[] for a in index['archives']}
    for entry in index['files']:
        path=PurePosixPath(entry['path'])
        if path.parts[0] not in ('data','runs'):continue
        target=destination.joinpath(*path.parts)
        if not target.resolve().is_relative_to(destination):raise ValueError('Archive target escapes destination')
        if target.exists() and (not target.is_file() or digest(target)!=entry['sha256']):raise ValueError(f'Existing file differs; refusing overwrite: {path}')
        groups[entry['archive']].append({k:entry[k] for k in ['path','bytes','sha256']})
    # Verify every archive identity before installing any data.
    for archive in index['archives']:
        path=index_path.parent/archive['name']
        if not path.exists():raise ValueError('Archive missing; run git lfs pull: '+archive['name'])
        with path.open('rb') as f:
            if f.read(128).startswith(b'version https://git-lfs.github.com/spec/v1'):raise ValueError('Archive is an LFS pointer; run git lfs pull')
        if path.stat().st_size!=archive['bytes'] or digest(path)!=archive['sha256']:raise ValueError('Archive integrity mismatch: '+archive['name'])
    total={'verified':0,'restored':0,'already_present':0}
    for name,entries in groups.items():
        if not entries:continue
        result=restore(index_path.parent/name,destination,verify_only,entries_override=entries)
        for k in total:total[k]+=result[k]
    return total


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path, nargs='?', default=Path(__file__).resolve().parents[1] / 'datasets/results-complete.zip')
    parser.add_argument('--destination', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    archive=args.archive
    default=Path(__file__).resolve().parents[1]/'datasets/results-complete.zip'
    if archive==default and default.with_name('results-index.json').exists():archive=default.with_name('results-index.json')
    print(json.dumps(restore(archive, args.destination, args.verify_only), indent=2))
