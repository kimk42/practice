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


def restore(archive, destination, verify_only=False):
    destination = Path(destination).resolve()
    counts = {'verified': 0, 'restored': 0, 'already_present': 0}
    with zipfile.ZipFile(archive) as z:
        entries = json.loads(z.read('FILES.json'))['files']
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--destination', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    print(json.dumps(restore(args.archive, args.destination, args.verify_only), indent=2))
