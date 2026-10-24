"""Resolve recorded dataset paths after relocating the project; no file mutation."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def tools_root():
    return Path(os.environ.get('FIBRINOGEN_TOOLS_DIR', ROOT / '.tools')).expanduser().resolve()


def amber_home():
    return Path(os.environ.get('AMBERHOME', tools_root() / 'ambertools')).expanduser().resolve()


def resolve_project_path(recorded, root=ROOT):
    """Prefer this checkout's dataset over an old producer's absolute path."""
    root = Path(root).resolve()
    path = Path(recorded).expanduser()
    if not path.is_absolute():
        candidate = (root / path).resolve()
        if candidate.is_relative_to(root) and candidate.exists():
            return candidate
    else:
        for index, part in enumerate(path.parts):
            if part in ('data', 'runs', 'artifacts'):
                candidate = (root / Path(*path.parts[index:])).resolve()
                if candidate.is_relative_to(root) and candidate.exists():
                    return candidate
        if path.exists():
            return path
    raise FileNotFoundError(f'Dataset not found: {recorded}. Restore the raw-results archive first.')
