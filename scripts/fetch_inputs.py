"""Fetch public source inputs; retain publisher checksums and coordinate hashes."""
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'research'

def fetch(url):
    with urllib.request.urlopen(url, timeout=40) as response:
        return response.read()

def main():
    OUT.mkdir(exist_ok=True)
    inventory = {}
    for article_id in [12906763, 7178012]:
        metadata = json.loads(fetch(f'https://api.figshare.com/v2/articles/{article_id}'))
        for item in metadata['files']:
            path = OUT / Path(item['name']).name
            content = path.read_bytes() if path.exists() else fetch(item['download_url'])
            if len(content) != item['size'] or hashlib.md5(content).hexdigest() != item['computed_md5']:
                raise ValueError(f'Publisher checksum/size mismatch: {path.name}; existing file preserved')
            if not path.exists():
                path.write_bytes(content)
            inventory[path.name] = {'source': item['download_url'], 'sha256': hashlib.sha256(content).hexdigest(), 'publisher_md5': item['computed_md5']}
    for pdb_id in ['3GHG', '8OS5']:
        path = OUT / f'{pdb_id}.cif'
        url = f'https://files.rcsb.org/download/{pdb_id}.cif'
        content = path.read_bytes() if path.exists() else fetch(url)
        if not content.lstrip().startswith(b'data_') or f'_entry.id'.encode() not in content:
            raise ValueError(f'Invalid mmCIF: {pdb_id}')
        if not path.exists():
            path.write_bytes(content)
        inventory[path.name] = {'source': url, 'sha256': hashlib.sha256(content).hexdigest()}
    (OUT / 'input-inventory.json').write_text(json.dumps(inventory, indent=2) + '\n')
    print('Source inputs verified:', ', '.join(inventory))

if __name__ == '__main__':
    main()
