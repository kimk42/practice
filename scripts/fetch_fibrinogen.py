"""Download the experimental 3GHG reference assembly, not a complete native model."""
import hashlib,json,urllib.request
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]/'research/proteins/fibrinogen'
def main():
    raw=OUT/'raw';raw.mkdir(parents=True,exist_ok=True);reports=[]
    for name,url in {
        '3GHG.cif':'https://files.rcsb.org/download/3GHG.cif',
        '3GHG-assembly1.cif':'https://files.rcsb.org/download/3GHG-assembly1.cif',
        '3GHG-entry.json':'https://data.rcsb.org/rest/v1/core/entry/3GHG',
    }.items():
        path=raw/name
        if path.exists():content=path.read_bytes();state='reused'
        else:
            with urllib.request.urlopen(url,timeout=60) as response:content=response.read()
            state='downloaded'
        if name.endswith('.cif') and (not content.lstrip().startswith(b'data_') or b'_atom_site.' not in content):raise ValueError('Invalid mmCIF')
        if name.endswith('.json'):json.loads(content)
        if not path.exists():path.write_bytes(content)
        reports.append({'file':name,'url':url,'sha256':hashlib.sha256(content).hexdigest(),'bytes':len(content),'state':state})
    (OUT/'downloads.json').write_text(json.dumps(reports,indent=2)+'\n')
if __name__=='__main__':main()
