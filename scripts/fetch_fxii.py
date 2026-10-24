"""Retrieve FXII references directly from RCSB; distinguish computed models."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/proteins/fxii'
URLS={
 '8OS5.cif':'https://files.rcsb.org/download/8OS5.cif',
 '7PRJ.cif':'https://files.rcsb.org/download/7PRJ.cif',
 '6X0T.cif':'https://files.rcsb.org/download/6X0T.cif',
}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--include-computed',action='store_true');args=parser.parse_args()
    raw=OUT/'raw';raw.mkdir(parents=True,exist_ok=True)
    urls=dict(URLS)
    if args.include_computed:urls['AF_AFP00748F1-modelserver.cif']='https://models.rcsb.org/v1/AF_AFP00748F1/full?encoding=cif'
    reports=[]
    for name,url in urls.items():
        path=raw/name
        try:
            if path.exists():content=path.read_bytes();state='reused'
            else:
                with urllib.request.urlopen(url,timeout=30) as response:content=response.read()
                state='downloaded'
            if not content.lstrip().startswith(b'data_') or b'_atom_site.' not in content:
                raise ValueError('Response is not an atomic-coordinate mmCIF')
            if not path.exists():path.write_bytes(content)
            row={'name':name,'source_url':url,'state':state,'size_bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()}
        except Exception as error:row={'name':name,'source_url':url,'state':'unavailable','error':str(error)}
        reports.append(row);print(json.dumps(row),flush=True)
    (OUT/'downloads.json').write_text(json.dumps(reports,indent=2)+'\n')
if __name__=='__main__':main()
