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
    if args.include_computed:
        with urllib.request.urlopen('https://alphafold.ebi.ac.uk/api/prediction/P00748',timeout=30) as response:
            metadata=response.read()
        (raw/'alphafold-api.json').write_bytes(metadata)
        item=json.loads(metadata)[0]
        urls[Path(item['pdbUrl']).name]=item['pdbUrl']
        urls[Path(item['paeDocUrl']).name]=item['paeDocUrl']
        urls['uniprot-features.json']='https://www.ebi.ac.uk/proteins/api/proteins/P00748'

    reports=[]
    for name,url in urls.items():
        path=raw/name
        try:
            if path.exists():content=path.read_bytes();state='reused'
            else:
                with urllib.request.urlopen(urllib.request.Request(url,headers={'Accept':'application/json' if name.endswith('.json') else '*/*'}),timeout=30) as response:content=response.read()
                state='downloaded'
            if name.endswith('.cif') and (not content.lstrip().startswith(b'data_') or b'_atom_site.' not in content):
                raise ValueError('Response is not an atomic-coordinate mmCIF')
            if name.endswith('.pdb') and b'ATOM  ' not in content:raise ValueError('No PDB atom records')
            if name.endswith('.json'):json.loads(content)
            if not path.exists():path.write_bytes(content)
            row={'name':name,'source_url':url,'state':state,'size_bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()}
        except Exception as error:row={'name':name,'source_url':url,'state':'unavailable','error':str(error)}
        reports.append(row);print(json.dumps(row),flush=True)
    (OUT/'downloads.json').write_text(json.dumps(reports,indent=2)+'\n')
if __name__=='__main__':main()
