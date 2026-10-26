"""Retrieve public glycosylation citations and fibrinogen chain annotations."""
import json,urllib.request,concurrent.futures,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    fx=ROOT/'research/proteins/fxii/raw';fg=ROOT/'research/proteins/fibrinogen/raw';features=json.loads((fx/'uniprot-features.json').read_text())['features']
    sites=[f for f in features if f['type']=='CARBOHYD'];pmids=sorted({e['source']['id'] for f in sites for e in f.get('evidences',[]) if e.get('source',{}).get('name')=='PubMed'})
    tasks=[(fx/f'pubmed-{p}.json',f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{p}%20AND%20SRC:MED&resultType=core&format=json') for p in pmids]
    tasks += [(fg/f'{p}-uniprot.json',f'https://www.ebi.ac.uk/proteins/api/proteins/{p}') for p in ['P02671','P02675','P02679']]
    def fetch(item):
        path,url=item
        try:
            if path.exists():data=path.read_bytes()
            else:
                with urllib.request.urlopen(urllib.request.Request(url,headers={'Accept':'application/json'}),timeout=30) as r:data=r.read()
                json.loads(data);path.write_bytes(data)
            return {'file':str(path.relative_to(ROOT)),'url':url,'sha256':hashlib.sha256(data).hexdigest(),'state':'retrieved'}
        except Exception as e:return {'file':str(path.relative_to(ROOT)),'url':url,'state':'unavailable','reason':str(e)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:reports=list(pool.map(fetch,tasks))
    out=ROOT/'research/proteins/modeling-evidence';out.mkdir(exist_ok=True)
    (out/'downloads.json').write_text(json.dumps(reports,indent=2)+'\n')
    (out/'fxii-glycosylation-sites.json').write_text(json.dumps([{'position':f['begin'],'annotation':f['description'],'pubmed_ids':[e['source']['id'] for e in f.get('evidences',[]) if e.get('source',{}).get('name')=='PubMed']} for f in sites],indent=2)+'\n')
    chains=[]
    for p in ['P02671','P02675','P02679']:
        file=fg/f'{p}-uniprot.json'
        if not file.exists():continue
        d=json.loads(file.read_text());chains.append({'accession':p,'canonical_precursor_residues':d['sequence']['length'],'signal_and_chain_annotations':[f for f in d['features'] if f['type'] in ['SIGNAL','CHAIN','VAR_SEQ']],'alternative_products':[c for c in d.get('comments',[]) if c.get('type')=='ALTERNATIVE_PRODUCTS'],'caution':'canonical sequence may not match the intended plasma isoform/experimental construct; do not fill missing structure blindly'})
    (out/'fibrinogen-chain-annotations.json').write_text(json.dumps(chains,indent=2)+'\n')
    isoforms=[]
    for accession in ['P02671-2','P02679-2']:
        path=fg/(accession+'.fasta');url=f'https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id={accession}&format=fasta&style=raw'
        try:
            if path.exists():data=path.read_bytes()
            else:
                with urllib.request.urlopen(url,timeout=30) as response:data=response.read()
            if not data.startswith(b'>') or accession not in data.decode().splitlines()[0]:raise ValueError('Invalid isoform FASTA')
            path.write_bytes(data);isoforms.append({'accession':accession,'url':url,'state':'retrieved','length':len(''.join(data.decode().splitlines()[1:])),'sha256':hashlib.sha256(data).hexdigest()})
        except Exception as error:isoforms.append({'accession':accession,'url':url,'state':'unavailable','reason':str(error)})
    (out/'isoform-downloads.json').write_text(json.dumps(isoforms,indent=2)+'\n')

if __name__=='__main__':main()
