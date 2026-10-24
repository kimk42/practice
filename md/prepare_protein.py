"""Fetch and audit a candidate fibrinogen structure without inventing missing regions."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import urllib.request
from datetime import datetime, timezone


def fetch(pdb_id, output):
    pdb_id=pdb_id.upper()
    if not re.fullmatch(r'[0-9][A-Z0-9]{3}',pdb_id):
        raise ValueError('Expected a four-character PDB identifier')
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    urls={'coordinates':f'https://files.rcsb.org/download/{pdb_id}.pdb',
          'metadata':f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'}
    report={'pdb_id':pdb_id,'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'urls':urls,'status':'fetching','candidate_only':True}
    try:
        for key,url in urls.items():
            with urllib.request.urlopen(url,timeout=30) as response: payload=response.read()
            if key=='coordinates':
                text=payload.decode()
                if not any(line.startswith('ATOM  ') for line in text.splitlines()):
                    raise ValueError('Downloaded file has no ATOM records')
                path=out/'source.pdb'
                path.write_bytes(payload)
                report['coordinate_sha256']=hashlib.sha256(payload).hexdigest()
            else:
                metadata=json.loads(payload)
                (out/'entry.json').write_text(json.dumps(metadata,indent=2))
                report['title']=metadata.get('struct',{}).get('title')
        report['status']='downloaded_requires_identity_and_assembly_review'
    except Exception as e:
        report.update(status='blocked',error=f'{type(e).__name__}: {e}')
    (out/'download.json').write_text(json.dumps(report,indent=2))
    return report


def audit(path, output):
    from pdbfixer import PDBFixer
    p=Path(path)
    fixer=PDBFixer(filename=str(p))
    fixer.findMissingResidues();fixer.findNonstandardResidues();fixer.findMissingAtoms()
    residues=list(fixer.topology.residues())
    def residue(r):return {'chain':r.chain.id,'id':r.id,'name':r.name}
    report={
        'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
        'atoms':fixer.topology.getNumAtoms(),
        'chains':[{'id':c.id,'residues':sum(1 for _ in c.residues())} for c in fixer.topology.chains()],
        'missing_segments':[{'chain_index':key[0],'insertion_index':key[1],'residues':names} for key,names in fixer.missingResidues.items()],
        'nonstandard_residues':[{'residue':residue(r),'proposed_replacement':name} for r,name in fixer.nonstandardResidues],
        'missing_heavy_atoms':[{'residue':residue(r),'atoms':[a.name for a in atoms]} for r,atoms in fixer.missingAtoms.items()],
        'missing_terminal_atoms':[{'residue':residue(r),'atoms':names} for r,names in fixer.missingTerminals.items()],
        'heterogen_names':sorted({r.name for r in residues if r.name not in fixer.templates}),
        'status':'audit_only_not_MD_ready',
        'required_review':['Confirm species and fibrinogen identity','Select biological assembly rather than assuming asymmetric unit is complete','Review unresolved regions, glycans, disulfides and protonation','Choose mutually compatible protein/polymer/water force fields'],
        'modifications_performed':[],
    }
    # Header records retain crystallographic missing-residue and assembly evidence.
    text=p.read_text()
    report['header_evidence']=[line for line in text.splitlines() if line.startswith(('TITLE ','COMPND','SOURCE','REMARK 350','REMARK 465'))]
    out=Path(output)
    with out.open('x') as f:json.dump(report,f,indent=2)
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    f=sub.add_parser('fetch');f.add_argument('pdb_id');f.add_argument('output')
    a=sub.add_parser('audit');a.add_argument('pdb');a.add_argument('output')
    args=parser.parse_args()
    report=fetch(args.pdb_id,args.output) if args.command=='fetch' else audit(args.pdb,args.output)
    print(json.dumps(report,indent=2))
    if report['status']=='blocked':raise SystemExit(2)
