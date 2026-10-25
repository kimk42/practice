"""Recover GRAFT AM1-BCC/topology from a converged SQM output, without rerunning SQM."""
import argparse,hashlib,json,os,shutil,subprocess,tempfile
from pathlib import Path
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('sqm_output',type=Path);args=parser.parse_args()
    source=args.sqm_output.resolve()
    if 'Calculation Completed' not in source.read_text():raise ValueError('SQM completion marker absent')
    out=Path(tempfile.mkdtemp(prefix='graft-completed-sqm-',dir=ROOT/'runs'));work=out/'GRAFT.acpype';work.mkdir()
    shutil.copy2(source,out/'completed-sqm.out')
    evidence={'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'method':'antechamber sqmout input reads converged AM1 charges then applies BCC and standard equivalent-atom averaging; no new quantum optimization'}
    (out/'provenance.json').write_text(json.dumps(evidence,indent=2)+'\n')
    def run(command,name,cwd):
        with (out/(name+'.log')).open('w') as log:subprocess.run(command,cwd=cwd,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['antechamber','-i',str(out/'completed-sqm.out'),'-fi','sqmout','-o','GRAFT_bcc_gaff2.mol2','-fo','mol2','-c','bcc','-at','gaff2','-nc','0','-m','1','-dr','no','-pf','n','-s','2'],'bcc-recovery',work)
    mol2=work/'GRAFT_bcc_gaff2.mol2';sections={};section=None
    for line in mol2.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):section=line[9:];sections[section]=[]
        elif section and line.strip():sections[section].append(line.split())
    original=ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol';m=Chem.MolFromMolFile(str(original),removeHs=False)
    expected={tuple(sorted((b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1))) for b in m.GetBonds()}
    actual={tuple(sorted(map(int,r[1:3]))) for r in sections['BOND']}
    if expected!=actual:raise ValueError('SQM coordinate-derived bond graph differs from original')
    numbers={'c':6,'h':1,'n':7,'o':8}
    if [a.GetAtomicNum() for a in m.GetAtoms()]!=[numbers[a[5][0]] for a in sections['ATOM']]:raise ValueError('Atom order/element mismatch')
    if abs(sum(float(a[8]) for a in sections['ATOM']))>1e-4:raise ValueError('BCC net charge inconsistent with neutral input')
    run(['acpype','-i',str(original),'-b','GRAFT','-a','gaff2','-c','bcc','-n','0','-m','1','-s','86400'],'GRAFT-finalize',out)
    # All audits must pass before this candidate replaces the active directory.
    manifest=ROOT/'runs/workflow-active-parameters.json';old=json.loads(manifest.read_text());updated=dict(old);updated['GRAFT']=str(out.relative_to(ROOT))
    manifest.write_text(json.dumps(updated,indent=2)+'\n')
    try:
        from surface_validation import parameter_gate
        if parameter_gate('GRAFT') is None:raise ValueError('Topology incomplete')
    except Exception:
        manifest.write_text(json.dumps(old,indent=2)+'\n');raise
    print(str(out),flush=True)
if __name__=='__main__':main()
