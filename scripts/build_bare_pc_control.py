"""Restore PC methyl C-H and remove coating using audited PC4 parameters.
No graft fragment charges are reused as bare-PC charges.
"""
import argparse,json,shutil
from pathlib import Path
import numpy as np
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'runs/surface-expanded-18x12/dry.gro');p.add_argument('--output',type=Path,default=ROOT/'runs/bare-pc-control-preparation');p.add_argument('--retained-grafts',type=int,choices=range(7),default=0);a=p.parse_args()
    if a.output.exists():raise RuntimeError('Existing bare-PC preparation preserved')
    a.output.mkdir()
    g=Chem.MolFromMolFile(str(ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol'),removeHs=False);pc=Chem.MolFromMolFile(str(ROOT/'models/pc-graft/bpa-pc-dp4-phenolic-ends.mol'),removeHs=False)
    offset=g.GetNumHeavyAtoms()-pc.GetNumHeavyAtoms();nh=pc.GetNumHeavyAtoms()
    for i in range(nh):
        if pc.GetAtomWithIdx(i).GetAtomicNum()!=g.GetAtomWithIdx(offset+i).GetAtomicNum():raise ValueError('PC atom mapping mismatch')
    for b in pc.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if max(i,j)<nh:
            mapped=g.GetBondBetweenAtoms(offset+i,offset+j)
            if mapped is None or mapped.GetBondType()!=b.GetBondType():raise ValueError('PC graph mismatch')
    lines=a.source.read_text().splitlines();raw=lines[2:2+32766]
    if len(raw)!=32766:raise ValueError('Expected six copies of5461-atom coating cell')
    xyz=np.array([[float(s[20:28]),float(s[28:36]),float(s[36:44])] for s in raw]);all_pc=[];new_hydrogens=0
    for cell in range(6):
        start=cell*5461;part=xyz[start:start+445];restored=np.zeros((132,3));restored[:nh]=part[offset:offset+nh]
        for i in range(nh):
            old_h=sorted(n.GetIdx() for n in pc.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum()==1)
            available=sorted(n.GetIdx() for n in g.GetAtomWithIdx(offset+i).GetNeighbors() if n.GetAtomicNum()==1)
            for h,q in zip(old_h,available):restored[h]=part[q]
            if len(old_h)!=len(available):
                if len(old_h)!=len(available)+1:raise ValueError('Unexpected H chemistry difference')
                external=[n.GetIdx() for n in g.GetAtomWithIdx(offset+i).GetNeighbors() if n.GetIdx()<offset]
                if len(external)!=1:raise ValueError('Missing graft partner')
                vector=part[external[0]]-part[offset+i];restored[old_h[-1]]=part[offset+i]+0.109*vector/np.linalg.norm(vector);new_hydrogens+=1
        lengths=[np.linalg.norm(restored[b.GetBeginAtomIdx()]-restored[b.GetEndAtomIdx()]) for b in pc.GetBonds()]
        if min(lengths)<0.07 or max(lengths)>0.25:raise ValueError('Restored PC bond geometry abnormal')
        all_pc.append((part,'GRAFT') if cell<a.retained_grafts else (restored,'PC4'))
        for free in range(38):all_pc.append((xyz[start+445+free*132:start+445+(free+1)*132],'PC4'))
    if new_hydrogens!=6:raise ValueError('Expected six restored C-H sites')
    templates={name:(ROOT/'models/forcefield'/f'{name}-{name}_GMX.gro').read_text().splitlines()[2:-1] for name in ['PC4','GRAFT']};records=[];atom=0
    for residue,(molecule,name) in enumerate(all_pc,1):
        for row,(x,y,z) in zip(templates[name],molecule):
            atom+=1;records.append(f'{residue:5d}{row[5:15]}{atom%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}')
    (a.output/'dry.gro').write_text('Bare PC restored chemistry, provisional control\n'+str(atom)+'\n'+'\n'.join(records)+'\n'+lines[-1]+'\n')
    for file in ['PC4.itp','GRAFT.itp','gaff2-types.itp']:shutil.copy2(ROOT/'runs/surface-nvt-10ps'/file,a.output/file)
    (a.output/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "gaff2-types.itp"\n#include "PC4.itp"\n#include "GRAFT.itp"\n#include "amber99sb-ildn.ff/tip3p.itp"\n[ system ]\nBare-PC provisional control\n[ molecules ]\n'+''.join('GRAFT 1\nPC4 38\n' if i<a.retained_grafts else 'PC4 39\n' for i in range(6)))
    report={'state':('coverage-variant coordinates/topology prepared; not minimized' if a.retained_grafts else 'bare-PC coordinates/topology prepared; not minimized'),'PC4_molecules':234-a.retained_grafts,'graft_molecules':a.retained_grafts,'dry_atoms':atom,'restored_methyl_hydrogens':6-a.retained_grafts,'source':str(a.source),'chemistry':'whole PC4 topologies/AM1-BCC charges reused; graft polymer removed and missing methyl H restored along original attachment direction','scope':'same PC heavy-atom conformations at preparation, not independently equilibrated bulk PC','limitations':['PC4 chain length and low density remain provisional','restore/minimize again from accepted coated surface before matched adsorption comparisons','no adsorption or equilibrium result']}
    (a.output/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
