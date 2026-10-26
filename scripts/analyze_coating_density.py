"""Descriptive z-density profile; no automatic equilibrium/brush-height claim."""
import argparse,json,hashlib,subprocess,csv
from pathlib import Path
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-dir',type=Path,default=ROOT/'runs/surface-nvt-10ps');p.add_argument('--prefix',default='initial-coating');a=p.parse_args();source=a.source_dir.resolve();out=ROOT/'models/analysis';out.mkdir(exist_ok=True)
    m=Chem.MolFromMolFile(str(ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol'),removeHs=False);reaction=json.loads((ROOT/'models/pc-graft/assumptions.json').read_text())['mapped_reaction_atoms_1based'];i,j=reaction['BPAA_C']-1,reaction['PC_C']-1
    if m.GetBondBetweenAtoms(i,j) is None:raise ValueError('Expected graft bond absent')
    seen={i};todo=[i]
    while todo:
        k=todo.pop()
        for neighbor in m.GetAtomWithIdx(k).GetNeighbors():
            n=neighbor.GetIdx()
            if {k,n}=={i,j}:continue
            if n not in seen:seen.add(n);todo.append(n)
    total=int((source/'nvt.gro').read_text().splitlines()[1]);groups={'Coating_polymer':[x+1 for x in sorted(seen)],'PC':[x+1 for x in range(445) if x not in seen]+list(range(446,5462)),'Water':list(range(5462,total+1))}
    if total!=27826:raise ValueError('This analysis expects the original small coating/water composition')
    index=''
    for name,atoms in groups.items():index+='[ '+name+' ]\n'+'\n'.join(' '.join(map(str,atoms[k:k+15])) for k in range(0,len(atoms),15))+'\n'
    (source/'analysis.ndx').write_text(index)
    with (source/'density-launch.log').open('w') as log:subprocess.run(['gmx','density','-s','nvt.tpr','-f','nvt.xtc','-n','analysis.ndx','-d','Z','-sl','80','-dens','mass','-ng','3','-o','density-profile.xvg'],cwd=source,input='0\n1\n2\n',text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    values=[[float(x) for x in line.split()] for line in (source/'density-profile.xvg').read_text().splitlines() if line and not line.startswith(('#','@'))]
    with (out/(a.prefix+'-density.csv')).open('w') as f:w=csv.writer(f,lineterminator="\n");w.writerow(['z_nm','coating_polymer_kg_m3','PC_kg_m3','water_kg_m3']);w.writerows(values)
    evidence={'source':str(source.relative_to(ROOT)),'trajectory_sha256':hashlib.sha256((source/'nvt.xtc').read_bytes()).hexdigest(),'topology_sha256':hashlib.sha256((source/'nvt.tpr').read_bytes()).hexdigest(),'groups_atoms':{k:len(v) for k,v in groups.items()},'bins':80,'scope':'descriptive profile only; not equilibrated bulk density or validated brush height'}
    (out/(a.prefix+'-density-provenance.json')).write_text(json.dumps(evidence,indent=2)+'\n')
if __name__=='__main__':main()
