"""Build a correlated 3x2 tiled surface candidate from the audited 10ps prototype.
No new experimental density/coverage or protein glycoform is asserted.
"""
import json,hashlib,subprocess,shutil
from pathlib import Path
import numpy as np
from rdkit import Chem
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'runs/surface-expanded-18x12'
    if out.exists():raise RuntimeError('Existing expanded candidate preserved')
    out.mkdir()
    source=ROOT/'runs/surface-nvt-10ps';result=json.loads((source/'result.json').read_text())
    if result['state']!='10ps surface numerical follow-up completed':raise ValueError('Source pilot incomplete')
    def run(args,name,stdin=None):
        with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['trjconv','-f',str(source/'nvt.gro'),'-s',str(source/'nvt.tpr'),'-o','whole.gro','-pbc','mol'],'whole','0\n')
    lines=(out/'whole.gro').read_text().splitlines();raw=lines[2:2+5461]
    xyz=np.array([[float(s[20:28]),float(s[28:36]),float(s[36:44])] for s in raw])
    molecules=[(0,445,'GRAFT')]+[(445+i*132,445+(i+1)*132,'PC4') for i in range(38)]
    models={'GRAFT':Chem.MolFromMolFile(str(ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.mol'),removeHs=False),'PC4':Chem.MolFromMolFile(str(ROOT/'models/pc-graft/bpa-pc-dp4-phenolic-ends.mol'),removeHs=False)}
    maximum=0.
    for first,last,name in molecules:
        part=xyz[first:last]
        lengths=[np.linalg.norm(part[b.GetBeginAtomIdx()]-part[b.GetEndAtomIdx()]) for b in models[name].GetBonds()]
        maximum=max(maximum,max(lengths))
        if min(lengths)<0.07 or max(lengths)>0.25:raise ValueError('Source molecule is not whole or has abnormal bond length')
        xyz[first:last,:2]-=np.floor(part[:,:2].mean(axis=0)/6)*6
    xyz[:,2]+=1.-xyz[:,2].min()
    records=[];number=0;residue=0
    for ix in range(3):
        for iy in range(2):
            for first,last,name in molecules:
                residue+=1
                for i in range(first,last):
                    number+=1;x,y,z=xyz[i]+[ix*6.,iy*6.,0.]
                    records.append(f'{residue%100000:5d}{raw[i][5:10]}{raw[i][10:15]}{number%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}')
    (out/'dry.gro').write_text('Correlated tiled provisional coating; no proteins\n'+str(len(records))+'\n'+'\n'.join(records)+'\n  18.00000 12.00000 18.00000\n')
    for file in ['GRAFT.itp','PC4.itp','gaff2-types.itp']:shutil.copy2(source/file,out/file)
    topology='#include "amber99sb-ildn.ff/forcefield.itp"\n#include "gaff2-types.itp"\n#include "PC4.itp"\n#include "GRAFT.itp"\n#include "amber99sb-ildn.ff/tip3p.itp"\n[ system ]\nProvisional correlated tiled coating18x12\n[ molecules ]\n'
    text=(source/'topol.top').read_text().split('[ molecules ]')[1].splitlines()
    rows=[s.split() for s in text if s.strip() and not s.lstrip().startswith(';')]
    graft=rows[0][0];pc=rows[1][0]
    topology+=''.join(f'{graft} 1\n{pc} 38\n' for _ in range(6))
    (out/'topol.top').write_text(topology)
    run(['solvate','-cp','dry.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
    atoms=int((out/'solvated.gro').read_text().splitlines()[1]);waters=(atoms-len(records))//3
    report={'state':'expanded solvated starting candidate','source_sha256':hashlib.sha256((source/'nvt.gro').read_bytes()).hexdigest(),'box_nm':[18,12,18],'tile_count':6,'graft_molecules':6,'free_PC4_molecules':228,'dry_atoms':len(records),'water_molecules':waters,'total_atoms':atoms,'source_max_direct_bond_length_nm':maximum,'protein_present':False,'limitations':['six correlated translations of one short pilot, not independent morphology/replicas','same provisional low-density PC and graft coverage, not measured coating','no experimental DP, crosslink distribution, glycans or physiological salt','new solvation requires EM before protein placement or dynamics']}
    (out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
    (out/'em.mdp').write_text('integrator = steep\nnsteps = 5000\nemtol = 1000\nemstep = 0.005\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
    run(['grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],'em-prep')
    run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','4'],'em')
    if 'converged to Fmax <' not in (out/'em.log').read_text():raise RuntimeError('Expanded EM threshold not met')
    report['state']='expanded solvated candidate minimized';report['EM']='Fmax below1000kJ/mol/nm';(out/'result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
