"""Pack a finite BPA-PC oligomer patch around one preconstructed graft.

This creates a density-targeted starting geometry, not an equilibrated bulk film.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--density',type=float,default=0.6);parser.add_argument('--resume',type=Path);parser.add_argument('--timeout',type=int,default=1800);args=parser.parse_args()
    if not 0.1 <= args.density <= 1.5: parser.error('Choose an explicitly provisional packing density between0.1 and1.5g/cm3')
    out=ROOT/'models/pc-patch';out.mkdir(parents=True,exist_ok=True)
    run=ROOT/'runs/pc-packing';run.mkdir(parents=True,exist_ok=True)
    pc=Chem.SDMolSupplier(str(ROOT/'models/pc-graft/bpa-pc-dp4-phenolic-ends.sdf'),removeHs=False)[0]
    graft=Chem.SDMolSupplier(str(ROOT/'models/pc-graft/cop10-graft-bpa-pc4-one-site.sdf'),removeHs=False)[0]
    cop_heavy=Chem.RemoveHs(Chem.SDMolSupplier(str(ROOT/'models/copolymer/cbaa1-bma-bpaa-dp10-seed2026.sdf'),removeHs=False)[0]).GetNumAtoms()
    pc_heavy=Chem.RemoveHs(pc).GetNumAtoms()
    pc_ids=[];cop_ids=[]
    for atom in graft.GetAtoms():
        owner=atom.GetIdx() if atom.GetAtomicNum()!=1 else atom.GetNeighbors()[0].GetIdx()
        (pc_ids if cop_heavy<=owner<cop_heavy+pc_heavy else cop_ids).append(atom.GetIdx())
    xyz=np.asarray(graft.GetConformer().GetPositions())
    z=xyz[cop_ids].mean(0)-xyz[pc_ids].mean(0);z/=np.linalg.norm(z)
    helper=np.array([0.,1.,0.]) if abs(z[1])<0.9 else np.array([1.,0.,0.])
    x=np.cross(helper,z);x/=np.linalg.norm(x);y=np.cross(z,x)
    basis=np.column_stack((x,y,z));assert np.linalg.det(basis)>0.99
    xyz=(xyz-xyz[pc_ids].mean(0))@basis
    thickness=max(30.,float(np.ptp(xyz[pc_ids,2]))+4.)
    width=max(60.,2.*float(np.abs(xyz[:,:2]).max())+6.)
    xyz+=np.array([width/2,width/2,thickness/2])
    for i,pos in enumerate(xyz):graft.GetConformer().SetAtomPosition(i,pos)
    Chem.MolToPDBFile(graft,str(run/'graft-oriented.pdb'))
    Chem.MolToPDBFile(pc,str(run/'pc4.pdb'))
    density=args.density # Low-density initial packing, not bulk or measured density.
    pc_mass=Descriptors.MolWt(pc)
    pc_total=round(density*width*width*thickness*1e-24*6.02214076e23/pc_mass)
    count=pc_total-1
    restart=''
    if args.resume:
        restart_path=args.resume.resolve()
        if len(restart_path.read_text().splitlines())!=count:raise ValueError('Restart pose count mismatch')
        restart=f'restart_from {restart_path}\n'
    text=f'''tolerance 2.0
filetype pdb
output patch.pdb
seed 2026
nloop 200
maxit 200
restart_to packing.restart
{restart}
structure graft-oriented.pdb
 number 1
 fixed 0. 0. 0. 0. 0. 0.
end structure
structure pc4.pdb
 number {count}
 inside box 0. 0. 0. {width:.6f} {width:.6f} {thickness:.6f}
end structure
'''
    (run/'packmol.inp').write_text(text)
    binary=shutil.which('packmol') or str(Path(os.environ.get('MD_TOOLS_DIR','/workspace/md-tools'))/'root/usr/bin/packmol')
    with (run/'packmol.inp').open() as inp,(run/'packmol.log').open('w') as log:
        subprocess.run([binary],cwd=run,stdin=inp,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
    content=(run/'packmol.log').read_text()
    if 'Success!' not in content:raise RuntimeError('Packmol did not report success')
    lines=[l for l in (run/'patch.pdb').read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))]
    expected=graft.GetNumAtoms()+count*pc.GetNumAtoms();assert len(lines)==expected
    output_xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lines])
    assert np.max(np.abs(output_xyz[:graft.GetNumAtoms()]-xyz))<0.003
    shutil.copy2(run/'patch.pdb',out/'pc4-grafted-cop10-patch.pdb')
    report={'status':'packed starting geometry only; not equilibrated or force-field validated',
            'substrate':'BPA-PC assumed; actual grade unknown', 'pc_oligomer_dp':4,'free_pc_oligomers':count,
            'pc_units_in_graft':4,'copolymer_repeats':{'CBAA1':7,'BMA':2,'BPAA':1},
            'packing_target_pc_density_g_cm3':density,'actual_mass_over_assigned_pc_region_g_cm3':pc_total*pc_mass/6.02214076e23/(width*width*thickness*1e-24),
            'patch_region_A':[width,width,thickness],'atom_count':expected,'geometry_tolerance_A':2.,'graft_position_validation':'fixed graft preserved within0.003A; no double translation',
            'molecule_order':['one covalently linked COP10-PC4 graft',f'{count} free PC4 oligomers'],
            'limitations':['finite OH-ended short oligomers, not bulk molecular weight','density is a packing target, not equilibrated density','single fully reacted BPAA in DP10 is a trial graft, not measured conversion','no interchain crosslinks or protein','no water or production dynamics yet']}
    (out/'patch.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
