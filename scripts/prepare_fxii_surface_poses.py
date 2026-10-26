"""Prepare three explicitly unglycosylated FXII pilot poses; no adsorption claim."""
import json,subprocess,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'runs/fxii-surface-poses'
    if out.exists():raise RuntimeError('Existing poses preserved')
    out.mkdir();source=ROOT/'runs/fxii-unglycosylated-topology'
    with (out/'whole-launch.log').open('w') as log:
        subprocess.run(['gmx','trjconv','-s',str(source/'md-recovery.tpr'),'-f',str(source/'md-recovery.gro'),'-o','whole.gro','-pbc','mol'],cwd=out,input='0\n',text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    lines=(out/'whole.gro').read_text().splitlines();raw=lines[2:2+8995]
    xyz=np.array([[float(x[20:28]),float(x[28:36]),float(x[36:44])] for x in raw])
    heavy=np.array([not x[10:15].strip().startswith('H') for x in raw]);center=xyz[heavy].mean(axis=0);centered=xyz-center
    eigenvalues,axes=np.linalg.eigh(np.cov(centered[heavy].T));axes=axes[:,np.argsort(eigenvalues)[::-1]]
    if np.linalg.det(axes)<0:axes[:,-1]*=-1
    principal=centered@axes
    rows=[]
    for index,degrees in enumerate([0,120,240]):
        angle=np.deg2rad(degrees);c,s=np.cos(angle),np.sin(angle);rotation=np.array([[1,0,0],[0,c,-s],[0,s,c]])
        pose=principal@rotation.T;extent=np.ptp(pose,axis=0)
        fits=bool(extent[0]+2.4<=18 and extent[1]+2.4<=12)
        name=f'pose-{index+1}.gro'
        records=[f'{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}' for line,(x,y,z) in zip(raw,pose)]
        (out/name).write_text('Unglycosylated FXII rigid pose; not yet placed on surface\n8995\n'+'\n'.join(records)+'\n18.00000 12.00000 18.00000\n')
        rows.append({'file':name,'rotation_about_principal_axis_degrees':degrees,'extents_nm':extent.tolist(),'fits_lateral_padding_1_2nm':fits})
    meta={'source':'predicted unglycosylated mature FXII after0.2ps water numerical pilot','source_sha256':hashlib.sha256((source/'md-recovery.gro').read_bytes()).hexdigest(),'protein_atoms':8995,'poses':rows,'protein_glycans_present':False,'status':'rigid orientation candidates only; placement/clashes/charge/EM still required','limitations':['correlated rigid orientations of one predicted conformation, not independent ensembles','no physiological glycoform, salt or protonation validation','not protein equilibrium or adsorption sampling']}
    (out/'poses.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta),flush=True)
if __name__=='__main__':main()
