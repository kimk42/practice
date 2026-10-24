"""Publish recorded PE/protein pilot observables and an actual-coordinate snapshot."""
import argparse,csv,json,os,shutil
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/pe-pilot-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from openmm.app import AmberPrmtopFile,PDBFile,element
from openmm import unit


def publish(bundle,run,output):
    bundle=Path(bundle);run=Path(run);out=Path(output)
    result=json.loads((run/'result.json').read_text())
    with (run/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    if result['status']!='completed' or int(rows[-1]['step'])!=result['steps']:raise ValueError('Incomplete selected run')
    prm=AmberPrmtopFile(str(bundle/'complex.prmtop'));final=PDBFile(str(run/'final.pdb'))
    xyz=np.asarray(final.positions.value_in_unit(unit.nanometer),dtype=float)
    ca=np.array([a.index for a in prm.topology.atoms() if a.residue.index<252 and a.name=='CA'])
    pe=np.array([a.index for a in prm.topology.atoms() if a.residue.name=='PEO' and a.element==element.carbon])
    fig=plt.figure(figsize=(13,7));gs=fig.add_gridspec(2,2,width_ratios=[1.15,1])
    ax=fig.add_subplot(gs[:,0],projection='3d')
    ax.scatter(*xyz[pe].T,s=3,color='#65a877',alpha=.55,label='PE carbon atoms')
    ax.plot(*xyz[ca].T,lw=.7,color='#8054a8',label='Fibrinogen gamma fragment (C-alpha)')
    bounds=np.ptp(xyz[np.concatenate([ca,pe])],axis=0);ax.set_box_aspect(np.maximum(bounds,.1));ax.view_init(elev=18,azim=-60)
    ax.set_xlabel('X (nm)');ax.set_ylabel('Y (nm)');ax.set_zlabel('Z (nm)');ax.set_title('Final actual coordinates (water/ions hidden)');ax.legend(loc='upper left',fontsize=8)
    t=[float(r['time_ps']) for r in rows]
    a=fig.add_subplot(gs[0,1]);a.plot(t,[float(r['min_protein_PE_heavy_distance_nm']) for r in rows],'o-',color='#218e98',ms=3);a.set_ylabel('Closest heavy-atom distance (nm)');a.grid(alpha=.2)
    b=fig.add_subplot(gs[1,1]);b.plot(t,[int(r['contact_residues_within_0.45nm']) for r in rows],'o-',color='#8054a8',ms=3);b.set_xlabel('Time (ps)');b.set_ylabel('Protein residues within 0.45 nm');b.grid(alpha=.2)
    fig.suptitle('Actual mixed-system pilot: fibrinogen fragment + finite PE patch',fontsize=16)
    fig.text(.5,.085,f"{result['atoms']:,} atoms | {result['model_manifest']['polymer_chains']} C24H50 oligomers | ff14SB + GAFF2/AM1-BCC + TIP3P",ha='center',fontsize=10)
    fig.text(.5,.045,'Finite oligomer surface; lower layer restrained. Short pre-positioned run, not an equilibrated adsorption experiment.',ha='center',fontsize=10,color='#a12d32')
    fig.text(.5,.012,'First contact: Thr393 next to artificial fragment terminus. No affinity conclusion; PTFE/PVDF not simulated.',ha='center',fontsize=10,color='#a12d32')
    fig.tight_layout(rect=(0,.12,1,.93));out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/'results.png',dpi=160);fig.savefig(out/'results.svg')
    for name in ('result.json','observables.csv','final.pdb'):shutil.copyfile(run/name,out/name)
    (out/'source-paths.json').write_text(json.dumps({'bundle':str(bundle),'run':str(run)},indent=2))
    print('Published',out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bundle');p.add_argument('run');p.add_argument('output');a=p.parse_args();publish(a.bundle,a.run,a.output)
