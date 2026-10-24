"""Publish plots from a completed real protein pilot, keeping scope visible."""
import argparse,csv,json,os,shutil
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/protein-pilot-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def publish(source,destination):
    source=Path(source);out=Path(destination)
    report=json.loads((source/'result.json').read_text())
    if report['status']!='completed':raise ValueError('Run is not complete')
    with (source/'observables.csv').open() as f:rows=list(csv.DictReader(f))
    if int(rows[-1]['step'])!=report['steps']:raise ValueError('Final frame missing')
    t=[float(r['time_ps']) for r in rows]
    if len(rows)==1:
        import numpy as np
        from openmm.app import PDBFile
        from openmm import unit
        def ca(name):
            p=PDBFile(str(source/name))
            return np.array([p.positions[a.index].value_in_unit(unit.angstrom) for a in p.topology.atoms() if a.name=='CA'])
        ref=ca('minimized.pdb');final=ca('final.pdb');a=final-final.mean(axis=0);b=ref-ref.mean(axis=0)
        u,_,vt=np.linalg.svd(a.T@b);correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
        displacement=np.linalg.norm(a@u@correction@vt-b,axis=1)
        fig,axes=plt.subplots(1,2,figsize=(12,4.8))
        axes[0].bar(['Prepared','After minimization','After 0.1 ps MD'],[report['initial_potential_kJ_mol'],report['minimized_potential_kJ_mol'],report['final_potential_kJ_mol']],color=['#98a7b1','#54a4a4','#7768af'])
        axes[0].set_ylabel('Total potential energy (kJ/mol)')
        axes[1].plot(np.arange(143,143+len(displacement)),displacement,color='#7768af',lw=1)
        axes[1].set_xlabel('Source gamma-chain sequence position');axes[1].set_ylabel('Aligned C-alpha displacement (Å)')
        fig.suptitle(f"Actual uploaded 3GHG fragment: {report['residues']} residues / {report['particles']:,} atoms",fontsize=16)
        fig.text(.5,.075,f"100 steps × 1 fs = 0.1 ps | final C-alpha RMSD: {report['final_CA_RMSD_from_minimized_A']:.3f} Å | Amber14 + GBn2",ha='center',fontsize=10)
        fig.text(.5,.025,'Protein only; no polymer. Minimization convergence and equilibration not established. No adsorption prediction.',ha='center',color='#a12d32',fontsize=10)
        fig.tight_layout(rect=(0,.13,1,.92));out.mkdir(parents=True,exist_ok=True)
        fig.savefig(out/'results.png',dpi=160);fig.savefig(out/'results.svg')
        for name in ('result.json','observables.csv','prepared.pdb','minimized.pdb','final.pdb'):shutil.copyfile(source/name,out/name)
        print('Published single-frame pilot (not a time series)',out)
        return
    fig,axes=plt.subplots(1,3,figsize=(14,4.8))
    for ax,col,label,color in zip(axes,['potential_kJ_mol','temperature_K','CA_RMSD_from_minimized_A'],['Potential energy (kJ/mol)','Temperature (K)','Aligned C-alpha RMSD (Å)'],['#228793','#b97835','#8155b3']):
        ax.plot(t,[float(r[col]) for r in rows],'o-',color=color,ms=3)
        ax.set_xlabel('Time (ps)');ax.set_ylabel(label);ax.grid(alpha=.2)
    axes[1].axhline(310,color='gray',ls='--',label='310 K target');axes[1].legend(fontsize=8)
    fig.suptitle(f"Actual 3GHG gamma-chain fragment: {report['residues']} residues / {report['particles']:,} atoms",fontsize=16)
    fig.text(.5,.07,'Amber14 + GBn2 implicit solvent | CPU | 1 fs timestep | source chain C, sequence 143–394',ha='center',fontsize=10)
    fig.text(.5,.02,'Protein preparation pilot only. No polymer present; no adsorption, affinity, or equilibration conclusion.',ha='center',color='#a12d32',fontsize=10)
    fig.tight_layout(rect=(0,.12,1,.92));out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/'results.png',dpi=160);fig.savefig(out/'results.svg')
    for name in ('result.json','observables.csv','prepared.pdb','final.pdb'):shutil.copyfile(source/name,out/name)
    print('Published',out)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('destination');a=p.parse_args();publish(a.source,a.destination)
