"""Checkpointed, bounded NVT follow-up; completion is not equilibrium proof."""
import argparse,json,subprocess,shutil,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--basename',required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--duration-ps',type=float,default=100);a=p.parse_args()
    source=a.source_dir.resolve();out=a.output_dir.resolve()
    if (out/'result.json').exists():print((out/'result.json').read_text());return
    out.mkdir(exist_ok=True);source_gro=source/(a.basename+'.gro');source_cpt=source/(a.basename+'.cpt')
    spec={'source_gro':str(source_gro),'source_sha256':hashlib.sha256(source_gro.read_bytes()).hexdigest(),'duration_ps':a.duration_ps,'ensemble':'NVT','temperature_K':300,'dt_ps':0.0005}
    if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text())!=spec:raise ValueError('Follow-up specification changed; refusing overwrite')
    (out/'spec.json').write_text(json.dumps(spec,indent=2)+'\n')
    def run(args,name,stdin=None):
        with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    if not (out/'nvt.tpr').exists():
        for file in ['topol.top',*[x.name for x in source.glob('*.itp')]]:shutil.copy2(source/file,out/file)
        shutil.copy2(source_gro,out/'source.gro');shutil.copy2(source_cpt,out/'source.cpt')
        (out/'nvt.mdp').write_text(f'integrator = md\ndt = 0.0005\nnsteps = {round(a.duration_ps/0.0005)}\ncontinuation = yes\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = no\nnstenergy = 1000\nnstlog = 10000\nnstxout-compressed = 10000\npbc = xyz\n')
        run(['grompp','-f','nvt.mdp','-c','source.gro','-t','source.cpt','-p','topol.top','-o','nvt.tpr'],'prep')
    args=['mdrun','-deffnm','nvt','-cpt','1','-maxh','23.5','-ntmpi','1','-ntomp','4']
    if (out/'nvt.cpt').exists():args+=['-cpi','nvt.cpt','-append']
    run(args,'nvt')
    run(['check','-e','nvt.edr'],'energy-check');run(['energy','-f','nvt.edr','-o','energies.xvg'],'energy','Potential\nTemperature\nPressure\n0\n')
    v=np.array([[float(x) for x in line.split()] for line in (out/'energies.xvg').read_text().splitlines() if line and not line.startswith(('#','@'))])
    if not np.isfinite(v).all() or v[-1,0]-v[0,0]<a.duration_ps-0.01:raise RuntimeError('Incomplete/nonfinite follow-up; inspect checkpoint before resuming')
    if 'LINCS WARNING' in (out/'nvt.log').read_text():raise RuntimeError('Constraint warnings')
    result={'state':'bounded numerical follow-up completed','duration_ps':a.duration_ps,'temperature_mean_K':float(v[:,2].mean()),'temperature_second_half_mean_K':float(v[v[:,0]>=np.median(v[:,0]),2].mean()),'potential_first_half_mean_kJ_mol':float(v[v[:,0]<np.median(v[:,0]),1].mean()),'potential_second_half_mean_kJ_mol':float(v[v[:,0]>=np.median(v[:,0]),1].mean()),'scope':'preliminary relaxation monitoring only; no automatic equilibrium, sampling convergence or adsorption conclusion'}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
