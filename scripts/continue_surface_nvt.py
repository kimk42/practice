"""Run a 10ps fixed-volume numerical surface follow-up, not production MD."""
import json,subprocess,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    source=ROOT/'runs/surface-validation';status=json.loads((ROOT/'runs/surface-validation-status.json').read_text())
    if status['state']!='passed_initialization_check':raise RuntimeError('Initial surface checks have not passed')
    out=ROOT/'runs/surface-nvt-10ps'
    if out.exists():raise RuntimeError('Existing run preserved; inspect rather than overwrite')
    out.mkdir()
    for name in ['topol.top','gaff2-types.itp','PC4.itp','GRAFT.itp','check.gro','check.cpt']:shutil.copy2(source/name,out/name)
    (out/'nvt.mdp').write_text('integrator = md\ndt = 0.0005\nnsteps = 20000\ncontinuation = yes\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = no\nnstenergy = 1000\nnstlog = 1000\nnstxout-compressed = 2000\npbc = xyz\n')
    def run(args,name,stdin=None):
        with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['grompp','-f','nvt.mdp','-c','check.gro','-t','check.cpt','-p','topol.top','-o','nvt.tpr'],'prepare')
    run(['mdrun','-deffnm','nvt','-ntmpi','1','-ntomp','4'],'nvt')
    run(['check','-e','nvt.edr'],'energy-check')
    run(['energy','-f','nvt.edr','-o','energies.xvg'],'energy','Potential\nTemperature\nPressure\n0\n')
    values=np.array([[float(x) for x in line.split()] for line in (out/'energies.xvg').read_text().splitlines() if line and not line.startswith(('#','@'))])
    if not np.isfinite(values).all():raise RuntimeError('Nonfinite recorded energy/temperature')
    if values[-1,0]-values[0,0]<9.99:raise RuntimeError('Trajectory duration incomplete')
    text=(out/'nvt.log').read_text()
    if 'LINCS WARNING' in text or 'constraint error' in text.lower():raise RuntimeError('Constraint issue')
    report={'state':'10ps surface numerical follow-up completed','duration_ps':10,'atoms':int((out/'nvt.gro').read_text().splitlines()[1]),'energy_frames':len(values),'temperature_mean_K':float(values[:,2].mean()),'protein_present':False,'scope':'fixed-volume low-density short-PC single-graft prototype; short numerical check, not equilibrated experimental coating or adsorption production'}
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
