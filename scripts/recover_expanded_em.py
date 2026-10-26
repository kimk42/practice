"""Recover expanded EM with periodic coordinates; never infer convergence from silence."""
import json,subprocess,re,shutil,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'runs/surface-expanded-18x12'
def main():
    if (BASE/'result.json').exists():return
    out=Path(tempfile.mkdtemp(prefix='expanded-em-recovery-',dir=ROOT/'runs'))
    coordinate=BASE/'solvated.gro';candidates=sorted((ROOT/'runs').glob('expanded-em-recovery-*/em.trr'),key=lambda p:p.stat().st_mtime)
    for trajectory in reversed(candidates):
        check=subprocess.run(['gmx','check','-f',str(trajectory)],capture_output=True,text=True)
        matches=re.findall(r'Last frame\s+\d+\s+time\s+([\d.eE+-]+)',check.stderr+check.stdout)
        if not matches:continue
        with (out/'extract.log').open('w') as log:
            subprocess.run(['gmx','trjconv','-f',str(trajectory),'-s',str(trajectory.with_suffix('.tpr')),'-dump',matches[-1],'-o',str(out/'restart.gro')],input='0\n',text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
        coordinate=out/'restart.gro';break
    for name in ['topol.top','GRAFT.itp','PC4.itp','gaff2-types.itp']:shutil.copy2(BASE/name,out/name)
    (out/'em.mdp').write_text('integrator = steep\nnsteps = 20000\nemtol = 1000\nemstep = 0.005\nnstxout = 100\nnstenergy = 100\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
    def run(args,name):
        with (out/(name+'.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    (BASE/'recovery.json').write_text(json.dumps({'directory':str(out.relative_to(ROOT)),'input':str(coordinate),'periodic_coordinate_steps':100,'state':'running; convergence not yet established'},indent=2)+'\n')
    run(['grompp','-f','em.mdp','-c',str(coordinate),'-p','topol.top','-o','em.tpr'],'prep')
    run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','4'],'run')
    if 'converged to Fmax <' not in (out/'em.log').read_text():raise RuntimeError('Expanded EM force criterion not met')
    shutil.copy2(out/'em.gro',BASE/'em.gro');shutil.copy2(out/'em.tpr',BASE/'em.tpr')
    report=json.loads((BASE/'preparation.json').read_text());report.update({'state':'expanded solvated candidate minimized','EM':'Fmax below1000kJ/mol/nm','recovery_directory':str(out.relative_to(ROOT))})
    (BASE/'result.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
