"""Build/minimize matched bare PC only after the coated surface EM gate."""
import json,subprocess
from pathlib import Path
from em_recovery import last_em_frame
ROOT=Path(__file__).resolve().parents[1]
def main():
    source=ROOT/'runs/surface-expanded-18x12';out=ROOT/'runs/bare-pc-matched'
    if (out/'result.json').exists():return
    if not (source/'result.json').exists():raise RuntimeError('Coated EM not completed')
    if not out.exists():subprocess.run(['python','scripts/build_bare_pc_control.py','--source',str(source/'em.gro'),'--output',str(out)],cwd=ROOT,check=True)
    def run(args,name):
        with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    if not (out/'em.tpr').exists():
        run(['solvate','-cp','dry.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
        (out/'em.mdp').write_text('integrator = steep\nnsteps = 20000\nemtol = 1000\nemstep = 0.005\nnstxout = 100\nnstenergy = 100\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
        run(['grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],'prep')
    restart=last_em_frame(out)
    if restart:run(['grompp','-f','em.mdp','-c',str(restart),'-p','topol.top','-o','em.tpr'],'recovery-prep')
    run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','4'],'em')
    if 'converged to Fmax <' not in (out/'em.log').read_text():raise RuntimeError('Bare-PC EM criterion not met')
    n=int((out/'em.gro').read_text().splitlines()[1]);report={'state':'surface candidate minimized','EM':'Fmax below1000kJ/mol/nm','dry_atoms':30888,'PC4_molecules':234,'graft_molecules':0,'total_atoms':n,'water_molecules':(n-30888)//3,'box_nm':[18,12,18],'scope':'matched initial PC conformations, restored methyl H, original whole-PC4 charges; not bulk equilibrium or experimental validation'}
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
