"""Repeat interrupted FXII minimization without overwriting original outputs."""
import json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'runs/fxii-unglycosylated-topology'
def main():
    if (OUT/'result.json').exists():return
    if (OUT/'em-recovery.log').exists():raise RuntimeError('Recovery output exists; inspect before repeating')
    def run(args,name):
        with (OUT/(name+'-launch.log')).open('w') as log:
            subprocess.run(['gmx',*args],cwd=OUT,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['mdrun','-s','em.tpr','-deffnm','em-recovery','-ntmpi','1','-ntomp','2'],'em-recovery')
    if 'converged to Fmax <' not in (OUT/'em-recovery.log').read_text():raise RuntimeError('EM force threshold not met')
    # Same numerical pilot conditions as validate_fxii_control.py.
    (OUT/'md-recovery.mdp').write_text('integrator = md\ndt = 0.0005\nnsteps = 400\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\nnstenergy = 20\npbc = xyz\n')
    run(['grompp','-f','md-recovery.mdp','-c','em-recovery.gro','-p','topol.top','-o','md-recovery.tpr'],'md-recovery-prep')
    run(['mdrun','-deffnm','md-recovery','-ntmpi','1','-ntomp','2'],'md-recovery')
    run(['check','-e','md-recovery.edr'],'energy-recovery-check')
    report={'state':'unglycosylated predicted FXII numerical pilot passed','duration_ps':0.2,'temperature_K':300,'scope':'predicted unglycosylated protein alone, neutralizing ions only; not native validation or surface adsorption','forcefield':'amber99sb-ildn / TIP3P','restart':'EM repeated from original TPR; interrupted output preserved'}
    (OUT/'result.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
