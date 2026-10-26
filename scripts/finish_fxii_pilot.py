"""Finish/recover strict EM and short MD phases of a prepared FXII pilot."""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
from em_recovery import last_em_frame
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);a=p.parse_args();out=a.directory.resolve()
    if (out/'result.json').exists():return
    spec=json.loads((out/'pilot-spec.json').read_text())
    def run(args,name,stdin=None):
        with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    em_pass=(out/'em.gro').exists() and 'converged to Fmax <' in (out/'em.log').read_text()
    if not em_pass:
        restart=last_em_frame(out)
        if restart:run(['grompp','-f','em.mdp','-c',str(restart),'-p','topol.top','-o','em.tpr'],'em-recovery-prep')
        run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','4'],'em')
        if 'converged to Fmax <' not in (out/'em.log').read_text():raise RuntimeError('Combined EM criterion not met')
    if not (out/'md.tpr').exists():
        (out/'md.mdp').write_text(f'integrator = md\ndt = 0.0005\nnsteps = 400\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = {spec["seed"]}\nnstenergy = 20\nnstxout-compressed = 100\npbc = xyz\n')
        run(['grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'md-prep')
    if not (out/'md.gro').exists():
        args=['mdrun','-deffnm','md','-cpt','1','-ntmpi','1','-ntomp','4']
        if (out/'md.cpt').exists():args+=['-cpi','md.cpt','-append']
        run(args,'md')
    run(['check','-e','md.edr'],'energy-check');run(['energy','-f','md.edr','-o','energies.xvg'],'energy','Potential\nTemperature\n0\n')
    v=np.array([[float(x) for x in l.split()] for l in (out/'energies.xvg').read_text().splitlines() if l and not l.startswith(('#','@'))])
    if not np.isfinite(v).all() or v[-1,0]-v[0,0]<0.199 or 'LINCS WARNING' in (out/'md.log').read_text():raise RuntimeError('Pilot energy/duration/constraint gate failed')
    result={'state':'unglycosylated_FXII_surface_initialization_passed','duration_ps':0.2,'protein_present':True,'protein_model':'predicted mature FXII, glycans absent','box_nm':[18,12,18],'total_atoms':int((out/'md.gro').read_text().splitlines()[1]),**spec,'salt':'neutralizing ions only','limitations':['initialization only, not adsorption, equilibrium, native glycoform or activation','correlated tiled provisional coating/PC;one rigid orientation per pilot']}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
