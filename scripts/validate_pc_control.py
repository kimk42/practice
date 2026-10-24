"""Run a numerical GAFF2/AM1-BCC PC4-in-water control, not a surface simulation."""
import json
from pathlib import Path
import shutil
import sys
import subprocess
sys.path.insert(0,str(Path(__file__).resolve().parent))
from surface_validation import parameter_gate
ROOT=Path(__file__).resolve().parents[1]

def main():
    item=parameter_gate('PC4')
    if item is None:raise RuntimeError('PC4 parameters not ready')
    out=ROOT/'runs/pc4-water-control'
    if out.exists():raise RuntimeError('Existing control preserved; inspect it before any rerun')
    out.mkdir()
    shutil.copy2(item['itp'],out/'PC4.itp')
    (out/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "PC4.itp"\n#include "amber99sb-ildn.ff/tip3p.itp"\n[ system ]\nExploratory PC4 water numerical control\n[ molecules ]\n'+item['molecule']+' 1\n')
    def run(args,name):
        with (out/(name+'.log')).open('w') as log:subprocess.run(args,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['gmx','editconf','-f',str(item['gro']),'-o','dry.gro','-d','1.2','-bt','cubic'],'box')
    run(['gmx','solvate','-cp','dry.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
    (out/'em.mdp').write_text('integrator = steep\nnsteps = 5000\nemtol = 1000\nemstep = 0.005\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
    run(['gmx','grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],'em-prep')
    run(['gmx','mdrun','-deffnm','em','-ntmpi','1','-ntomp','2'],'em-run')
    if 'converged to Fmax <' not in (out/'em.log').read_text():raise RuntimeError('EM force threshold not met')
    (out/'md.mdp').write_text('integrator = md\ndt = 0.0005\nnsteps = 400\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\nnstenergy = 20\npbc = xyz\n')
    run(['gmx','grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'md-prep')
    run(['gmx','mdrun','-deffnm','md','-ntmpi','1','-ntomp','2'],'md-run')
    run(['gmx','check','-e','md.edr'],'energy-check')
    if not (out/'md.gro').exists():raise RuntimeError('Missing final coordinates')
    report={'state':'PC4 water numerical control passed','PC4_atoms':len(item['atoms']),'net_PC4_charge':item['charge'],
            'EM':'configured force threshold met','MD_steps':400,'MD_dt_ps':0.0005,'MD_duration_ps':0.2,
            'temperature_K':300,'scope':'one PC4 oligomer in water; not the coating, substrate slab, adsorption or full proteins',
            'limitations':'Numerical initialization only, not force-field accuracy or equilibrium validation'}
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
