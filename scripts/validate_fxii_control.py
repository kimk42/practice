"""Unglycosylated predicted FXII numerical pilot; not native adsorption MD."""
import json,re,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/fxii-unglycosylated-topology'
def main():
    audit=json.loads((ROOT/'models/proteins/fxii-audit.json').read_text())
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'em.tpr').exists():raise RuntimeError('Existing dynamics preserved; inspect before rerun')
    def run(args,name,stdin=None):
        with (OUT/(name+'.log')).open('w') as log:
            subprocess.run(['gmx',*args],cwd=OUT,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['pdb2gmx','-f',str(ROOT/audit['candidate_file']),'-o','processed.gro','-p','topol.top','-ff','amber99sb-ildn','-water','tip3p','-ignh'],'preparation')
    links={tuple(sorted(map(int,pair))) for pair in re.findall(r'Linking CYS-(\d+) SG-\d+ and CYS-(\d+) SG-',(OUT/'preparation.log').read_text())}
    expected={tuple(x['residues']) for x in audit['disulfide_annotations']}
    if links!=expected:raise RuntimeError('Disulfide mismatch')
    run(['editconf','-f','processed.gro','-o','dry.gro','-d','1.2','-bt','cubic'],'box')
    run(['solvate','-cp','dry.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
    # This cutoff-only TPR is solely for neutralizing-ion placement, never dynamics.
    common='integrator = steep\nnsteps = 5000\nemtol = 1000\nemstep = 0.005\ncutoff-scheme = Verlet\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n'
    (OUT/'ions.mdp').write_text(common+'coulombtype = Cut-off\n')
    run(['grompp','-f','ions.mdp','-c','solvated.gro','-p','topol.top','-o','ions.tpr'],'ions-prep')
    run(['genion','-s','ions.tpr','-o','neutralized.gro','-p','topol.top','-pname','NA','-nname','CL','-neutral'],'ions','SOL\n')
    (OUT/'em.mdp').write_text(common+'coulombtype = PME\n')
    run(['grompp','-f','em.mdp','-c','neutralized.gro','-p','topol.top','-o','em.tpr'],'em-prep')
    run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','2'],'em-run')
    if 'converged to Fmax <' not in (OUT/'em.log').read_text():raise RuntimeError('EM force threshold not met')
    (OUT/'md.mdp').write_text('integrator = md\ndt = 0.0005\nnsteps = 400\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\nnstenergy = 20\npbc = xyz\n')
    run(['grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'md-prep')
    run(['mdrun','-deffnm','md','-ntmpi','1','-ntomp','2'],'md-run')
    run(['check','-e','md.edr'],'energy-check')
    report={'state':'unglycosylated predicted FXII numerical pilot passed','disulfides_verified':len(links),'duration_ps':0.2,'temperature_K':300,'salt':'neutralizing ions only, not physiological saline','forcefield':'amber99sb-ildn / TIP3P','protonation':'pdb2gmx standard ionization and automatic histidine assignment; provisional','scope':'protein alone, predicted conformation, glycans absent; not native validation or surface adsorption'}
    (OUT/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
if __name__=='__main__':main()
