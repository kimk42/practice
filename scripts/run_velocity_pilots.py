"""Three velocity seeds per accepted pilot geometry; not production replicas."""
import json,subprocess,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs'
def main():
    root=RUN/'velocity-pilots';root.mkdir(exist_ok=True);records=[]
    for case,surface in [('fxii-coating-near-pilot','surface-expanded-18x12'),('fxii-bare-near-pilot','bare-pc-matched')]:
        if not (RUN/case/'result.json').exists():subprocess.run(['python',str(ROOT/'scripts/run_fxii_surface_pilot.py'),'--surface-dir',str(RUN/surface),'--output-dir',str(RUN/case),'--gap-nm','0.5'],cwd=ROOT,check=True)
    cases=['fxii-surface-pilot','fxii-bare-pc-pilot','fxii-coating-pose3-pilot','fxii-bare-pc-pose3-pilot','fxii-coating-near-pilot','fxii-bare-near-pilot']
    for case in cases:
        source=RUN/case
        if not (source/'result.json').exists():raise RuntimeError('Source pilot incomplete: '+case)
        records.append({'case':case,'seed':2026,'source_result':str((source/'result.json').relative_to(ROOT))})
        for seed in [2027,2028]:
            out=root/f'{case}-{seed}';out.mkdir(exist_ok=True)
            if not (out/'result.json').exists():
                def run(args,name,stdin=None):
                    with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
                if not (out/'md.tpr').exists():
                    for f in ['topol.top','em.gro',*[p.name for p in source.glob('*.itp')]]:shutil.copy2(source/f,out/f)
                    mdp=(source/'md.mdp').read_text()
                    import re
                    mdp=re.sub(r'gen-seed\s*=\s*\d+',f'gen-seed = {seed}',mdp)
                    (out/'md.mdp').write_text(mdp)
                    run(['grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'prep')
                args=['mdrun','-deffnm','md','-cpt','1','-ntmpi','1','-ntomp','4']
                if (out/'md.cpt').exists():args+=['-cpi','md.cpt','-append']
                run(args,'md');run(['check','-e','md.edr'],'check');run(['energy','-f','md.edr','-o','energies.xvg'],'energy','Potential\nTemperature\n0\n')
                v=np.array([[float(x) for x in l.split()] for l in (out/'energies.xvg').read_text().splitlines() if l and not l.startswith(('#','@'))])
                if not np.isfinite(v).all() or v[-1,0]-v[0,0]<0.199 or 'LINCS WARNING' in (out/'md.log').read_text():raise RuntimeError('Velocity pilot did not pass')
                (out/'result.json').write_text(json.dumps({'state':'velocity numerical pilot passed','case':case,'seed':seed,'duration_ps':0.2,'scope':'same minimized geometry, changed velocities; not independent conformational production sampling'},indent=2)+'\n')
            records.append({'case':case,'seed':seed,'source_result':str((out/'result.json').relative_to(ROOT))})
    (root/'result.json').write_text(json.dumps({'state':'eighteen short velocity pilots passed','records':records,'scope':'numerical sensitivity only, not replicated adsorption evidence'},indent=2)+'\n')
if __name__=='__main__':main()
