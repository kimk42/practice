"""After expanded EM, build and test one unglycosylated FXII surface pilot.
This numerical pilot is not native adsorption, equilibration or activation MD.
"""
import json,subprocess,time,shutil,re
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    surface=ROOT/'runs/surface-expanded-18x12';poses=ROOT/'runs/fxii-surface-poses';out=ROOT/'runs/fxii-surface-pilot'
    state=ROOT/'runs/fxii-surface-pilot-status.json'
    def status(value):state.write_text(json.dumps(value,indent=2)+'\n')
    deadline=time.monotonic()+43200
    while not (surface/'result.json').exists():
        if time.monotonic()>deadline:raise RuntimeError('Expanded EM wait limit reached; do not bypass its gate')
        status({'state':'waiting_for_expanded_surface_EM','scope':'unglycosylated predicted FXII numerical pilot'})
        time.sleep(15)
    report=json.loads((surface/'result.json').read_text())
    if report['state']!='expanded solvated candidate minimized':raise ValueError('Expanded EM not accepted')
    if out.exists():raise RuntimeError('Existing protein-surface pilot preserved')
    out.mkdir();status({'state':'building_unglycosylated_pilot'})
    def run(args,name,stdin=None):
        with (out/(name+'-launch.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['trjconv','-s',str(surface/'em.tpr'),'-f',str(surface/'em.gro'),'-o','whole-surface.gro','-pbc','mol'],'whole','0\n')
    lines=(out/'whole-surface.gro').read_text().splitlines();n=report['dry_atoms'];dry=lines[2:2+n];water=lines[2+n:-1]
    xyz=np.array([[float(s[20:28]),float(s[28:36]),float(s[36:44])] for s in dry])
    candidate=next(x for x in json.loads((poses/'poses.json').read_text())['poses'] if x['fits_lateral_padding_1_2nm'])
    protein=(poses/candidate['file']).read_text().splitlines()[2:-1]
    pxyz=np.array([[float(s[20:28]),float(s[28:36]),float(s[36:44])] for s in protein])
    pxyz[:,:2]+=[9.,6.]-(pxyz[:,:2].min(axis=0)+pxyz[:,:2].max(axis=0))/2
    pxyz[:,2]+=xyz[:,2].max()+2.-pxyz[:,2].min()
    # Check actual rounded coordinates, all atoms, before accepting placement.
    if np.any(pxyz.min(axis=0)<1.2) or np.any(pxyz.max(axis=0)>np.array([18.,12.,18.])-1.2):raise ValueError('Protein does not fit requested periodic clearance')
    combined=list(dry)
    for index,(row,(x,y,z)) in enumerate(zip(protein,pxyz),start=n+1):combined.append(f'{row[:15]}{index%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}')
    box='18.00000 12.00000 18.00000\n'
    (out/'dry.gro').write_text('Unglycosylated predicted FXII above provisional coating\n'+str(len(combined))+'\n'+'\n'.join(combined)+'\n'+box)
    (out/'water-box.gro').write_text('Existing expanded surface waters\n'+str(len(water))+'\n'+'\n'.join(water)+'\n'+box)
    for file in ['GRAFT.itp','PC4.itp','gaff2-types.itp']:shutil.copy2(surface/file,out/file)
    protein_source=ROOT/'runs/fxii-unglycosylated-topology';top=(protein_source/'topol.top').read_text()
    body=top[top.index('[ moleculetype ]'):top.index('; Include water topology')]
    (out/'protein.itp').write_text(body);shutil.copy2(protein_source/'posre.itp',out/'posre.itp')
    st=(surface/'topol.top').read_text();prefix,molecules=st.split('[ molecules ]');molecule_rows=[line for line in molecules.splitlines() if line.strip() and not line.lstrip().startswith(';') and not line.split()[0]=='SOL']
    # All molecule definitions must precede [ system ].
    system=prefix.index('[ system ]');prefix=prefix[:system]+'#include "protein.itp"\n#include "amber99sb-ildn.ff/ions.itp"\n'+prefix[system:]
    (out/'topol.top').write_text(prefix+'[ molecules ]\n'+'\n'.join(molecule_rows)+'\nProtein_chain_A 1\n')
    run(['solvate','-cp','dry.gro','-cs','water-box.gro','-p','topol.top','-o','solvated.gro'],'solvate')
    common='integrator = steep\nnsteps = 5000\nemtol = 1000\nemstep = 0.005\ncutoff-scheme = Verlet\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n'
    (out/'ions.mdp').write_text(common+'coulombtype = Cut-off\n')
    run(['grompp','-f','ions.mdp','-c','solvated.gro','-p','topol.top','-o','ions.tpr'],'ions-prep')
    run(['genion','-s','ions.tpr','-o','neutralized.gro','-p','topol.top','-pname','NA','-nname','CL','-neutral','-seed','2026'],'ions','SOL\n')
    (out/'em.mdp').write_text(common+'coulombtype = PME\n')
    run(['grompp','-f','em.mdp','-c','neutralized.gro','-p','topol.top','-o','em.tpr'],'em-prep')
    status({'state':'minimizing_unglycosylated_FXII_surface_pilot','pose':candidate,'minimum_initial_vertical_separation_nm':2.0})
    run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','4'],'em')
    if 'converged to Fmax <' not in (out/'em.log').read_text():raise RuntimeError('Combined EM threshold not met')
    (out/'md.mdp').write_text('integrator = md\ndt = 0.0005\nnsteps = 400\nconstraints = h-bonds\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\nnstenergy = 20\npbc = xyz\n')
    run(['grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'md-prep')
    status({'state':'running_0_2ps_numerical_pilot'})
    run(['mdrun','-deffnm','md','-ntmpi','1','-ntomp','4'],'md')
    run(['check','-e','md.edr'],'energy-check')
    run(['energy','-f','md.edr','-o','energies.xvg'],'energy','Potential\nTemperature\n0\n')
    values=np.array([[float(x) for x in line.split()] for line in (out/'energies.xvg').read_text().splitlines() if line and not line.startswith(('#','@'))])
    if not np.isfinite(values).all() or values[-1,0]-values[0,0]<0.199:raise RuntimeError('Incomplete/nonfinite energy record')
    if 'LINCS WARNING' in (out/'md.log').read_text():raise RuntimeError('Constraint warnings detected')
    result={'state':'unglycosylated_FXII_surface_initialization_passed','duration_ps':0.2,'protein_present':True,'protein_model':'predicted mature FXII, glycans absent','box_nm':[18,12,18],'total_atoms':int((out/'md.gro').read_text().splitlines()[1]),'initial_pose':candidate,'salt':'neutralizing ions only','limitations':['initialization only; not adsorption, equilibrium, native-glycoform validation or FXII activation','correlated tiled low-density short-PC coating','one rigid initial orientation sampled; others prepared but not simulated']}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');status(result)
if __name__=='__main__':
    try:main()
    except Exception as e:
        (ROOT/'runs/fxii-surface-pilot-status.json').write_text(json.dumps({'state':'blocked','reason':str(e)},indent=2)+'\n');raise
