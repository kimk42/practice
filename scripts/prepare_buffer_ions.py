"""Convert published JC/TIP3P monovalent parameters without altering active baselines."""
import json,hashlib,math,subprocess
from pathlib import Path
from acpype.utils import bundled_amber_dir
ROOT=Path(__file__).resolve().parents[1]
def main():
    if (ROOT/'runs/jc-ion-water-control').exists():raise RuntimeError('Existing ion control preserved; no output overwritten')
    source=Path(bundled_amber_dir())/'dat/leap/parm/frcmod.ionsjc_tip3p';data={};section=None
    selected={'Na+':('JC_NA','NAJC',11,1),'Cl-':('JC_CL','CLJC',17,-1),'I-':('JC_I','IJC',53,-1)}
    for line in source.read_text().splitlines():
        row=line.split()
        if not row:continue
        if row[0] in ['MASS','NONBON']:section=row[0];continue
        if row[0] not in selected:continue
        entry=data.setdefault(row[0],{})
        if section=='MASS':entry['mass']=float(row[1])
        if section=='NONBON':entry.update({'Rmin_half_A':float(row[1]),'epsilon_kcal_mol':float(row[2])})
    out=ROOT/'models/buffer-parameters';out.mkdir(exist_ok=True);itp='; Joung & Cheatham2008 Ewald/TIP3P,from Amber bundled frcmod\n[ atomtypes ]\n';molecules=''
    for name,(atom,molecule,z,charge) in selected.items():
        e=data[name];sigma=2*e['Rmin_half_A']*0.1/2**(1/6);epsilon=e['epsilon_kcal_mol']*4.184;e.update({'sigma_nm':sigma,'epsilon_kJ_mol':epsilon,'formal_charge':charge,'atomtype':atom,'molecule':molecule})
        rmin=2*e['Rmin_half_A']*.1;energy=4*epsilon*((sigma/rmin)**12-(sigma/rmin)**6)
        if abs(energy+epsilon)>1e-10:raise ValueError('LJ conversion failed')
        itp+=f'{atom} {z} {e["mass"]} 0 A {sigma:.10f} {epsilon:.10f}\n'
        molecules+=f'\n[ moleculetype ]\n{molecule} 1\n[ atoms ]\n1 {atom} 1 {molecule} {atom} 1 {charge} {e["mass"]}\n'
    (out/'jc-tip3p-monovalent.itp').write_text(itp+molecules)
    proof={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'citation':'Joung and Cheatham,J Phys Chem B2008,monovalent Ewald/TIP3P parameters','species':data,'conversion':'sigma_nm=2*(Rmin/2 in A)*0.1/2^(1/6);epsilon_kJ=epsilon_kcal*4.184','scope':'consistent Na/Cl/I proposal for future buffer;not substituted into ongoing legacy-ion numerical baselines','full_buffer_ready':False}
    (out/'provenance.json').write_text(json.dumps(proof,indent=2)+'\n')
    control=ROOT/'runs/jc-ion-water-control'
    if control.exists():raise RuntimeError('Existing ion control preserved')
    control.mkdir();rows=[]
    for i,(res,atom,q) in enumerate([('NAJC','JC_NA',(.4,.4,.4)),('NAJC','JC_NA',(1.6,1.6,1.6)),('CLJC','JC_CL',(1.4,.4,1.4)),('IJC','JC_I',(.4,1.4,.6))],1):
        rows.append(f'{i:5d}{res:<5}{atom:>5}{i:5d}{q[0]:8.3f}{q[1]:8.3f}{q[2]:8.3f}')
    (control/'ions.gro').write_text('JC ion import numerical control\n4\n'+'\n'.join(rows)+'\n3.0 3.0 3.0\n')
    (control/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "'+str(out/'jc-tip3p-monovalent.itp')+'"\n#include "amber99sb-ildn.ff/tip3p.itp"\n[ system ]\nJC ion numerical control\n[ molecules ]\nNAJC 2\nCLJC 1\nIJC 1\n')
    def run(args,name):
        with (control/(name+'.log')).open('w') as log:subprocess.run(['gmx',*args],cwd=control,stdout=log,stderr=subprocess.STDOUT,check=True)
    run(['solvate','-cp','ions.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
    (control/'em.mdp').write_text('integrator = steep\nnsteps = 5000\nemtol = 1000\nemstep = 0.005\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\npbc = xyz\n')
    run(['grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],'prep');run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','1'],'run')
    if 'converged to Fmax <' not in (control/'em.log').read_text():raise RuntimeError('Ion-control EM criterion not met')
    proof['numerical_control']='neutral2Na/1Cl/1I water import and EM passed;not physical solvation validation or experimental buffer';(out/'provenance.json').write_text(json.dumps(proof,indent=2)+'\n');print(proof['numerical_control'])
if __name__=='__main__':main()
