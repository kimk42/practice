"""Audited refined-citrate topology and isolated neutral-water numerical checks.

Two fixed protonation candidates, no assigned pH population. No existing system
is modified and no cPBSzI mixture/protein simulation is launched.
"""
from decimal import Decimal
import fcntl
import json
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
from rdkit import Chem
from audit_citrate_representation import review
from pc_length_parameters import ROOT,sha
from run_pc_sqm_checkpointed import amber_environment,write_json
from surface_validation import section_rows

BASE=ROOT/'runs/citrate-refined-controls-v2'


def execute(species,name,charge,lock_fd):
    qm=ROOT/'runs/citrate-gradient-refinement-v2'/species
    completed=json.loads((qm/'completed.json').read_text())
    source=qm/'bcc-compatible-v2/candidate.mol2'
    evidence=json.loads((qm/'candidate-audit.json').read_text())
    if sha(source)!=evidence['candidate_mol2_sha256']:raise ValueError('Source charge candidate changed')
    audit=review(species,source,ROOT/completed['final_output'])
    if not audit['QM_log_review']['final_gradient_numerically_certified']:raise ValueError('Missing final gradient gate')
    out=BASE/name
    if (out/'result.json').exists():
        result=json.loads((out/'result.json').read_text())
        if sha(out/'md.cpt')!=result['checkpoint_sha256']:raise ValueError('Completed checkpoint changed')
        return result
    if out.exists():raise RuntimeError('Partial control retained;inspect checkpoint/log before recovery')
    out.mkdir();shutil.copy2(source,out/'input.mol2');write_json(out/'representation-audit.json',audit)
    env=amber_environment()
    def run(command,label,stdin=None):
        with (out/(label+'.log')).open('x') as log:
            subprocess.run(command,cwd=out,env=env,input=stdin,text=True,stdout=log,stderr=subprocess.STDOUT,
                           check=True,timeout=1200,pass_fds=(lock_fd,))
    run(['parmchk2','-i','input.mol2','-f','mol2','-o','candidate.frcmod','-s','gaff2'],'parmchk')
    if re.search(r'ATTN|need\s+revision',(out/'candidate.frcmod').read_text(),re.I):raise ValueError('Unassigned bonded terms')
    (out/'leap.in').write_text('source leaprc.gaff2\nloadAmberParams candidate.frcmod\nM = loadMol2 input.mol2\n'
                             'check M\nsaveAmberParm M candidate.prmtop candidate.inpcrd\nquit\n')
    run(['tleap','-f','leap.in'],'tleap')
    leap=(out/'tleap.log').read_text()
    if 'Errors = 0; Warnings = 2;' not in leap or leap.count(f'The unperturbed charge of the unit ({charge:.6f}) is not zero.')!=2 or leap.count('Warning!')!=2:
        raise ValueError('TLEaP output differs from the two expected isolated-anion charge warnings')
    run(['acpype','-p','candidate.prmtop','-x','candidate.inpcrd','-b',name,'-a','gaff2'],'acpype')
    text=(out/f'{name}.amb2gmx/{name}_GMX.top').read_text()
    types=text[text.index('[ atomtypes ]'):text.index('[ moleculetype ]')]
    molecule=text[text.index('[ moleculetype ]'):text.index('[ system ]')]
    mol2=source.read_text();original=[r.split() for r in mol2.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0].splitlines() if r.strip()]
    if sum(Decimal(r[8]) for r in original)!=Decimal(charge):raise ValueError('Source is not exactly integer charge;no automatic correction')
    topatoms=list(section_rows(molecule,'atoms'))
    if [r[4] for r in topatoms]!=[r[1] for r in original]:raise ValueError('Topology atom names/order differ')
    if any(Decimal(a[6])!=Decimal(b[8]) for a,b in zip(topatoms,original)):raise ValueError('Converter changed citrate charges')
    if {tuple(sorted(map(int,r[:2]))) for r in section_rows(molecule,'bonds')}!={tuple(sorted(map(int,r[1:3]))) for r in [s.split() for s in mol2.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0].splitlines() if s.strip()]}:
        raise ValueError('Topology graph differs')
    renamed=[]
    for line in types.splitlines():
        row=line.split(';',1)[0].split()
        if row and not line.lstrip().startswith(('[',';')):
            row[0]='CBUF_'+row[0];row[1]='CBUF_'+row[1];line=' '.join(row)
        renamed.append(line)
    (out/'types.itp').write_text('\n'.join(renamed)+'\n')
    section=None;renamed=[]
    for line in molecule.splitlines():
        row=line.split(';',1)[0].split()
        if line.startswith('['):section=line.split(']')[0].strip('[ ')
        elif row and section=='atoms':row[1]='CBUF_'+row[1];line=' '.join(row)
        renamed.append(line)
    (out/'citrate.itp').write_text('\n'.join(renamed)+'\n')
    molname=list(section_rows(molecule,'moleculetype'))[0][0]
    gro=(out/f'{name}.amb2gmx/{name}_GMX.gro').read_text().splitlines();raw=gro[2:-1]
    xyz=np.array([[float(s[j:j+8]) for j in (20,28,36)] for s in raw]);xyz+=1-xyz.mean(axis=0)
    rows=[f'{s[:20]}{x:8.3f}{y:8.3f}{z:8.3f}' for s,(x,y,z) in zip(raw,xyz)]
    for i,(x,y,z) in enumerate([(2.4,2.4,2.4),(.4,2.4,.4),(2.4,.4,.4)][: -charge]):
        rows.append(f'{i+2:5d}{"NAJC":<5}{"JC_NA":>5}{len(raw)+i+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}')
    (out/'ions.gro').write_text('Refined citrate numerical control\n'+str(len(rows))+'\n'+'\n'.join(rows)+'\n3 3 3\n')
    (out/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "types.itp"\n'
        '#include "'+str(ROOT/'models/buffer-parameters/jc-tip3p-monovalent.itp')+'"\n#include "citrate.itp"\n'
        '#include "amber99sb-ildn.ff/tip3p.itp"\n'+f'[ system ]\nRefined {name} numerical control\n[ molecules ]\n{molname} 1\nNAJC {-charge}\n')
    run(['gmx','solvate','-cp','ions.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
    common='cutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1\nrvdw = 1\npbc = xyz\n'
    (out/'em.mdp').write_text(common+'integrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.002\n')
    run(['gmx','grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],'em-prep')
    run(['gmx','mdrun','-deffnm','em','-ntmpi','1','-ntomp','1'],'em-launch')
    em=(out/'em.log').read_text()
    if 'converged to Fmax <' not in em:raise ValueError('EM did not meet force criterion')
    (out/'md.mdp').write_text(common+'integrator = md\ndt = 0.0005\nnsteps = 4000\nconstraints = h-bonds\n'
        'tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\n'
        'gen-seed = 2026\nnstenergy = 100\nnstxout = 100\nnstxout-compressed = 100\nnstlog = 1000\n')
    run(['gmx','grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'md-prep')
    run(['gmx','mdrun','-deffnm','md','-ntmpi','1','-ntomp','1','-cpt','1'],'md-launch')
    run(['gmx','energy','-f','md.edr','-o','energy.xvg'],'energy','Potential\nTemperature\n0\n')
    data=np.array([[float(v) for v in s.split()] for s in (out/'energy.xvg').read_text().splitlines() if s and not s.startswith(('#','@'))])
    if data.shape!=(41,3) or not np.isfinite(data).all() or data[-1,0]!=2:raise ValueError('Missing/nonfinite energy frames')
    if 'Finished mdrun' not in (out/'md.log').read_text() or 'LINCS WARNING' in (out/'md.log').read_text():raise ValueError('MD completion/constraints failed')
    for stage in ('em','md'):
        if re.search(r'^WARNING',(out/(stage+'-prep.log')).read_text(),re.M):raise ValueError('GROMACS preprocessing warning')
    n=len(original);(out/'solute.ndx').write_text('[ Solute ]\n'+' '.join(map(str,range(1,n+1)))+'\n')
    run(['gmx','traj','-f','md.trr','-s','md.tpr','-n','solute.ndx','-ox','solute-xyz.xvg','-xvg','none'],'geometry','0\n')
    positions=np.loadtxt(out/'solute-xyz.xvg')
    if positions.shape!=(41,1+3*n):raise ValueError('Unexpected coordinate-frame count')
    xyz=positions[:,1:].reshape(41,n,3)
    formal=Chem.MolFromMolFile(str(ROOT/f'models/buffer-parameters/{species}.mol'),removeHs=False)
    lengths=[]
    for bond in formal.GetBonds():
        delta=xyz[:,bond.GetBeginAtomIdx()]-xyz[:,bond.GetEndAtomIdx()];delta-=3*np.rint(delta/3)
        lengths.extend(np.linalg.norm(delta,axis=1).tolist())
    if not np.isfinite(lengths).all() or min(lengths)<.07 or max(lengths)>.21:raise ValueError('Bond geometry outside numerical envelope')
    result={'species':species,'state':'refined candidate neutral-water EM and2ps numerical control completed',
            'duration_ps':2.,'energy_frames':41,'geometry_frames':41,'temperature_mean_K':float(data[:,2].mean()),
            'EM_Fmax_kJ_mol_nm':float(re.search(r'Maximum force\s*=\s*([\deE+.-]+)',em).group(1)),
            'bond_length_range_nm':[min(lengths),max(lengths)],'solute_charge_e':charge,'sodium_ions':-charge,
            'nominal_solute_concentration_mM':1/6.02214076e23/(27e-24)*1000,
            'charge_conversion':'verbatim source charges retained;no projection or hidden ACPYPE balancing',
            'warnings_reviewed':'only two expected isolated-anion TLEaP charge warnings;strict GROMACS without warning override',
            'source_mol2_sha256':sha(source),'checkpoint_sha256':sha(out/'md.cpt'),
            'itp_sha256':sha(out/'citrate.itp'),'types_sha256':sha(out/'types.itp'),'topology_sha256':sha(out/'topol.top'),
            'scope':'fixed microstate,300K,small TIP3P/JC-Na neutral box;concentration differs from10mM citrate in cPBSzI',
            'physical_approval':False,'full_buffer_ready':False,'run_directory':str(out.relative_to(ROOT))}
    write_json(out/'result.json',result);print(json.dumps(result),flush=True);return result


if __name__=='__main__':
    BASE.mkdir(exist_ok=True)
    with (BASE/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        results=[execute(s,n,q,lock.fileno()) for s,n,q in [('citrate3-','CIT3',-3),('hydrogen_citrate2-','HCIT2',-2)]]
        write_json(ROOT/'models/buffer-parameters/refined-citrate-water-audit.json',{'species':results,'physical_approval':False,'full_buffer_ready':False})
