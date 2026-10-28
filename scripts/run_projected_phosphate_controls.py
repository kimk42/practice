"""Strict isolated ion/water controls of projected-charge GAFF2 candidates."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    reports=[]
    for species,name,target in [('dihydrogen_phosphate-','PH1',-1),('hydrogen_phosphate2-','PH2',-2)]:
        param=ROOT/'runs/phosphate-projected-parameterization-v1'/name
        out=ROOT/'runs'/f'{name.lower()}-projected-water-control'
        if out.exists():raise RuntimeError('Existing control preserved:'+str(out))
        leap=(param/'tleap.log').read_text()
        assert 'Errors = 0; Warnings = 2;' in leap
        assert leap.count(f'The unperturbed charge of the unit ({target:.6f}) is not zero.')==2
        assert 'Warning!' not in leap.replace('Warning!','',2)
        # The two warnings refer solely to the deliberately charged single ion.
        # Its MD box is explicitly neutralized below;no grompp warning bypass.
        mol2=(param/'input.mol2').read_text();atoms=[x.split() for x in mol2.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0].splitlines() if x.strip()]
        charge_by_name={x[1]:x[8] for x in atoms}
        assert sum(Decimal(q) for q in charge_by_name.values())==Decimal(target)
        original=(param/f'{name}.amb2gmx/{name}_GMX.top').read_text()
        types=original[original.index('[ atomtypes ]'):original.index('[ moleculetype ]')]
        molecule=original[original.index('[ moleculetype ]'):original.index('[ system ]')]
        out.mkdir();newtypes=[]
        for line in types.splitlines():
            row=line.split(';',1)[0].split()
            if row and not line.lstrip().startswith(('[',';')):
                row[0]='PBUF_'+row[0];row[1]='PBUF_'+row[1];line=' '.join(row)
            newtypes.append(line)
        section=None;newmol=[];seen=[]
        for line in molecule.splitlines():
            row=line.split(';',1)[0].split()
            if line.startswith('['):section=line.split(']')[0].strip('[ ')
            elif row and section=='atoms':
                assert row[4] in charge_by_name
                row[1]='PBUF_'+row[1];row[6]=charge_by_name[row[4]]
                seen.append(row[4]);line=' '.join(row)
            newmol.append(line)
        assert seen==[a[1] for a in atoms]
        # Bond graph must agree with the preserved typed input;do not infer new bonds.
        block=molecule.split('[ bonds ]')[1].split('[ pairs ]')[0]
        bonded={tuple(sorted(map(int,s.split()[:2]))) for s in block.splitlines() if s.strip() and not s.lstrip().startswith(';')}
        expected={tuple(sorted(map(int,s.split()[1:3]))) for s in mol2.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0].splitlines() if s.strip()}
        assert bonded==expected
        (out/'types.itp').write_text('\n'.join(newtypes)+'\n');(out/'phosphate.itp').write_text('\n'.join(newmol)+'\n')
        gro=(param/f'{name}.amb2gmx/{name}_GMX.gro').read_text().splitlines();raw=gro[2:-1]
        xyz=np.array([[float(s[20:28]),float(s[28:36]),float(s[36:44])] for s in raw]);xyz+=1-xyz.mean(axis=0)
        rows=[f'{s[:20]}{x:8.3f}{y:8.3f}{z:8.3f}' for s,(x,y,z) in zip(raw,xyz)]
        for i in range(-target):
            rows.append(f'{i+2:5d}{"NAJC":<5}{"JC_NA":>5}{len(raw)+i+1:5d}{2.:8.3f}{2.:8.3f}{2.-1.5*i:8.3f}')
        (out/'ions.gro').write_text('Projected phosphate numerical control\n'+str(len(rows))+'\n'+'\n'.join(rows)+'\n3 3 3\n')
        (out/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "types.itp"\n'
            '#include "'+str(ROOT/'models/buffer-parameters/jc-tip3p-monovalent.itp')+'"\n'
            '#include "phosphate.itp"\n#include "amber99sb-ildn.ff/tip3p.itp"\n'
            f'[ system ]\nProjected {name} numerical control\n[ molecules ]\n{name} 1\nNAJC {-target}\n')
        def run(args,log,stdin=None):
            with (out/(log+'.log')).open('w') as stream:
                subprocess.run(['gmx',*args],cwd=out,input=stdin,text=True,stdout=stream,stderr=subprocess.STDOUT,check=True)
        common='cutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1\nrvdw = 1\npbc = xyz\n'
        run(['solvate','-cp','ions.gro','-cs','spc216.gro','-p','topol.top','-o','solvated.gro'],'solvate')
        (out/'em.mdp').write_text(common+'integrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.002\n')
        run(['grompp','-f','em.mdp','-c','solvated.gro','-p','topol.top','-o','em.tpr'],'em-prep')
        run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','1'],'em-launch')
        em=(out/'em.log').read_text();assert 'converged to Fmax <' in em
        (out/'md.mdp').write_text(common+'integrator = md\ndt = 0.0005\nnsteps = 4000\nconstraints = h-bonds\ntcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\npcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\nnstenergy = 100\nnstxout = 100\nnstxout-compressed = 100\nnstlog = 1000\n')
        run(['grompp','-f','md.mdp','-c','em.gro','-p','topol.top','-o','md.tpr'],'md-prep')
        run(['mdrun','-deffnm','md','-ntmpi','1','-ntomp','1','-cpt','1'],'md-launch')
        run(['energy','-f','md.edr','-o','energy.xvg'],'energy','Potential\nTemperature\n0\n')
        energy=np.array([[float(v) for v in s.split()] for s in (out/'energy.xvg').read_text().splitlines() if s and not s.startswith(('#','@'))])
        assert energy.shape==(41,3) and np.isfinite(energy).all() and energy[-1,0]==2
        assert 'Finished mdrun' in (out/'md.log').read_text() and 'LINCS WARNING' not in (out/'md.log').read_text()
        for stage in ['em','md']:assert not re.search(r'^WARNING',(out/(stage+'-prep.log')).read_text(),re.M)
        report={'species':species,'state':'neutral ion-water EM and2ps numerical control passed',
                'duration_ps':2,'finite_energy_frames':len(energy),'temperature_mean_K':float(energy[:,2].mean()),
                'EM_Fmax_kJ_mol_nm':float(re.search(r'Maximum force\s*=\s*([\deE+.-]+)',em).group(1)),
                'net_solute_charge_e':target,'sodium_ions':-target,
                'charge_precision_policy':'GROMACS atom charges restored verbatim from projected mol2 after matching names/order;no second charge fit',
                'typing_policy':'previously typed GAFF2 p5/o/oh/ho topology preserved;nonprotonated O are equivalent resonance sites',
                'parameter_generation':'parmchk2 -s gaff2;tleap leaprc.gaff2;ACPYPE amb2gmx;source files preserved',
                'tleap_warnings_reviewed':'two warnings that isolated anion is not neutral;expected integer charge and explicitly neutral MD box',
                'source_mol2_sha256':sha(param/'input.mol2'),'checkpoint_sha256':sha(out/'md.cpt'),
                'prmtop_sha256':sha(param/'candidate.prmtop'),'types_sha256':sha(out/'types.itp'),
                'itp_sha256':sha(out/'phosphate.itp'),'run_directory':str(out.relative_to(ROOT)),
                'physical_force_field_approved':False,'full_buffer_ready':False}
        (out/'result.json').write_text(json.dumps(report,indent=2)+'\n');reports.append(report)
    (ROOT/'models/buffer-parameters/projected-phosphate-water-audit.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))


if __name__=='__main__':main()
