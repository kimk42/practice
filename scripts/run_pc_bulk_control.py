"""Strict EM,2ps NVT and20ps NPT dry-PC numerical controls with MD checkpoints.

No equilibrium-PC or experimental-density claim follows from this short run.
300K,1bar and barostat compressibility are numerical protocol assumptions.
"""
import argparse
import fcntl
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from pc_length_parameters import ROOT, sha
from surface_validation import section_rows
from run_pc_sqm_checkpointed import write_json
from run_pc_packing_margin import accepted_folder

BASE = ROOT/'runs/pc-bulk-controls-v2'


def execute(case_name):
    packed = accepted_folder(case_name)
    if not (packed/'geometry-audit.json').exists():
        return None
    case=json.loads((packed/'case.json').read_text()); dp=case['DP']
    param=ROOT/f'runs/pc-length-forcefield-v2/pc{dp}'
    if not (param/'parameter-audit.json').exists():
        return None
    evidence=json.loads((param/'parameter-audit.json').read_text())
    if evidence['state']!='numerical-control candidate;physical forcefield not approved':
        raise ValueError('Forcefield has not passed numerical assignment review')
    folder=ROOT/evidence['topology_folder']; itp=folder/f'PC{dp}_GMX.itp'
    if sha(itp)!=evidence['outputs_sha256'][itp.name]:
        raise ValueError('Parameter hash changed')
    geo=json.loads((packed/'geometry-audit.json').read_text())
    if sha(packed/'packed.pdb')!=geo['packed_sha256']:
        raise ValueError('Packed coordinates changed after audit')
    BASE.mkdir(exist_ok=True)
    with (BASE/(case_name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        out=BASE/case_name; out.mkdir(exist_ok=True)
        spec={'case':case_name,'packed_sha256':geo['packed_sha256'],'itp_sha256':sha(itp),
              'protocol':'EM Fmax<1000;NVT2ps;isotropic stochastic-cell-rescaling NPT20ps',
              'temperature_K':300,'pressure_bar':1,'dt_ps':.0005,
              'compressibility_bar_inverse':4.5e-5,
              'assumptions':'300K numerical baseline,ambient-pressure assumption;compressibility is barostat tuning,not measured PC compressibility',
              'scope':'finite OH-ended oligomer dry-cell numerical/density response;not an equilibrated or hydrated substrate'}
        if (out/'spec.json').exists() and json.loads((out/'spec.json').read_text())!=spec:
            raise ValueError('Existing protocol/input spec differs')
        if (out/'result.json').exists():
            return json.loads((out/'result.json').read_text())
        write_json(out/'spec.json',spec)
        atoms=list(section_rows(itp.read_text(),'atoms'))
        molecule=list(section_rows(itp.read_text(),'moleculetype'))[0][0]
        mass=sum(float(a[7]) for a in atoms)*case['chains']
        if not (out/'packed.gro').exists():
            if (out/'em.tpr').exists():
                raise ValueError('Existing MD without preparation coordinates')
            shutil.copy2(itp,out/'pc.itp')
            rows=[r for r in (packed/'packed.pdb').read_text().splitlines() if r.startswith(('ATOM  ','HETATM'))]
            if len(rows)!=len(atoms)*case['chains']:
                raise ValueError('Coordinate/topology atom counts differ')
            text=['Periodic dry PC numerical control',str(len(rows))]
            for i,row in enumerate(rows):
                a=atoms[i%len(atoms)];xyz=[float(row[j:j+8])/10 for j in (30,38,46)]
                text.append(f'{i//len(atoms)+1:5d}{a[3]:<5.5s}{a[4]:>5.5s}{i+1:5d}'+''.join(f'{x:10.5f}' for x in xyz))
            text.append(' '.join(f'{x:.9f}' for x in case['box_nm']))
            (out/'packed.gro').write_text('\n'.join(text)+'\n')
            (out/'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n#include "pc.itp"\n'
                                         f'[ system ]\nDry PC{dp} numerical control\n[ molecules ]\n{molecule} {case["chains"]}\n')
        if sha(out/'pc.itp')!=spec['itp_sha256']:
            raise ValueError('Local topology changed')
        def run(command,label,stdin=None):
            with (out/(label+'.log')).open('a') as stream:
                subprocess.run(['gmx',*command],cwd=out,input=stdin,text=True,stdout=stream,stderr=subprocess.STDOUT,check=True)
        common='cutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\nDispCorr = EnerPres\npbc = xyz\n'
        if not (out/'em.tpr').exists():
            (out/'em.mdp').write_text(common+'integrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.002\n')
            run(['grompp','-f','em.mdp','-c','packed.gro','-p','topol.top','-o','em.tpr'],'em-prep')
        if not (out/'em.gro').exists():
            if (out/'em.log').exists():
                raise RuntimeError('Interrupted EM requires preserved-coordinate recovery;no silent reset')
            run(['mdrun','-deffnm','em','-ntmpi','1','-ntomp','1'],'em-launch')
        emlog=(out/'em.log').read_text()
        if 'converged to Fmax <' not in emlog:
            raise ValueError('EM force criterion not met')
        fmax=float(re.search(r'Maximum force\s*=\s*([\deE+.-]+)',emlog).group(1))
        stages={}
        for name,duration,previous in [('nvt',2,'em'),('npt',20,'nvt')]:
            selections='Potential\nTemperature\n0\n' if name=='nvt' else 'Potential\nTemperature\nPressure\nVolume\nDensity\n0\n'
            def energies():
                run(['energy','-f',name+'.edr','-o',name+'-energy.xvg'],name+'-energy',selections)
                data=np.array([[float(x) for x in s.split()] for s in (out/(name+'-energy.xvg')).read_text().splitlines()
                               if s and not s.startswith(('#','@'))])
                if data.ndim!=2 or data.shape[1]!=(3 if name=='nvt' else 6) or not np.isfinite(data).all():
                    raise ValueError('Missing or nonfinite energy data')
                return data
            done=False
            if (out/(name+'.edr')).exists():
                data=energies();done=abs(data[-1,0]-duration)<.001
            if not done:
                if not (out/(name+'.tpr')).exists():
                    mdp=common+f'integrator = md\ndt = 0.0005\nnsteps = {int(duration/.0005)}\nconstraints = h-bonds\n'
                    mdp+='tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\nnstenergy = 100\nnstlog = 1000\nnstxout-compressed = 200\n'
                    if name=='nvt':
                        mdp+=f'pcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = {case["seed"]}\n'
                    else:
                        mdp+='continuation = yes\ngen-vel = no\npcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 5\nref-p = 1\ncompressibility = 4.5e-5\n'
                    (out/(name+'.mdp')).write_text(mdp)
                    args=['grompp','-f',name+'.mdp','-c',previous+'.gro','-p','topol.top','-o',name+'.tpr']
                    if name=='npt':args+=['-t','nvt.cpt']
                    run(args,name+'-prep')
                args=['mdrun','-deffnm',name,'-ntmpi','1','-ntomp','1','-cpt','1','-maxh','1']
                if (out/(name+'.cpt')).exists():args+=['-cpi',name+'.cpt','-append']
                elif (out/(name+'.log')).exists():raise RuntimeError('Partial MD has no checkpoint;no reset')
                run(args,name+'-launch');data=energies()
                if abs(data[-1,0]-duration)>.001:
                    raise RuntimeError('Bounded MD incomplete;resume from its existing checkpoint')
            if 'LINCS WARNING' in (out/(name+'.log')).read_text():
                raise ValueError('Constraint warnings require review')
            if re.search(r'^WARNING',(out/(name+'-prep.log')).read_text(),re.M):
                raise ValueError('Preprocessing warnings were not accepted')
            stages[name]={'duration_ps':float(data[-1,0]),'finite_energy_frames':len(data),'temperature_mean_K':float(data[:,2].mean())}
            if name=='npt':
                early=data[data[:,0]<=5,5];late=data[data[:,0]>=15,5]
                stages[name].update(density_first_5ps_mean_kg_m3=float(early.mean()),density_last_5ps_mean_kg_m3=float(late.mean()),
                                    density_final_kg_m3=float(data[-1,5]),volume_final_nm3=float(data[-1,4]))
        result={'state':'EM/2ps NVT/20ps NPT numerical control completed','case':case_name,'PC_DP':dp,
                'atoms':case['atoms'],'chain_count':case['chains'],'topology_mass_sum_Da':mass,
                'initial_density_from_topology_kg_m3':mass/6.02214076e23/np.prod(case['box_nm'])*1e24,
                'EM_Fmax_kJ_mol_nm':fmax,'stages':stages,'physical_forcefield_approved':False,
                'equilibrium_claim':False,'limitation':'short dry finite-oligomer relaxation;no hydrated slab,coating,FXII or cPBSzI',
                'checkpoints_sha256':{s+'.cpt':sha(out/(s+'.cpt')) for s in ['nvt','npt']}}
        write_json(out/'result.json',result);print(json.dumps(result),flush=True);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('case');a=p.parse_args();execute(a.case)
