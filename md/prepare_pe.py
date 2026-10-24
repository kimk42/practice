"""Prepare a neutral C24H50 finite PE oligomer with GAFF2/AM1-BCC using AmberTools.
The oligomer and later finite slab are exploratory models, not a validated bulk PE surface.
"""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolTransforms

try:
    from .project_paths import amber_home, resolve_project_path
except ImportError:
    from project_paths import amber_home, resolve_project_path


def prepare(output,amberhome):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    amber=Path(amberhome).resolve();env={**os.environ,'AMBERHOME':str(amber),'OMP_NUM_THREADS':'2'}
    mol=Chem.AddHs(Chem.MolFromSmiles('C'*24))
    params=AllChem.ETKDGv3();params.randomSeed=42
    if AllChem.EmbedMolecule(mol,params)!=0:raise ValueError('Conformer generation failed')
    for i in range(21):rdMolTransforms.SetDihedralDeg(mol.GetConformer(),i,i+1,i+2,i+3,180)
    converged=AllChem.UFFOptimizeMolecule(mol,maxIters=2000)
    if converged!=0:raise ValueError('Initial oligomer conformer optimization did not converge')
    mol.SetProp('_Name','PE24_GAFF2_AM1BCC')
    with Chem.SDWriter(str(out/'oligomer.sdf')) as writer:writer.write(mol)
    commands=[
        [str(amber/'bin/antechamber'),'-i','oligomer.sdf','-fi','sdf','-o','PE.mol2','-fo','mol2','-at','gaff2','-c','bcc','-nc','0','-rn','PEO','-s','2'],
        [str(amber/'bin/parmchk2'),'-i','PE.mol2','-f','mol2','-o','PE.frcmod','-s','gaff2'],
    ]
    for i,cmd in enumerate(commands):
        print('Running',Path(cmd[0]).name,flush=True)
        with (out/f'command-{i}.log').open('w') as log:subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    # SQM charge output is rounded; preserve raw output and distribute only a small residual
    # uniformly to enforce the molecule's known net charge, as in standard charge normalization.
    raw=(out/'PE.mol2').read_text();(out/'PE.raw.mol2').write_text(raw)
    lines=raw.splitlines();start=lines.index('@<TRIPOS>ATOM')+1
    stop=next(i for i in range(start,len(lines)) if lines[i].startswith('@<TRIPOS>'))
    atom_lines=[i for i in range(start,stop) if lines[i].strip()]
    raw_charge=sum(float(lines[i].split()[-1]) for i in atom_lines)
    if abs(raw_charge)>.01:raise ValueError('Charge residual too large for rounding normalization')
    correction=-raw_charge/len(atom_lines)
    for i in atom_lines:
        fields=lines[i].split();fields[-1]=f'{float(fields[-1])+correction:.10f}';lines[i]=' '.join(fields)
    (out/'PE.mol2').write_text('\n'.join(lines)+'\n')
    if 'ATTN' in (out/'PE.frcmod').read_text():raise ValueError('Unresolved GAFF2 parameters in frcmod')
    leap='source leaprc.gaff2\nloadamberparams PE.frcmod\nPEO = loadmol2 PE.mol2\nsaveamberparm PEO polymer.prmtop polymer.inpcrd\nsavepdb PEO polymer.pdb\nquit\n'
    (out/'isolated.leap').write_text(leap)
    with (out/'tleap.log').open('w') as log:subprocess.run([str(amber/'bin/tleap'),'-f','isolated.leap'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    from openmm.app import AmberPrmtopFile,AmberInpcrdFile,NoCutoff
    from openmm import NonbondedForce,unit
    p=AmberPrmtopFile(str(out/'polymer.prmtop'));x=AmberInpcrdFile(str(out/'polymer.inpcrd'))
    system=p.createSystem(nonbondedMethod=NoCutoff,constraints=None)
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
    charge=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
    if abs(charge)>1e-5 or system.getNumParticles()!=74:raise ValueError('Oligomer charge/atom count differs from neutral C24H50')
    manifest={'model':'finite hydrogen-terminated PE oligomer','formula':'C24H50','repeat_units':12,'atoms':74,'forcefield':'GAFF2','charge_method':'AM1-BCC (AmberTools antechamber / sqm)',
              'total_charge_e':charge,'raw_AM1BCC_charge_sum_e':raw_charge,'uniform_charge_normalization_per_atom_e':correction,'amberhome':str(amber),'sha256':{name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['oligomer.sdf','PE.raw.mol2','PE.mol2','PE.frcmod','polymer.prmtop','polymer.inpcrd']},
              'limitations':['Finite PE oligomer; terminal methyl groups differ from infinite bulk PE','No claim of validated PE crystallinity, density, surface energy or protein affinity','UFF used only for initial conformer; MD parameters are GAFF2/AM1-BCC']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output');p.add_argument('--amberhome',default=str(amber_home()));a=p.parse_args();prepare(a.output,a.amberhome)
