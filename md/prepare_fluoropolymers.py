"""Generate finite fluoropolymer candidates; no crystal/surface validation implied.

One common PVDF charge set is used across conformers so a nominal phase change
cannot silently change chemical parameters. Torsion motifs are starting-chain
constructions, not experimental unit cells or validated polymorphs.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolDescriptors,rdMolTransforms
from openmm import unit,NonbondedForce,Platform,Context,VerletIntegrator
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,NoCutoff
from .project_paths import amber_home
from workbench.manager import atomic

MOTIFS={'pvdf-alpha':[180,60,180,-60],'pvdf-beta':[180],'pvdf-gamma':[180,180,180,60,180,180,180,-60],'ptfe':[180]}


def molecule(family,repeats=12):
 if family=='PVDF':smiles='CC(F)(F)'*repeats
 elif family=='PTFE':smiles='FC(F)(F)'+'C(F)(F)'*(2*repeats-2)+'C(F)(F)F'
 else:raise ValueError('Family must be PTFE or PVDF')
 mol=Chem.AddHs(Chem.MolFromSmiles(smiles));p=AllChem.ETKDGv3();p.randomSeed=42
 if AllChem.EmbedMolecule(mol,p)!=0:raise ValueError('Molecular embedding failed')
 return mol


def conformer(mol,name):
 result=Chem.Mol(mol);carbon=[a.GetIdx() for a in result.GetAtoms() if a.GetSymbol()=='C'];motif=MOTIFS[name]
 if any(result.GetBondBetweenAtoms(a,b) is None for a,b in zip(carbon,carbon[1:])):raise ValueError('Carbon atom order is not a continuous backbone')
 for i in range(len(carbon)-3):rdMolTransforms.SetDihedralDeg(result.GetConformer(),*carbon[i:i+4],motif[i%len(motif)])
 result.SetProp('_Name',name+'_FINITE_UNVALIDATED');return result


def save_molecule(mol,path):
 with Chem.SDWriter(str(path)) as w:w.write(mol)


def command(argv,cwd,env,timeout=1800):
 log=Path(cwd)/(Path(argv[0]).name+'.log')
 with log.open('w') as f:
  p=subprocess.Popen(argv,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:code=p.wait(timeout=timeout)
  except subprocess.TimeoutExpired:
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
   raise TimeoutError(f'{Path(argv[0]).name} exceeded {timeout}s; logs retained')
 if code:raise RuntimeError(f'{Path(argv[0]).name} exited {code}; inspect {log}')
 return {'argv':argv,'returncode':code,'log':log.name}


def normalize_charge(path):
 raw=Path(path).read_text();Path(path).with_suffix('.raw.mol2').write_text(raw);lines=raw.splitlines();a=lines.index('@<TRIPOS>ATOM')+1;b=next(i for i in range(a,len(lines)) if lines[i].startswith('@<TRIPOS>'));indices=[i for i in range(a,b) if lines[i].strip()];total=sum(float(lines[i].split()[-1]) for i in indices)
 if abs(total)>.01:raise ValueError('Charge residual exceeds rounding-normalization limit')
 for i in indices:
  fields=lines[i].split();fields[-1]=f'{float(fields[-1])-total/len(indices):.10f}';lines[i]=' '.join(fields)
 Path(path).write_text('\n'.join(lines)+'\n');return {'raw_charge_sum_e':total,'uniform_rounding_correction_per_atom_e':-total/len(indices)}


def parameterize(family,out,amber):
 out=Path(out);out.mkdir(parents=True,exist_ok=False);amber=Path(amber).resolve();env={**os.environ,'AMBERHOME':str(amber),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'};mol=molecule(family);names=['ptfe'] if family=='PTFE' else ['pvdf-alpha','pvdf-beta','pvdf-gamma'];source=conformer(mol,'ptfe' if family=='PTFE' else 'pvdf-beta');save_molecule(source,out/'charge-reference.sdf')
 for name in names:save_molecule(conformer(mol,name),out/(name+'.sdf'))
 status={'family':family,'status':'parameterizing','atoms':mol.GetNumAtoms(),'formula':rdMolDescriptors.CalcMolFormula(mol),'crystal_validated':False,'surface_validated':False,'adsorption_MD_allowed':False};atomic(out/'status.json',status)
 commands=[];start=time.monotonic()
 try:
  commands.append(command([str(amber/'bin/antechamber'),'-i','charge-reference.sdf','-fi','sdf','-o','polymer.mol2','-fo','mol2','-at','gaff2','-c','bcc','-nc','0','-rn','PTF' if family=='PTFE' else 'PVF','-s','2'],out,env));correction=normalize_charge(out/'polymer.mol2')
  commands.append(command([str(amber/'bin/parmchk2'),'-i','polymer.mol2','-f','mol2','-o','polymer.frcmod','-s','gaff2'],out,env))
  if 'ATTN' in (out/'polymer.frcmod').read_text():raise ValueError('Unresolved/estimated parameters flagged ATTN; requires review')
  (out/'isolated.leap').write_text('source leaprc.gaff2\nloadamberparams polymer.frcmod\nPOL = loadmol2 polymer.mol2\nsaveamberparm POL polymer.prmtop polymer.inpcrd\nsavepdb POL polymer.pdb\nquit\n')
  commands.append(command([str(amber/'bin/tleap'),'-f','isolated.leap'],out,env))
  prm=AmberPrmtopFile(str(out/'polymer.prmtop'));inp=AmberInpcrdFile(str(out/'polymer.inpcrd'));system=prm.createSystem(nonbondedMethod=NoCutoff,constraints=None);nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce));q=np.array([nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())]);symbols=[a.element.symbol for a in prm.topology.atoms()]
  if symbols!=[a.GetSymbol() for a in mol.GetAtoms()] or abs(q.sum())>1e-5:raise ValueError('Amber atom mapping/neutrality mismatch')
  integrator=VerletIntegrator(.5*unit.femtosecond);context=Context(system,integrator,Platform.getPlatformByName('Reference'));checks=[]
  for name in names:
   candidate=conformer(mol,name);xyz=np.array(candidate.GetConformer().GetPositions())*.1;context.setPositions(xyz*unit.nanometer);s=context.getState(getEnergy=True,getForces=True);energy=s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole);forces=s.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer)
   if not np.isfinite(energy) or not np.isfinite(forces).all():raise ValueError('Nonfinite isolated-candidate force or energy')
   check={'candidate':name,'formula':rdMolDescriptors.CalcMolFormula(candidate),'initial_torsion_motif_deg':MOTIFS[name],'atoms':len(xyz),'net_charge_e':float(q.sum()),'finite_potential_kJ_mol':float(energy),'max_force_kJ_mol_nm':float(np.linalg.norm(forces,axis=1).max()),'molecular_point_charge_dipole_D':(np.sum(q[:,None]*(xyz-xyz.mean(axis=0)),axis=0)*48.0320427).tolist(),'shared_charge_reference':'charge-reference.sdf','status':'finite_single_molecule_force_evaluation_only','crystal_validated':False,'surface_validated':False,'MD_steps':0};checks.append(check)
  status.update(status='technical_parameter_checks_passed_physical_validation_pending',charge_method='GAFF2 / AM1-BCC candidate only',charge_normalization=correction,checks=checks,commands=commands,elapsed_wall_seconds=time.monotonic()-start,limitations=['Finite oligomer, not a crystalline surface or validated polymer force field.','PVDF motif names identify starting torsions only; no lattice, packing, crystallographic face or domain has been validated.','PTFE trans starting candidate does not represent its bulk helical chain.','Shared PVDF atom mapping and charge set across all motifs.','Finite static energy/force and molecular point-charge dipole are technical diagnostics, not adsorption strength or bulk polarization.'])
  status['sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!='status.json'}
 except Exception as e:status.update(status='needs_review',error=str(e),elapsed_wall_seconds=time.monotonic()-start)
 atomic(out/'status.json',status);return status


def main():
 p=argparse.ArgumentParser();p.add_argument('--family',required=True,choices=['PTFE','PVDF']);p.add_argument('--output',required=True);p.add_argument('--amberhome',default=str(amber_home()));a=p.parse_args();print(json.dumps(parameterize(a.family,a.output,a.amberhome),indent=2),flush=True)

if __name__=='__main__':main()
