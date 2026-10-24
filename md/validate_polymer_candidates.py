"""Independent finite-candidate graph/charge/overlap checks, not phase validation."""
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from openmm import NonbondedForce,unit
from openmm.app import AmberPrmtopFile,NoCutoff
from workbench.manager import atomic


def validate(directory):
 directory=Path(directory);status=json.loads((directory/'status.json').read_text())
 if status['status']!='technical_parameter_checks_passed_physical_validation_pending':return {'status':'waiting_or_parameter_review','parameter_status':status['status'],'surface_validated':False}
 reference=next(m for m in Chem.SDMolSupplier(str(directory/'charge-reference.sdf'),removeHs=False) if m is not None);prm=AmberPrmtopFile(str(directory/'polymer.prmtop'));system=prm.createSystem(nonbondedMethod=NoCutoff,constraints=None);nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce));atoms=list(prm.topology.atoms());expected={tuple(sorted([b.GetBeginAtomIdx(),b.GetEndAtomIdx()])) for b in reference.GetBonds()};actual={tuple(sorted([a.index,b.index])) for a,b in prm.topology.bonds()}
 if expected!=actual or [a.GetSymbol() for a in reference.GetAtoms()]!=[a.element.symbol for a in atoms]:raise ValueError('Amber graph/atom mapping differs from charge-reference SDF')
 q=sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(len(atoms)))
 if abs(q)>1e-5:raise ValueError('Candidate is not neutral within rounding tolerance')
 excluded={tuple(sorted([int(nb.getExceptionParameters(i)[0]),int(nb.getExceptionParameters(i)[1])])) for i in range(nb.getNumExceptions())};pairs=np.array([(i,j) for i in range(len(atoms)) for j in range(i+1,len(atoms)) if (i,j) not in excluded]);heavy=np.array([a.element.symbol!='H' for a in atoms]);checks=[]
 for file in sorted(directory.glob('*.sdf')):
  if file.name=='charge-reference.sdf':continue
  mol=next(m for m in Chem.SDMolSupplier(str(file),removeHs=False) if m is not None)
  if expected!={tuple(sorted([b.GetBeginAtomIdx(),b.GetEndAtomIdx()])) for b in mol.GetBonds()}:raise ValueError('Conformer graph differs')
  if [a.GetSymbol() for a in mol.GetAtoms()]!=[a.element.symbol for a in atoms]:raise ValueError('Conformer element order differs')
  xyz=mol.GetConformer().GetPositions()*.1;d=np.linalg.norm(xyz[pairs[:,0]]-xyz[pairs[:,1]],axis=1);mask=heavy[pairs[:,0]]&heavy[pairs[:,1]];hard=int(np.count_nonzero((d<.08)|(mask&(d<.12))));checks.append({'candidate':file.stem,'atom_graph_matches_Amber':True,'minimum_nonexception_allatom_distance_nm':float(d.min()),'hard_overlap_pairs':hard,'initial_geometry_hard_overlap_screen_passed':hard==0})
 result={'status':'independent_candidate_mapping_charge_checks_passed','atoms':len(atoms),'element_counts':dict(Counter(a.element.symbol for a in atoms)),'charge_e':q,'Amber_graph_matches_reference':True,'checks':checks,'crystal_validated':False,'surface_validated':False,'production_MD_allowed':False,'interpretation':'Checks molecular identity, common charge mapping and hard initial overlaps only; neither crystal structure nor polymer/water-interface accuracy is established.'};atomic(directory/'independent-validation.json',result);return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('directory');a=p.parse_args();print(json.dumps(validate(a.directory),indent=2))
