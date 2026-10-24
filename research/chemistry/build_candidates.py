"""Generate chemistry candidates, NOT simulation-ready force-field topologies."""
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
out=Path(__file__).resolve().parent
branch='C(=O)NCCC[N+](C)(C)CC(=O)[O-]'
specs={
 'cbaa1_monomer': 'C=CC(=O)NCCC[N+](C)(C)CC(=O)[O-]',
 'cbaa1_capped_trimer_example': 'C[C@H]('+branch+')C[C@@H]('+branch+')C[C@H]('+branch+')C',
 'bma_n_butyl_methacrylate': 'C=C(C)C(=O)OCCCC',
 'bpaa_acrylamide_literature_candidate': 'C=CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1',
 'paper_small_molecule_betaine': 'C[N+](C)(C)CC(=O)[O-]',
}
report={}
for name,smiles in specs.items():
 mol=Chem.MolFromSmiles(smiles)
 if mol is None: raise ValueError(name)
 assert len(Chem.GetMolFrags(mol))==1
 assert Chem.GetFormalCharge(mol)==0
 positives=sum(a.GetFormalCharge()==1 for a in mol.GetAtoms())
 negatives=sum(a.GetFormalCharge()==-1 for a in mol.GetAtoms())
 expected=3 if 'trimer' in name else (0 if name.startswith(('bma_', 'bpaa_')) else 1)
 assert positives==negatives==expected
 if name=='cbaa1_monomer':assert rdMolDescriptors.CalcMolFormula(mol)=='C10H18N2O3'
 mol=Chem.AddHs(mol)
 params=AllChem.ETKDGv3();params.randomSeed=2026
 if AllChem.EmbedMolecule(mol,params)!=0:raise RuntimeError('Embedding failed: '+name)
 # Generic MMFF relaxation is only a geometry preparation step.
 if not AllChem.MMFFHasAllMoleculeParams(mol):raise RuntimeError('MMFF unsupported: '+name)
 status=AllChem.MMFFOptimizeMolecule(mol,maxIters=2000)
 if status!=0:raise RuntimeError('MMFF did not converge: '+name)
 mol.SetProp('_Name',name+' | chemistry candidate only')
 mol.SetProp('SOURCE', '10.1021/acs.langmuir.8b02540 SI (BMA)' if name.startswith('bma_') else ('10.1080/14686996.2022.2119883 Table 1, citing Lin 2020 (BPAA candidate)' if name.startswith('bpaa_') else '10.1039/d3sc01977b Fig. 1 / Scheme 1; reconstructed connectivity'))
 mol.SetProp('STATUS','No production charges, atom types, graft anchor, or validated torsions')
 with Chem.SDWriter(str(out/(name+'.sdf'))) as w:w.write(mol)
 # Normalize the writer timestamp so regeneration does not dirty versioned inputs.
 sdf_path=out/(name+'.sdf')
 lines=sdf_path.read_text().splitlines()
 lines[1]='     RDKit          3D'
 sdf_path.write_text('\n'.join(lines)+'\n')
 Chem.MolToPDBFile(mol,str(out/(name+'.pdb')))
 report[name]={'input_smiles':smiles,'canonical_smiles':Chem.MolToSmiles(Chem.RemoveHs(mol)),
 'formula':rdMolDescriptors.CalcMolFormula(mol),'molar_mass_g_mol':Descriptors.MolWt(mol),
 'formal_charge':Chem.GetFormalCharge(mol),'positive_centers':positives,'negative_centers':negatives,
 'atom_count_with_hydrogen':mol.GetNumAtoms(),'geometry_MMFF_converged':status==0}
(out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
