"""Formal-chemistry inventory only; no invented partial charges or speciation fractions."""
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'models/buffer-parameters';out.mkdir(exist_ok=True)
    candidates=[('citrate3-','O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]',-3),('hydrogen_citrate2-','O=C(O)CC(O)(CC(=O)[O-])C(=O)[O-]',-2),('dihydrogen_phosphate-','O=P([O-])(O)O',-1),('hydrogen_phosphate2-','O=P([O-])([O-])O',-2),('azide-','[N-]=[N+]=[N-]',-1)]
    rows=[]
    for name,smiles,charge in candidates:
        mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
        if Chem.GetFormalCharge(mol)!=charge:raise ValueError('Formal chemistry mismatch')
        mol.SetProp('_Name',name);mol.SetProp('STATUS','formal graph candidate only;no assigned force-field partial charges or production coordinates')
        Chem.MolToMolFile(mol,str(out/(name+'.mol')))
        rows.append({'name':name,'smiles':smiles,'formal_charge':charge,'formula':rdMolDescriptors.CalcMolFormula(mol),'atoms':mol.GetNumAtoms(),'partial_charges_assigned':False,'topology_assigned':False,'population_fraction':None})
    inventory={'candidates':rows,'hydrogen_citrate_caution':'single formal protonation microstate candidate,not a validated native microstate/population','speciation_status':'fractions atpH7.4 require sourced pKa/ionic-condition model or documented sensitivity;not selected here','finite_box':{'volume_nm3':3888,'analytical_molecules_per_mM':3888*6.02214076e-4,'interpretation':'rounding/sodium charge balance must follow accepted speciation;no final integer mixture assigned'},'buffer_ready_for_MD':False}
    (out/'polyatomic-species-audit.json').write_text(json.dumps(inventory,indent=2)+'\n')
if __name__=='__main__':main()
