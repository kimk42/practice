"""Analytical-count alternatives only; no invented pH populations or final MD mixture."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    source=ROOT/'config/experimental-buffer.json';recipe=json.loads(source.read_text());volume=18*12*18;factor=volume*6.02214076e-4
    fixed=[];phosphate=None
    for item in recipe['components']:
        if 'conditional_total_mM' in item:
            phosphate=[{'hydration_hypothesis':key,'target_total_mM':v,'nearest_total_molecules':round(v*factor),'rounded_total_mM':round(v*factor)/factor} for key,v in item['conditional_total_mM'].items()]
        else:
            value=item.get('nominal_mM',item.get('nominal_total_mM'));number=round(value*factor)
            fixed.append({'component':item['name'],'target_analytical_mM':value,'nearest_formula_units':number,'rounded_analytical_mM':number/factor})
    manifest={'state':'analytical preparation alternatives;not a runnable matched-buffer system','source_recipe_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'box_nm':[18,12,18],'volume_nm3':volume,'formula_units_per_mM':factor,'fixed_analytical_counts':fixed,'phosphate_alternatives':phosphate,'unassigned_populations':['N_Cit3+N_HCit2=N_citrate_total','N_HPO4+N_H2PO4=N_phosphate_total'],'charge_balance_constraint':'N_Na+Q_macromolecule=N_Cl+N_I+3*N_Cit3+2*N_HCit2+2*N_HPO4+N_H2PO4+N_azide','counterion_policy':'not selected;document whether protein-stock counterions are added versus adjusting buffer Na,then use matched policy across bare/coated systems','free_NaOH_added_as_species':False,'finite_box_protein_density':{'protein_count':1,'nominal_full_box_molarity_M':1/(6.02214076e23*volume*1e-24),'experiment_reference_target_nM':3.75,'scope':'full-box count density only,not accessible-water-phase concentration and not a reproduction of the experimental dilute assay'},'ready_for_grompp':False,'blocking_validation':['azide automatic bond-typing warning;no approved azide partial charges/topology','citrate/phosphate bond-order representation and charge/speciation review','species populations atpH7.4 not assigned','phosphate hydration alternatives not resolved as a single experimental value']}
    (ROOT/'config/cpbszi-system-preparation.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
