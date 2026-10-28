"""Water-referenced planning alternatives;no final mixture or pH population assigned."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
NA=6.02214076e23
WATER_REFERENCE_M=55.5  # explicit approximate solvent reference,not a measurement


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def counts(recipe,factor):
    fixed=[];phosphate=[]
    for item in recipe['components']:
        if 'conditional_total_mM' in item:
            phosphate=[{'hydration_hypothesis':key,'target_total_mM':value,
                        'nearest_total_molecules':round(value*factor),
                        'rounded_water_reference_mM':round(value*factor)/factor}
                       for key,value in item['conditional_total_mM'].items()]
        else:
            value=item.get('nominal_mM',item.get('nominal_total_mM'));number=round(value*factor)
            fixed.append({'component':item['name'],'target_analytical_mM':value,
                          'nearest_formula_units':number,'rounded_water_reference_mM':number/factor})
    return fixed,phosphate


def main():
    source=ROOT/'config/experimental-buffer.json';recipe=json.loads(source.read_text())
    target=ROOT/'config/cpbszi-system-preparation.json'
    if target.exists() and 'schema_version' not in json.loads(target.read_text()):
        old=ROOT/'models/buffer-parameters/cpbszi-whole-box-planning-legacy.json'
        if not old.exists():old.write_bytes(target.read_bytes())
    cases=[]
    for name in ['fxii-coating-near-10ps','fxii-bare-near-10ps']:
        path=ROOT/'runs'/name/'topol.top';lines=path.read_text().split('[ molecules ]',1)[1].splitlines()
        composition={}
        for line in lines:
            row=line.split(';',1)[0].split()
            if len(row)==2:composition[row[0]]=composition.get(row[0],0)+int(row[1])
        water=composition['SOL'];factor=water/(WATER_REFERENCE_M*1000)
        fixed,phosphate=counts(recipe,factor)
        cases.append({'reference_case':name,'source_topology_sha256':sha(path),
                      'reference_water_count':water,'water_reference_molarity_M':WATER_REFERENCE_M,
                      'water_equivalent_volume_nm3':water/(WATER_REFERENCE_M*NA)*1e24,
                      'formula_units_per_mM':factor,'macromolecule_charge_e':6,
                      'legacy_ions':{k:v for k,v in composition.items() if k in ['NA','CL']},
                      'legacy_ion_policy':'do not add a new recipe on top of these ions;define final ion inventory separately',
                      'fixed_analytical_counts':fixed,'phosphate_alternatives':phosphate})
    result={'schema_version':2,'state':'provisional water-referenced count alternatives;not a runnable matched-buffer system',
            'source_recipe_sha256':sha(source),'cases':cases,
            'count_basis':'N_component=target_mM/(1000*55.5M)*N_water;55.5M is an explicit reference assumption',
            'interpretation':'water-equivalent volume is not measured accessible or bulk-solution volume;accounts approximately for solvent inventory only',
            'final_preparation_requirements':['recompute after final solvation/replacement because water count changes',
                  'measure a defined bulk-water region after relaxation and report local ion partitioning separately',
                  'state consistent counterion and protonation policies for coated/bare systems'],
            'unassigned_populations':['N_Cit3+N_HCit2=N_citrate_total','N_HPO4+N_H2PO4=N_phosphate_total'],
            'charge_balance_constraint':'N_Na+Q_macro=N_Cl+N_I+3*N_Cit3+2*N_HCit2+2*N_HPO4+N_H2PO4+N_azide',
            'counterion_policy':'not selected;explicitly choose stock-counterions or buffer-Na adjustment',
            'free_NaOH_added_as_species':False,'ready_for_grompp':False,
            'blocking_validation':['citrate/phosphate representation and aqueous parameter review',
                  'choose sourced pH population model or explicit sensitivity scenarios',
                  'document hydrate sensitivity scenario;single measured hydrate not mandatory for exploratory comparison',
                  'assess azide and multianion aqueous transfer;source access and rigid azide implementation already resolved'],
            'legacy_whole_box_reference':'models/buffer-parameters/cpbszi-whole-box-planning-legacy.json'}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps([{'case':x['reference_case'],'waters':x['reference_water_count'],
                       'water_equivalent_nm3':x['water_equivalent_volume_nm3'],
                       'fixed_counts':[r['nearest_formula_units'] for r in x['fixed_analytical_counts']]}
                      for x in cases],indent=2))


if __name__=='__main__':main()
