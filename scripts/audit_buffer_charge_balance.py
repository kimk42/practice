"""Enumerate formal-charge alternatives, without selecting pH populations or charges."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT / 'config/cpbszi-system-preparation.json'
    manifest = json.loads(source.read_text())
    assert manifest['ready_for_grompp'] is False
    alternatives = []
    for case in manifest['cases']:
        chloride, iodide, citrate, azide = [
            item['nearest_formula_units'] for item in case['fixed_analytical_counts']
        ]
        for hypothesis in case['phosphate_alternatives']:
            phosphate = hypothesis['nearest_total_molecules']
            for macrocharge in (0, 6):
                sodium_counts = []
                for cit3 in range(citrate + 1):
                    for hpo4 in range(phosphate + 1):
                        anion_charge = (chloride + iodide + azide + 3 * cit3
                                        + 2 * (citrate - cit3) + 2 * hpo4
                                        + phosphate - hpo4)
                        sodium = anion_charge - macrocharge
                        assert sodium >= 0
                        assert sodium + macrocharge - anion_charge == 0
                        sodium_counts.append(sodium)
                alternatives.append({
                    'reference_case': case['reference_case'],
                    'phosphate_hydration_hypothesis': hypothesis['hydration_hypothesis'],
                    'phosphate_total': phosphate,
                    'macromolecule_formal_charge_e': macrocharge,
                    'enumerated_population_pairs': len(sodium_counts),
                    'neutralizing_total_Na_count_range': [min(sodium_counts), max(sodium_counts)],
                    'all_formal_charge_residuals_e': 0,
                })
    result = {
        'state': 'formal integer bookkeeping only; no mixture selected',
        'source_manifest_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'restricted_species_basis': ['Cit3-', 'HCit2-', 'HPO4²-', 'H2PO4-', 'N3-'],
        'scope': 'all integer partitions within this restricted candidate basis; not exhaustive acid-base chemistry or predicted pH7.4 speciation',
        'macrocharge_basis': 'formal Q=0/+6 sensitivity at each reference water inventory;these are not two prepared coating-only/FXII systems',
        'sodium_policy': 'conditional total for global neutrality; stock-counterion versus buffer adjustment policy remains unselected',
        'alternatives': alternatives,
        'partial_charges_used': False,
        'MD_topology_exported': False,
        'ready_for_grompp': False,
    }
    target = ROOT / 'models/buffer-parameters/formal-charge-balance-audit.json'
    if target.exists():
        legacy = target.with_name('formal-charge-balance-legacy-whole-box.json')
        previous = json.loads(target.read_text())
        if 'reference_case' not in previous['alternatives'][0] and not legacy.exists():
            legacy.write_bytes(target.read_bytes())
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
