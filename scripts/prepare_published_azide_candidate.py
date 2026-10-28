"""Transcribe published azide values; never export an approved MD topology."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    inventory = ROOT / 'research/proteins/modeling-evidence/azide-uploaded-sources.json'
    for item in json.loads(inventory.read_text())['sources']:
        assert hashlib.sha256((ROOT / item['private_local_pdf']).read_bytes()).hexdigest() == item['pdf_sha256']
    sigma = 3.32 * 0.1
    epsilon = 0.08485 * 4.184
    charges = [-0.9215, 0.8430, -0.9215]
    assert math.isclose(sum(charges), -1, abs_tol=1e-12)
    # Independent reproduction of the LJ functional-form landmarks.
    def lj(r):
        return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
    assert abs(lj(sigma)) < 1e-12
    assert math.isclose(lj(2**(1/6) * sigma), -epsilon, abs_tol=1e-12)
    candidate = {
        'state': 'published numerical candidate;not approved for matched-buffer MD',
        'source_inventory_sha256': hashlib.sha256(inventory.read_bytes()).hexdigest(),
        'atomic_order': ['terminal N', 'central N', 'terminal N'],
        'partial_charges_e': charges,
        'charge_sum_e': sum(charges),
        'charge_source': 'Morita1998 p5515 sectionIII.C (uploaded PDF page6);independently repeated Olschewski2011 p214504-7 (PDF page8)',
        'equal_NN_distance_nm': 1.1626 * 0.1,
        'linear_angle_deg': 180,
        'geometry_source': 'Morita1998 p5514 sectionIII.A;Olschewski2011 p214504-7',
        'all_N_sites_LJ': {
            'sigma_nm': sigma,
            'epsilon_kJ_mol_from_1998': epsilon,
            'epsilon_kJ_mol_2011_rounded': 0.355,
            'source': 'Morita1998 p5518 sectionIV.A (PDF page9),citing Ferrario;Olschewski2011 p214504-7 corroborates',
            'sigma_interpretation': 'explicit sigma in1998 resolves ambiguous word radius in2011;not Rmin or half sigma',
            'cross_mixing': 'arithmetic sigma,geometric epsilon;Morita1998 equation29',
            'LJ_at_sigma_kJ_mol': lj(sigma),
            'LJ_at_minimum_kJ_mol': lj(2**(1/6) * sigma),
        },
        'source_model_limits': [
            '1998 uses flexible SPC solvent,rigid solute distance/angle SHAKE,and field-dependent solute charges;fixed charges alone do not reproduce that full model',
            '2011 uses SPC water and Ewald;does not provide a transferable table of harmonic azide bond/angle constants in the inspected MD-method passage',
            '2025 uses OPLS-AA in DMSO;SI explicitly says neither bonds nor angles constrained,without an explicit azide bonded-constant table in section2',
        ],
        'remaining_gates': [
            'choose and verify a nonsingular GROMACS representation of rigid linear azide,or obtain a complete supported flexible model',
            'water/counterion interaction validation for Amber/TIP3P and Joung-Cheatham sodium;Lorentz-Berthelot form matching alone is not physical validation',
            'numerical import/energy/geometry checks before any coupled-buffer use',
        ],
        'original_1993_pdf_required_to_read_these_repeated_values': False,
        'harmonic_bond_or_angle_constants_invented': False,
        'MD_topology_exported': False,
        'approved_for_MD': False,
    }
    target = ROOT / 'models/buffer-parameters/published-azide-candidate.json'
    target.write_text(json.dumps(candidate, indent=2) + '\n')
    print(json.dumps({'charges_e': charges, 'sum_e': sum(charges), 'sigma_nm': sigma,
                      'epsilon_kJ_mol': epsilon, 'numerical_conversion_checks': 'passed',
                      'approved_for_MD': False}, indent=2))


if __name__ == '__main__':
    main()
