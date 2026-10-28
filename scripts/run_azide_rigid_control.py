"""Separate numerical control of a fixed-charge rigid azide candidate, not validation."""
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/azide-rigid-water-control'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if OUT.exists():
        raise RuntimeError('Existing control preserved; inspect its logs/checkpoint instead of restarting')
    candidate_path = ROOT / 'models/buffer-parameters/published-azide-candidate.json'
    candidate = json.loads(candidate_path.read_text())
    r = candidate['equal_NN_distance_nm']
    # Two noninteracting mass carriers at +/-d reproduce mass, COM and both
    # transverse principal moments of three equal-mass nuclei at -r,0,+r.
    nitrogen_mass = 14.0067  # ordinary natural-abundance N mass convention, dalton
    carrier_mass = 1.5 * nitrogen_mass
    d = r * math.sqrt(2 / 3)
    coefficients = [(1 - r / d) / 2, 0.5, (1 + r / d) / 2]
    assert math.isclose(2 * carrier_mass, 3 * nitrogen_mass)
    assert math.isclose(2 * carrier_mass * d**2, 2 * nitrogen_mass * r**2)
    assert math.isclose(sum(candidate['partial_charges_e']), -1, abs_tol=1e-12)
    OUT.mkdir()
    (OUT / 'status.json').write_text(json.dumps({'state': 'preparing separate numerical control'}) + '\n')
    sigma = candidate['all_N_sites_LJ']['sigma_nm']
    epsilon = candidate['all_N_sites_LJ']['epsilon_kJ_mol_from_1998']
    itp = ('; Exploratory fixed-charge rigid ion; not a full reproduction of1998 polarization\n'
           '[ atomtypes ]\nAZMC 0 0 0 A 0 0\n'
           f'AZN 7 {nitrogen_mass} 0 A {sigma:.12f} {epsilon:.12f}\n'
           '[ moleculetype ]\nAZR 0\n[ atoms ]\n')
    for i in (1, 2):
        itp += f'{i} AZMC 1 AZR MC{i} {i} 0 {carrier_mass:.12f}\n'
    for i, (name, charge) in enumerate(zip(['N1', 'NC', 'N2'], candidate['partial_charges_e']), 3):
        itp += f'{i} AZN 1 AZR {name} {i} {charge:.12f} 0\n'
    itp += f'[ constraints ]\n1 2 1 {2*d:.15f}\n[ virtual_sites2 ]\n'
    for i, a in enumerate(coefficients, 3):
        itp += f'{i} 1 2 1 {a:.15f}\n'
    itp += '[ exclusions ]\n' + ''.join(
        f'{i} ' + ' '.join(str(j) for j in range(i + 1, 6)) + '\n' for i in range(1, 5))
    types, molecule = itp.split('[ moleculetype ]', 1)
    (OUT / 'azide-types.itp').write_text(types)
    (OUT / 'azide.itp').write_text('[ moleculetype ]' + molecule)
    provenance = {
        'state': 'numerical-control candidate only',
        'published_candidate_sha256': sha(candidate_path),
        'mass_convention_N_Da': nitrogen_mass,
        'carrier_mass_Da': carrier_mass,
        'carrier_positions_from_COM_nm': [-d, d],
        'linear_virtual_site_coefficients': coefficients,
        'virtual_site_rule': 'x_s=(1-a)*x_1+a*x_2;GROMACS2025.2 local reference manual,interaction-methods eq284',
        'total_mass_Da': 2 * carrier_mass,
        'physical_transverse_inertia_Da_nm2': 2 * nitrogen_mass * r**2,
        'carrier_transverse_inertia_Da_nm2': 2 * carrier_mass * d**2,
        'longitudinal_inertia_Da_nm2': 0,
        'COM_offset_nm': 0,
        'intramolecular_interactions': 'all ten atom pairs explicitly excluded;only carrier-carrier distance constrained',
        'no_harmonic_azide_terms': True,
        'limitations': ['fixed published reference charges;not polarizable1998 model',
                       'TIP3P/JC sodium transfer unvalidated;not a reproduction of SPC or DMSO studies',
                       '2ps300K NVT is only a numerical check,not experimental temperature/calibration or hydration validation'],
        'azide_itp_sha256': sha(OUT / 'azide.itp'),
        'approved_for_matched_buffer_MD': False,
    }
    (OUT / 'representation.json').write_text(json.dumps(provenance, indent=2) + '\n')
    rows = []
    for i, (name, x) in enumerate(zip(['MC1', 'MC2', 'N1', 'NC', 'N2'],
                                    [1-d, 1+d, 1-r, 1, 1+r]), 1):
        rows.append(f'{1:5d}{"AZR":<5}{name:>5}{i:5d}{x:13.8f}{1.:13.8f}{1.:13.8f}')
    rows.append(f'{2:5d}{"NAJC":<5}{"JC_NA":>5}{6:5d}{2.:13.8f}{2.:13.8f}{2.:13.8f}')
    (OUT / 'ions.gro').write_text('Rigid azide candidate and JC sodium\n6\n' + '\n'.join(rows) + '\n3 3 3\n')
    ions_path = ROOT / 'models/buffer-parameters/jc-tip3p-monovalent.itp'
    (OUT / 'topol.top').write_text('#include "amber99sb-ildn.ff/forcefield.itp"\n'
        f'#include "azide-types.itp"\n#include "{ions_path}"\n#include "azide.itp"\n'
        '#include "amber99sb-ildn.ff/tip3p.itp"\n'
        '[ system ]\nRigid azide numerical TIP3P control\n[ molecules ]\nAZR 1\nNAJC 1\n')

    def run(args, name, stdin=None):
        with (OUT / (name + '.stdout.log')).open('w') as log:
            subprocess.run(['gmx', *args], cwd=OUT, input=stdin, text=True,
                           stdout=log, stderr=subprocess.STDOUT, check=True)

    common = ('cutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\n'
              'rvdw = 1.0\npbc = xyz\nconstraints = none\nlincs-order = 6\n')
    run(['solvate', '-cp', 'ions.gro', '-cs', 'spc216.gro', '-radius', '0.17',
         '-p', 'topol.top', '-o', 'solvated.gro'], 'solvate')
    (OUT / 'em.mdp').write_text(common + 'integrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.002\n')
    run(['grompp', '-f', 'em.mdp', '-c', 'solvated.gro', '-p', 'topol.top',
         '-o', 'em.tpr', '-pp', 'em-processed.top'], 'em-grompp')
    run(['mdrun', '-deffnm', 'em', '-ntmpi', '1', '-ntomp', '2'], 'em-run')
    if 'converged to Fmax <' not in (OUT / 'em.log').read_text():
        raise RuntimeError('EM target not met; no MD started')
    mdp = (common + 'integrator = md\ndt = 0.0005\nnsteps = 4000\n'
           'tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\n'
           'pcoupl = no\ngen-vel = yes\ngen-temp = 300\ngen-seed = 2026\n'
           'nstenergy = 100\nnstlog = 1000\nnstxout = 100\nnstvout = 1000\n'
           'nstxout-compressed = 100\ncontinuation = no\n')
    (OUT / 'md.mdp').write_text(mdp)
    run(['grompp', '-f', 'md.mdp', '-c', 'em.gro', '-p', 'topol.top',
         '-o', 'md.tpr', '-pp', 'md-processed.top'], 'md-grompp')
    (OUT / 'status.json').write_text(json.dumps({'state': '2ps numerical MD running'}) + '\n')
    run(['mdrun', '-deffnm', 'md', '-ntmpi', '1', '-ntomp', '2', '-cpt', '1'], 'md-run')
    run(['energy', '-f', 'md.edr', '-o', 'energy.xvg'], 'energy', 'Potential\nTemperature\n0\n')
    (OUT / 'sites.ndx').write_text('[ AZIDE ]\n1 2 3 4 5\n')
    run(['traj', '-f', 'md.trr', '-s', 'md.tpr', '-n', 'sites.ndx', '-ox', 'sites.xvg'],
        'sites', 'AZIDE\n')
    run(['trjconv', '-f', 'md.trr', '-s', 'md.tpr', '-n', 'sites.ndx', '-o',
         'sites.g96', '-ndec', '9'], 'sites-g96', 'AZIDE\n')
    (OUT / 'status.json').write_text(json.dumps({'state': 'MD finished;independent numerical analysis pending'}) + '\n')
    print(OUT)


if __name__ == '__main__':
    main()
