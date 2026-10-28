"""Read the completed control; check geometry/energy without modifying trajectories."""
import hashlib
import json
import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'runs/azide-rigid-water-control'


def table(name):
    return np.array([[float(x) for x in line.split()] for line in
                     (RUN / name).read_text().splitlines()
                     if line.strip() and not line.startswith(('#', '@'))])


def main():
    energy = table('energy.xvg')
    # gmx traj XVG uses six significant digits; use nine-decimal G96 exported
    # directly from the same TRR rather than loosening geometry tolerances.
    lines = (RUN / 'sites.g96').read_text().splitlines()
    frames, time = [], None
    for index, line in enumerate(lines):
        if line == 'TIMESTEP':
            time = float(lines[index + 1].split()[1])
        if line == 'POSITIONRED':
            positions = [float(v) for row in lines[index + 1:index + 6] for v in row.split()]
            assert lines[index + 6] == 'END' and len(positions) == 15
            frames.append([time, *positions])
    sites = np.array(frames)
    assert energy.shape == (41, 3) and sites.shape == (41, 16)
    assert np.isfinite(energy).all() and np.isfinite(sites).all()
    assert np.allclose(energy[:, 0], sites[:, 0], atol=1e-7)
    assert energy[-1, 0] == 2 and energy[0, 0] == 0
    representation = json.loads((RUN / 'representation.json').read_text())
    candidate = json.loads((ROOT / 'models/buffer-parameters/published-azide-candidate.json').read_text())
    for stem in ['em', 'md']:
        log = (RUN / (stem + '.log')).read_text()
        assert 'Finished mdrun' in log
        assert not re.search(r'LINCS WARNING|constraint error|\bnan\b', log, re.I)
        assert not re.search(r'^WARNING', (RUN / (stem + '-grompp.stdout.log')).read_text(), re.M)
    em = (RUN / 'em.log').read_text()
    fmax = float(re.search(r'Maximum force\s*=\s*([\deE+.-]+)', em).group(1))
    assert fmax < 1000
    xyz = sites[:, 1:].reshape((-1, 5, 3))
    # Reimage every site relative to the first carrier in the cubic3nm box.
    delta = xyz - xyz[:, :1]
    delta -= 3 * np.rint(delta / 3)
    xyz = xyz[:, :1] + delta
    axis = xyz[:, 1] - xyz[:, 0]
    expected_sites = (xyz[:, :1] + np.array(representation['linear_virtual_site_coefficients'])
                      [None, :, None] * axis[:, None, :])
    vsite_error = float(np.max(np.abs(xyz[:, 2:] - expected_sites)))
    center_error = float(np.max(np.linalg.norm(xyz[:, 3] - xyz[:, :2].mean(axis=1), axis=1)))
    left, right = xyz[:, 2] - xyz[:, 3], xyz[:, 4] - xyz[:, 3]
    lengths = np.stack([np.linalg.norm(left, axis=1), np.linalg.norm(right, axis=1)], axis=1)
    length_error = float(np.max(np.abs(lengths - candidate['equal_NN_distance_nm'])))
    angles = np.degrees(np.arccos(np.clip((left * right).sum(axis=1)
                        / np.linalg.norm(left, axis=1) / np.linalg.norm(right, axis=1), -1, 1)))
    carrier_error = float(np.max(np.abs(np.linalg.norm(axis, axis=1)
                        - 2 * representation['carrier_positions_from_COM_nm'][1])))
    # Tolerance includes mixed-precision integration and nine-decimal export.
    assert max(vsite_error, center_error, length_error, carrier_error) < 5e-6
    assert np.max(np.abs(angles - 180)) < 0.01
    hashes = {name: hashlib.sha256((RUN / name).read_bytes()).hexdigest() for name in
              ['md.tpr', 'md.trr', 'md.edr', 'md.cpt', 'md.gro', 'md.mdp',
               'em.log', 'md.log', 'azide.itp', 'azide-types.itp', 'representation.json', 'sites.g96']}
    result = {
        'state': '2ps rigid-azide/TIP3P/JC-sodium numerical control completed',
        'run_directory': str(RUN.relative_to(ROOT)),
        'water_molecules': int(re.search(r'^SOL\s+(\d+)', (RUN / 'topol.top').read_text(), re.M).group(1)),
        'azide_ions': 1, 'sodium_ions': 1,
        'total_coordinate_sites_including_mass_carriers': int((RUN / 'md.gro').read_text().splitlines()[1]),
        'duration_ps': float(energy[-1, 0]), 'saved_frames_checked': len(sites),
        'EM_Fmax_kJ_mol_nm': fmax,
        'temperature_mean_K': float(energy[:, 2].mean()),
        'temperature_endpoint_K': float(energy[-1, 2]),
        'potential_range_kJ_mol': [float(energy[:, 1].min()), float(energy[:, 1].max())],
        'max_NN_length_error_nm': length_error,
        'NNN_angle_range_deg': [float(angles.min()), float(angles.max())],
        'max_carrier_constraint_error_nm': carrier_error,
        'max_virtual_site_construction_error_nm': vsite_error,
        'max_COM_vs_central_site_error_nm': center_error,
        'constraint_warning_detected': False,
        'input_and_output_sha256': hashes,
        'interpretation': 'finite energy and saved-frame rigid geometry passed;does not validate hydration,ion pairing,concentration,polarization or protein interaction accuracy',
        'physical_force_field_approved': False,
        'full_cpbszi_ready': False,
    }
    target = ROOT / 'models/buffer-parameters/azide-rigid-water-control-audit.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    (RUN / 'status.json').write_text(json.dumps({'state': result['state'], 'audit': str(target.relative_to(ROOT))}, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
