"""Pack independent dry PC starting cells, with explicit periodic clash checks.

Packing is independent of charge assignment. No MD is run here. Cubic bulk
preparation avoids imposing a thin-film shape before a bulk density review;
the planned 6x6 nm lateral area applies only to the later slab stage.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
from rdkit import Chem
from pc_length_parameters import ROOT, preflight, source_path, sha

BASE = ROOT / 'runs/pc-dry-packing-v2'
MANIFEST = ROOT / 'models/surface-revision-v2/dry-packing-plan-v2.json'


def write_json(path, value):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def prepare():
    if BASE.exists() or MANIFEST.exists():
        raise RuntimeError('Existing packing inputs preserved')
    design = json.loads((ROOT / 'config/surface-revision-design.json').read_text())['substrate_screen']
    total = design['matched_total_BPA_units_example']
    BASE.mkdir()
    cases = []
    for dp in design['PC_BPA_units_per_chain']:
        mol, info = preflight(dp)
        if total % dp:
            raise ValueError('Noninteger number of chains')
        nchain = total // dp
        for density in design['initial_packing_density_sensitivity_g_cm3']:
            volume_nm3 = nchain * info['molecular_mass_g_mol'] / 6.02214076e23 / density * 1e21
            side_A = (volume_nm3 ** (1/3)) * 10
            for seed in design['independent_packing_seeds_planned']:
                name = f'pc{dp}-rho{density:.1f}-seed{seed}'
                out = BASE / name
                out.mkdir()
                Chem.MolToPDBFile(mol, str(out / 'molecule.pdb'))
                restart = ''
                old = ROOT / 'runs/pc-dry-packing-v1' / name / 'packing.restart'
                restart_evidence = None
                if old.exists():
                    if len(old.read_text().splitlines()) != nchain:
                        raise ValueError('Incomplete packing restart')
                    shutil.copy2(old, out / 'source.restart')
                    restart = 'restart_from source.restart\n'
                    restart_evidence = {'source': str(old.relative_to(ROOT)), 'sha256': sha(old),
                                        'reason': 'continue existing poses after removing whole-chain boundary restraint; no random reset'}
                (out / 'packmol.inp').write_text(
                    f'tolerance 2.0\nfiletype pdb\noutput packed.pdb\nseed {seed}\n'
                    f'pbc 0. 0. 0. {side_A:.9f} {side_A:.9f} {side_A:.9f}\n'
                    'nloop 200\nmaxit 200\nrestart_to packing.restart\n' + restart +
                    f'structure molecule.pdb\n number {nchain}\n'
                    ' atoms 1\n'
                    f'  inside box 0. 0. 0. {side_A:.9f} {side_A:.9f} {side_A:.9f}\n'
                    ' end atoms\n'
                    'end structure\n')
                case = {
                    'name': name, 'DP': dp, 'chains': nchain, 'BPA_units': total,
                    'atoms': nchain * info['atoms'], 'phenolic_ends': 2*nchain,
                    'carbonate_links': nchain*(dp-1),
                    'molecular_mass_g_mol': info['molecular_mass_g_mol'],
                    'cell_molecular_mass_sum_Da': nchain*info['molecular_mass_g_mol'],
                    'starting_density_g_cm3': density, 'box_nm': [side_A/10]*3,
                    'seed': seed, 'source_sha256': info['sha256'],
                    'input_sha256': sha(out / 'packmol.inp'),
                    'initial_PDB_sha256': sha(out / 'molecule.pdb'),
                    'restart': restart_evidence,
                    'density_basis': 'RDKit standard atomic weights; recompute using actual topology masses before MD',
                }
                write_json(out / 'case.json', case)
                write_json(out / 'status.json', {'state': 'prepared', 'MD_started': False})
                cases.append(case)
    report = {
        'state': 'geometry packing plan; not equilibrated substrates', 'cases': cases,
        'matched_basis': '192 BPA units per dry cubic cell; chain-end/carbonate counts and cell mass differ',
        'boundary': 'three-dimensional periodic dry cubic cells; only atom1 lies inside primary cell, allowing chains to cross faces; no slab/free interface yet',
        'later_slab_area_nm2': design['lateral_area_nm2'],
        'limitations': ['0.6 and 1.0 g/cm3 are starting hypotheses, not experimental or equilibrium densities',
                       'independent packing seeds share one intrachain conformer per length; not independent conformational ensembles',
                       'PC8/16 forcefield assignments pending; all new MD requires parameter/charge review',
                       'room-temperature short MD alone cannot establish equilibration of glassy polymer'],
    }
    write_json(MANIFEST, report)
    print(f'Prepared {len(cases)} distinct periodic packing inputs', flush=True)


def audit(out, case):
    mol, _ = preflight(case['DP'])
    if sha(source_path(case['DP'])) != case['source_sha256']:
        raise ValueError('Source molecule changed')
    if sha(out / 'packmol.inp') != case['input_sha256'] or sha(out / 'molecule.pdb') != case['initial_PDB_sha256']:
        raise ValueError('Packing input changed')
    if 'Success!' not in (out / 'packmol.log').read_text():
        raise ValueError('Packmol did not report success; preserve restart before any further attempt')
    rows = [x for x in (out / 'packed.pdb').read_text().splitlines() if x.startswith(('ATOM  ', 'HETATM'))]
    if len(rows) != case['atoms']:
        raise ValueError('Packed atom count mismatch')
    expected = [a.GetSymbol().upper() for a in mol.GetAtoms()] * case['chains']
    if [r[76:78].strip().upper() for r in rows] != expected:
        raise ValueError('Packed element order mismatch')
    xyz = np.array([[float(r[30:38]), float(r[38:46]), float(r[46:54])] for r in rows])
    box = np.array(case['box_nm']) * 10
    if not np.isfinite(xyz).all():
        raise ValueError('Nonfinite coordinates')
    n = mol.GetNumAtoms()
    source = mol.GetConformer().GetPositions()
    bond_error = 0.
    own_image_min = float('inf')
    # Reconstruct each chain across PBC before checking its nonzero images.
    for block in xyz.reshape(case['chains'], n, 3):
        whole = block.copy()
        visited, queue = {0}, [0]
        while queue:
            i = queue.pop()
            for neighbor in mol.GetAtomWithIdx(i).GetNeighbors():
                j = neighbor.GetIdx()
                if j in visited:
                    continue
                delta = block[j] - block[i]
                delta -= box * np.rint(delta / box)
                whole[j] = whole[i] + delta
                visited.add(j); queue.append(j)
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond_error = max(bond_error, abs(np.linalg.norm(whole[i]-whole[j]) - np.linalg.norm(source[i]-source[j])))
        delta = whole[:, None, :] - whole[None, :, :]
        shifts = np.rint(delta / box)
        mask = np.any(shifts != 0, axis=2)
        nearest_d2 = np.sum((delta-shifts*box)**2, axis=2)
        # If the nearest translation is zero, require a shift on one axis.
        alternate_d2 = np.min(np.sum(delta*delta, axis=2)[:, :, None] + box**2 - 2*box*np.abs(delta), axis=2)
        nonzero_d2 = np.where(mask, nearest_d2, alternate_d2)
        own_image_min = min(own_image_min, float(np.sqrt(nonzero_d2.min())))
    if bond_error > .0035:
        raise ValueError(f'Rigid bond geometry changed beyond PDB rounding: {bond_error} A')
    # Independent minimum-image check, not just trust in the optimizer label.
    owners = np.arange(len(xyz)) // n
    closest = float('inf')
    for start in range(0, len(xyz), 128):
        delta = xyz[start:start+128, None, :] - xyz[None, :, :]
        delta -= box * np.rint(delta / box)
        d2 = np.sum(delta*delta, axis=2)
        d2[owners[start:start+128, None] == owners[None, :]] = np.inf
        closest = min(closest, float(np.sqrt(d2.min())))
    # Output PDB has 0.001 A precision. At most sqrt(3)*0.001 A pair error.
    floor = 2.0 - .002
    if min(closest, own_image_min) < floor:
        raise ValueError(f'Periodic clash: intermolecular {closest}, self-image {own_image_min} A')
    return {
        'state': 'periodic packing geometry audit passed', 'case': case['name'],
        'atoms': len(xyz), 'minimum_interchain_distance_A': closest,
        'minimum_chain_to_own_periodic_image_distance_A': own_image_min,
        'max_bond_length_difference_from_source_A': bond_error,
        'tolerance_A': 2., 'PDB_serialization_allowance_A': .002,
        'packed_sha256': sha(out / 'packed.pdb'), 'MD_started': False,
        'scope': 'starting geometry only; no forcefield, density-equilibration or physical approval',
    }


def run_all():
    with (BASE / 'queue.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan = json.loads(MANIFEST.read_text())
        for case in plan['cases']:
            out = BASE / case['name']
            state = json.loads((out / 'status.json').read_text())
            if state['state'] != 'prepared':
                # Never reset or duplicate a running, completed, failed, or interrupted case.
                print(case['name'], 'preserved:', state['state'], flush=True)
                continue
            if (out / 'packmol.log').exists() or (out / 'packed.pdb').exists():
                raise RuntimeError('Existing unaccounted output requires inspection: ' + str(out))
            state = {'state': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(), 'MD_started': False}
            write_json(out / 'status.json', state)
            try:
                with (out / 'packmol.inp').open() as inp, (out / 'packmol.log').open('x') as log:
                    subprocess.run([shutil.which('packmol') or '/workspace/md-tools/root/usr/bin/packmol'],
                                   cwd=out, stdin=inp, stdout=log, stderr=subprocess.STDOUT,
                                   check=True, timeout=3600)
                result = audit(out, case)
                write_json(out / 'geometry-audit.json', result)
                state['state'] = 'geometry_passed'
            except Exception as error:
                state.update(state='blocked_preserved', reason=str(error))
            state['finished_utc'] = datetime.now(timezone.utc).isoformat()
            write_json(out / 'status.json', state)
            print(case['name'], state['state'], state.get('reason', ''), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run'])
    args = parser.parse_args()
    prepare() if args.action == 'prepare' else run_all()
