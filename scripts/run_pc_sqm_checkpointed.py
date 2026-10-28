"""Restartable native-SQM geometry segments with an explicit gradient gate.

Each finished segment retains SQM's final coordinates (printed at0.0001A).
Continuation restarts the optimizer/SCF from those coordinates; it is not a
bitwise optimizer-state checkpoint. No charge is accepted from a maxcyc stop
unless the reported gradient, including its rounding bound, meets the target.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
from acpype.utils import bundled_amber_dir
from pc_length_parameters import ROOT, preflight, sha

BASE = ROOT / 'runs/pc-length-sqm-checkpointed-v2'
TARGET_GRMS = .0005


def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n'); temp.replace(path)


def amber_environment():
    home = bundled_amber_dir()
    env = dict(os.environ, AMBERHOME=home, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    env['PATH'] = home + '/bin:' + env['PATH']
    env['LD_LIBRARY_PATH'] = '/workspace/md-tools/root/usr/lib/x86_64-linux-gnu/openblas-pthread' + (
        ':' + env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    return env


def read_initial(path):
    text = path.read_text()
    rows = text.split('\n /\n')[1].splitlines()
    atoms = [line.split() for line in rows if line.strip()]
    return [int(a[0]) for a in atoms], [a[1] for a in atoms], np.array([[float(x) for x in a[2:5]] for a in atoms])


def parse_completed(path, elements):
    text = path.read_text()
    if 'Calculation Completed' not in text or 'Final Structure' not in text:
        raise ValueError('No complete coordinate-bearing SQM output')
    if re.search(r'SCF.*(?:failed|not converg)|XMIN ERROR|FATAL', text, re.I):
        raise ValueError('SQM reports a convergence/runtime error')
    final = text.rsplit('Final Structure', 1)[1]
    rows = re.findall(r'^\s*QMMM:\s+(\d+)\s+(\d+)\s+([A-Z][a-z]?)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', final, re.M)
    symbols = {1: 'H', 6: 'C', 8: 'O'}
    if len(rows) != len(elements) or [int(r[0]) for r in rows] != list(range(1, len(elements)+1)):
        raise ValueError('Incomplete final coordinate table')
    if [r[2] for r in rows] != [symbols[e] for e in elements]:
        raise ValueError('Final atom order/elements changed')
    xyz = np.array([[float(x) for x in row[3:6]] for row in rows])
    if not np.isfinite(xyz).all():
        raise ValueError('Nonfinite final coordinates')
    records = re.findall(r'^xmin\s+(\d+)\s+([-\d.]+) kcal/mol\s+([-\d.]+) kcal/\(mol\*A\)', text, re.M)
    if not records:
        raise ValueError('Missing numerical gradient; completion marker alone is insufficient')
    iteration, energy, gradient = records[-1]
    gradient = float(gradient)
    if not np.isfinite(gradient) or gradient < 0:
        raise ValueError('Invalid gradient')
    return xyz, {'iterations': int(iteration), 'heat_of_formation_kcal_mol_printed': float(energy),
                 'rms_gradient_kcal_mol_A_printed': gradient,
                 'gradient_print_rounding_upper_bound': gradient+.00005,
                 'gradient_gate_passed': gradient+.00005 <= TARGET_GRMS,
                 'coordinate_precision_A': .0001, 'output_sha256': sha(path)}


def input_text(elements, names, xyz, cycles):
    header = ("Checkpointed PC AM1 geometry segment\n &qmmm\n"
              " qm_theory='AM1', qmcharge=0, scfconv=1.d-10, ndiis_attempts=700,\n"
              f" grms_tol=0.0004, maxcyc={cycles}, ntpr=1,\n /\n")
    return header + ''.join(f'{z:3d} {name:5s} {r[0]:15.8f} {r[1]:15.8f} {r[2]:15.8f}\n'
                            for z, name, r in zip(elements, names, xyz))


def run(dp):
    mol, evidence = preflight(dp)
    BASE.mkdir(exist_ok=True)
    with (BASE / f'pc{dp}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out = BASE / f'pc{dp}'; out.mkdir(exist_ok=True)
        initial = ROOT / ('runs/parameterization-long/pc4/PC4.acpype/sqm.in' if dp == 4 else
                          f'runs/pc-length-parameters-v1/pc{dp}/PC{dp}.acpype/sqm.in')
        elements, names, xyz = read_initial(initial)
        coordinate_source = initial
        if dp == 4:
            coordinate_source = initial.with_name('sqm.out')
            xyz, _ = parse_completed(coordinate_source, elements)
        if elements != [a.GetAtomicNum() for a in mol.GetAtoms()]:
            raise ValueError('Original SQM input atom ordering differs from PC graph')
        spec = {'input': evidence, 'initial_sqm_input_sha256': sha(initial),
                'sqm_binary_sha256': sha(Path(bundled_amber_dir())/'bin/sqm'),
                'AM1_SCF_threshold': 1e-10, 'geometry_target_GRMS': TARGET_GRMS,
                'requested_solver_GRMS': .0004, 'coordinate_restart_precision_A': .0001,
                'restart_method': 'native SQM segments;coordinate continuation with optimizer/SCF state reset',
                'initial_rerun_authorized': 'user requested run through completion after interruption was explained'}
        if dp == 4:
            spec.update(initial_coordinate_source=str(coordinate_source.relative_to(ROOT)),
                        initial_coordinate_source_sha256=sha(coordinate_source))
        if (out/'spec.json').exists():
            if json.loads((out/'spec.json').read_text()) != spec:
                raise ValueError('Existing run spec differs; do not overwrite it')
        else:
            write_json(out/'spec.json', spec); shutil.copy2(initial, out/'original-sqm.in')
        completed, total_iterations = [], 0
        # Recover completed segments even if the parent died before writing status.
        for folder in sorted(out.glob('segment-*-attempt-*')):
            if not (folder/'sqm.out').exists():
                continue
            text = (folder/'sqm.out').read_text()
            if 'Calculation Completed' not in text:
                continue
            candidate, info = parse_completed(folder/'sqm.out', elements)
            lengths = [float(np.linalg.norm(candidate[b.GetBeginAtomIdx()]-candidate[b.GetEndAtomIdx()])) for b in mol.GetBonds()]
            if min(lengths) < .7 or max(lengths) > 2.0:
                raise ValueError('Recovered coordinate checkpoint fails bond-geometry gate')
            segment = int(folder.name.split('-')[1])
            if segment != len(completed)+1:
                raise ValueError('Nonconsecutive or duplicate completed segments require inspection')
            xyz = candidate; total_iterations += info['iterations']; completed.append(folder)
            info['bond_lengths_A'] = [min(lengths), max(lengths)]
            write_json(folder/'segment-audit.json', info)
            if info['gradient_gate_passed']:
                write_json(out/'completed.json', {'state': 'AM1_gradient_gate_passed', 'final_output': str((folder/'sqm.out').relative_to(ROOT)),
                                                'total_logged_iterations': total_iterations, **info})
                write_json(out/'status.json', {'state': 'AM1_gradient_gate_passed', 'final_segment': folder.name})
                return
        started = time.monotonic()
        env = amber_environment()
        while total_iterations < 4000:
            if shutil.disk_usage(out).free < 2*1024**3:
                raise RuntimeError('Less than2GiB free; outputs preserved')
            segment = len(completed)+1
            attempt = len(list(out.glob(f'segment-{segment:04d}-attempt-*')))+1
            folder = out/f'segment-{segment:04d}-attempt-{attempt:02d}'; folder.mkdir()
            last_gradient = info['rms_gradient_kcal_mol_A_printed'] if completed else float('inf')
            # Near a minimum, retain the native optimizer memory for longer.
            # Repeated very short segments can otherwise spend most work
            # relaxing the0.0001A coordinate serialization perturbation.
            cycles = 300 if dp == 4 else (150 if dp == 8 and last_gradient < .05 else
                                         30 if dp == 16 and last_gradient < .05 else
                                         50 if dp == 8 else 10)
            cycles = min(cycles, 4000-total_iterations)
            (folder/'sqm.in').write_text(input_text(elements, names, xyz, cycles))
            state = {'state': 'running', 'worker_pid': os.getpid(), 'segment': segment, 'attempt': attempt,
                     'completed_segments': len(completed), 'total_logged_iterations': total_iterations,
                     'input_sha256': sha(folder/'sqm.in'), 'started_utc': datetime.now(timezone.utc).isoformat()}
            write_json(out/'status.json', state)
            try:
                with (folder/'launch.log').open('x') as log:
                    subprocess.run([env['AMBERHOME']+'/bin/sqm', '-O', '-i', 'sqm.in', '-o', 'sqm.out'],
                                   cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=7200)
                xyz, info = parse_completed(folder/'sqm.out', elements)
                lengths = [float(np.linalg.norm(xyz[b.GetBeginAtomIdx()]-xyz[b.GetEndAtomIdx()])) for b in mol.GetBonds()]
                if min(lengths) < .7 or max(lengths) > 2.0:
                    raise ValueError('QM geometry changed chemical bond lengths outside the preflight envelope')
                info['bond_lengths_A'] = [min(lengths), max(lengths)]
                write_json(folder/'segment-audit.json', info)
                total_iterations += info['iterations']; completed.append(folder)
                write_json(out/'coordinate-checkpoint.json', {'segment': segment, 'atoms': elements, 'names': names,
                           'coordinates_A': xyz.tolist(), 'total_logged_iterations': total_iterations, **info})
                print(f'PC{dp} segment{segment} RMS={info["rms_gradient_kcal_mol_A_printed"]} iterations={total_iterations}', flush=True)
                if info['gradient_gate_passed']:
                    write_json(out/'completed.json', {'state': 'AM1_gradient_gate_passed', 'final_output': str((folder/'sqm.out').relative_to(ROOT)),
                               'total_logged_iterations': total_iterations, **info})
                    write_json(out/'status.json', {'state': 'AM1_gradient_gate_passed', 'final_segment': folder.name})
                    return
                if info['iterations'] < 1:
                    raise ValueError('No optimization progress')
            except Exception as error:
                state.update(state='blocked_preserved', reason=str(error)); write_json(out/'status.json', state); raise
            if time.monotonic()-started > 86400:
                write_json(out/'status.json', {'state': 'coordinate_checkpoint_saved_at_24h_bound', 'completed_segments': len(completed)})
                return
        write_json(out/'status.json', {'state': 'blocked_preserved', 'reason': '4000 iterations without gradient convergence'})
        raise RuntimeError('AM1 geometry did not converge within4000 logged iterations')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dp', type=int, choices=[4, 8, 16], required=True)
    args = parser.parse_args(); run(args.dp)
