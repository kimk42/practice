"""Read-only directional derivative diagnostics at the preserved PC8 geometry.

Every displaced calculation is isolated from the optimization checkpoints.
No diagnostic coordinate, charge or gradient is approved for parameter export.
"""
from datetime import datetime, timezone
import argparse
import fcntl
import json
from pathlib import Path
import subprocess

import numpy as np

from pc_length_parameters import ROOT, sha
from pc_sqm_lbfgs import parse_singlepoint, validate_coordinates
from run_pc_sqm_checkpointed import amber_environment, read_initial, write_json

OUT = ROOT / 'runs/pc8-gradient-diagnostic-20261030'


def main(tight_density=False):
    global OUT
    if tight_density:
        OUT = ROOT / 'runs/pc8-tight-p-gradient-diagnostic-20261030'
    source = ROOT / 'runs/pc-length-sqm-lbfgs-v3/pc8'
    with (ROOT/'runs/pc-length-sqm-checkpointed-v2/pc8.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        saved = json.loads((source/'accepted-checkpoint.json').read_text())
        raw = source/saved['evaluation']/'sqm.out'
        if sha(raw) != saved['output_sha256']:
            raise ValueError('Accepted evidence changed')
        elements, names, _ = read_initial(ROOT/'runs/pc-length-sqm-checkpointed-v2/pc8/original-sqm.in')
        xyz = np.array(saved['coordinates_A'])
        _, original_gradient, _ = parse_singlepoint(raw, len(elements))
        direction = original_gradient / np.linalg.norm(original_gradient)
        OUT.mkdir(exist_ok=True)
        spec = {'source_checkpoint_sha256': sha(source/'accepted-checkpoint.json'),
                'source_evaluation': saved['evaluation'], 'source_output_sha256': sha(raw),
                'direction': 'unit Euclidean norm of original analytic gradient',
                'displacement_norms_A': ([0., .00001, -.00001, .00003, -.00003, .0001, -.0001]
                                         if tight_density else [0., .0001, -.0001, .001, -.001, .01, -.01]),
                'scf_energy_thresholds': [1e-10] if tight_density else [1e-10, 1e-12],
                'engine_sha256': sha(Path(amber_environment()['AMBERHOME'])/'bin/sqm'),
                'scope': 'diagnostic single points only; original optimizer state unchanged'}
        if tight_density:
            spec['density_threshold'] = 1e-10
        if (OUT/'spec.json').exists() and json.loads((OUT/'spec.json').read_text()) != spec:
            raise ValueError('Diagnostic source changed')
        write_json(OUT/'spec.json', spec)
        results = {}
        env = amber_environment()
        for exponent in ([10] if tight_density else [10, 12]):
            records = {}
            for h in spec['displacement_norms_A']:
                step_label = f'{h:+.5f}' if tight_density else f'{h:+.4f}'
                case = OUT/f'scf{exponent}-step{step_label}'
                trial = xyz + h*direction
                if (case/'record.json').exists():
                    result = json.loads((case/'record.json').read_text())
                    if sha(case/'sqm.out') != result['output_sha256']:
                        raise ValueError('Diagnostic output changed')
                else:
                    case.mkdir()  # Interrupted/incomplete directories require review.
                    text = ('PC8 directional derivative diagnostic\n &qmmm\n'
                            f" qm_theory='AM1', qmcharge=0, scfconv=1.d-{exponent}, ndiis_attempts=700,\n"
                            ' maxcyc=0, verbosity=4, diag_routine=6, pseudo_diag=0,\n'
                            + (' tight_p_conv=1,\n' if tight_density else '') + ' /\n')
                    text += ''.join(f'{z:3d} {name:5s} '+ ' '.join(f'{v:.17g}' for v in row)+'\n'
                                    for z, name, row in zip(elements, names, trial))
                    (case/'sqm.in').write_text(text)
                    np.save(case/'input-coordinates.npy', trial)
                    with (case/'launch.log').open('x') as log:
                        subprocess.run([env['AMBERHOME']+'/bin/sqm', '-O', '-i', 'sqm.in', '-o', 'sqm.out'],
                                       cwd=case, env=env, stdout=log, stderr=subprocess.STDOUT,
                                       check=True, timeout=300, pass_fds=(lock.fileno(),))
                    energy, gradient, info = parse_singlepoint(case/'sqm.out', len(elements))
                    validate_coordinates(case/'sqm.out', trial, elements)
                    result = {'energy_kcal_mol': energy, 'step_norm_A': h,
                              'projected_gradient_kcal_mol_A': float(np.sum(gradient*direction)),
                              'gradient_kcal_mol_A': gradient.tolist(), **info}
                    write_json(case/'record.json', result)
                records[h] = result
                print(case.name, result['energy_kcal_mol'], result['rms_gradient_kcal_mol_A'], flush=True)
            rows = []
            for h in ([.00001, .00003, .0001] if tight_density else [.0001, .001, .01]):
                finite = (records[h]['energy_kcal_mol'] - records[-h]['energy_kcal_mol'])/(2*h)
                analytic = records[0.]['projected_gradient_kcal_mol_A']
                rows.append({'step_norm_A': h, 'central_difference': finite, 'analytic': analytic,
                             'absolute_error': abs(finite-analytic),
                             'descent_energy_change': records[-h]['energy_kcal_mol']-records[0.]['energy_kcal_mol']})
            results[str(exponent)] = rows
            write_json(OUT/'partial-analysis.json', results)
        write_json(OUT/'result.json', {'checked_utc': datetime.now(timezone.utc).isoformat(),
                                     'finite_difference': results, 'spec': spec,
                                     'physical_approval': False})
        print(json.dumps(results), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tight-density', action='store_true')
    main(parser.parse_args().tight_density)
