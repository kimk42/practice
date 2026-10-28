"""Continue preserved packing poses with stricter optimizer precision.

The 2.0 A separation and independent output audit are unchanged. Revision2 used
Packmol's default precision=0.01 (squared-distance violation), which can accept
about1.9975 A before PDB rounding. Tighten optimizer precision, not the audit.
This script never starts or restarts quantum calculations or MD.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

import prepare_pc_bulk_screen as packing
from pc_length_parameters import ROOT, sha

SOURCE = ROOT / 'runs/pc-dry-packing-v2'
OUT = ROOT / 'runs/pc-dry-packing-v3'
PLAN = ROOT / 'models/surface-revision-v2/dry-packing-recovery-v3.json'


def prepare():
    if subprocess.run(['pgrep', '-x', 'packmol'], capture_output=True).returncode == 0:
        raise RuntimeError('Inspect active Packmol process before recovery')
    if OUT.exists() or PLAN.exists():
        raise RuntimeError('Existing recovery preserved; inspect it before any new attempt')
    plan = json.loads((ROOT / 'models/surface-revision-v2/dry-packing-plan-v2.json').read_text())
    OUT.mkdir()
    for case in plan['cases']:
        source, dest = SOURCE / case['name'], OUT / case['name']
        state = json.loads((source / 'status.json').read_text())['state']
        if state == 'running':
            raise RuntimeError('Unreconciled running status in ' + str(source))
        if sha(source / 'packmol.inp') != case['input_sha256'] or sha(source / 'molecule.pdb') != case['initial_PDB_sha256']:
            raise RuntimeError('Revision2 source input changed')
        restart = source / 'packing.restart'
        if state != 'prepared' and not restart.exists():
            raise RuntimeError('Previously attempted case has no recoverable poses; no random reset')
        dest.mkdir()
        shutil.copy2(source / 'molecule.pdb', dest / 'molecule.pdb')
        text = (source / 'packmol.inp').read_text()
        text = re.sub(r'^restart_from .*\n', '', text, flags=re.M)
        text = 'precision 0.0001\n' + text
        case['recovery_source_state'] = state
        case['recovery_source_input_sha256'] = case['input_sha256']
        if restart.exists():
            if len(restart.read_text().splitlines()) != case['chains']:
                raise RuntimeError('Incomplete rigid-body restart')
            shutil.copy2(restart, dest / 'source.restart')
            text = 'restart_from source.restart\n' + text
            # Near-converged low-density cases benefit from continuous local
            # refinement instead of moving bad molecules to new trial positions.
            if case['starting_density_g_cm3'] == .6:
                text = 'disable_movebad\n' + text
                text = text.replace('maxit 200\n', 'maxit 1000\n')
            case['restart'] = {'source': str(restart.relative_to(ROOT)), 'sha256': sha(restart),
                               'reason': 'continue saved rigid-body poses; no random reset'}
        else:
            case['restart'] = None
        (dest / 'packmol.inp').write_text(text)
        case['input_sha256'] = sha(dest / 'packmol.inp')
        case['packmol_precision'] = .0001
        packing.write_json(dest / 'case.json', case)
        packing.write_json(dest / 'status.json', {'state': 'prepared', 'MD_started': False})
    plan.update(state='recovery inputs prepared; check live worker status',
                previous_plan='models/surface-revision-v2/dry-packing-plan-v2.json',
                recovery_policy='4 saved-pose cases continued;8 previously unstarted cases prepared;no QM restart',
                precision_change={'old_default': .01, 'new': .0001,
                                  'quantity': 'maximum squared-distance/restraint violation;not the 2.0A target distance'},
                independent_audit_unchanged=True)
    packing.write_json(PLAN, plan)
    print('Prepared12 recovery cases:4 from stored poses,8 not previously started', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run'])
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    else:
        packing.BASE, packing.MANIFEST = OUT, PLAN
        packing.run_all()
