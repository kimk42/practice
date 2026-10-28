"""Recover periodic packing with a conservative generator-distance margin.

Packmol keeps precision0.01 and targets2.004A; the independent audit still uses
2.0A with the existing0.002A PDB rounding allowance. No distance gate is relaxed.
Two case workers use disjoint directories and preserve all source pose files.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import re
import shutil
import subprocess

from prepare_pc_bulk_screen import audit, write_json
from pc_length_parameters import ROOT, sha

BASE = ROOT / 'runs/pc-dry-packing-v4'
PLAN = ROOT / 'models/surface-revision-v2/dry-packing-margin-v4.json'


def prepare():
    if BASE.exists() or PLAN.exists():
        raise RuntimeError('Existing preparation preserved')
    states = subprocess.run(['ps', '-C', 'packmol', '-o', 'stat='], capture_output=True, text=True).stdout.split()
    if any(not s.startswith('Z') for s in states):
        raise RuntimeError('Inspect active Packmol before preparing recovery')
    plan = json.loads((ROOT/'models/surface-revision-v2/dry-packing-recovery-v3.json').read_text())
    BASE.mkdir()
    for case in plan['cases']:
        source = ROOT/'runs/pc-dry-packing-v3'/case['name']; dest = BASE/case['name']; dest.mkdir()
        state = json.loads((source/'status.json').read_text())['state']
        if state == 'running':
            raise RuntimeError('Unreconciled source worker')
        shutil.copy2(source/'molecule.pdb', dest/'molecule.pdb')
        if sha(source/'packmol.inp') != case['input_sha256']:
            raise RuntimeError('Source input differs from manifest')
        text = (source/'packmol.inp').read_text()
        text = re.sub(r'^precision .*$', 'precision 0.01', text, flags=re.M)
        text = re.sub(r'^tolerance .*$', 'tolerance 2.004', text, flags=re.M)
        text = re.sub(r'^restart_from .*\n', '', text, flags=re.M)
        checkpoint = source/'packing.restart'
        if not checkpoint.exists():
            checkpoint = source/'source.restart'
        if checkpoint.exists():
            if len(checkpoint.read_text().splitlines()) != case['chains']:
                raise RuntimeError('Incomplete saved poses')
            shutil.copy2(checkpoint, dest/'source.restart')
            text = 'restart_from source.restart\n'+text
            case['restart'] = {'source': str(checkpoint.relative_to(ROOT)), 'sha256': sha(checkpoint),
                               'scope': 'latest available saved rigid-body poses; unsaved optimizer state not recoverable'}
        elif state != 'prepared':
            raise RuntimeError('Previously attempted case lacks a restart; no random reset')
        else:
            case['restart'] = None
        (dest/'packmol.inp').write_text(text)
        case.update(input_sha256=sha(dest/'packmol.inp'), packmol_precision=.01,
                    generator_target_distance_A=2.004, independent_audit_distance_A=2.0)
        write_json(dest/'case.json', case); write_json(dest/'status.json', {'state': 'prepared', 'MD_started': False})
    # Obtain low-density candidates first; density changes remain explicit cases.
    plan['cases'].sort(key=lambda c:(c['starting_density_g_cm3'], c['DP'], c['seed']))
    plan.update(state='prepared;live job state must be checked',
                recovery_policy='latest available source poses retained;unstarted cases not duplicated',
                precision_change={'optimizer_precision': .01, 'generator_target_A': 2.004,
                                  'minimum_unrounded_distance_at_optimizer_limit_A': (2.004**2-.01)**.5,
                                  'independent_audit_target_A': 2.0, 'PDB_allowance_A': .002},
                worker_concurrency=2, per_case_timeout_seconds=10800)
    write_json(PLAN, plan)


def accepted_folder(case_name):
    root = BASE/case_name
    state = json.loads((root/'status.json').read_text())
    path = root/state.get('accepted_directory', '.')
    if root.resolve() not in (path.resolve(), *path.resolve().parents):
        raise ValueError('Packing output is outside its case directory')
    return path


def recover_attempt(root, state, case):
    """Retain all interrupted artifacts and resume their last saved poses."""
    source = root/state.get('attempt_directory', '.')
    if 'Success!' in (source/'packmol.log').read_text():
        previous_case=json.loads((source/'case.json').read_text())
        result=audit(source,previous_case)
        write_json(source/'geometry-audit.json',result)
        state.update(state='geometry_passed',accepted_directory=str(source.relative_to(root)))
        write_json(root/'status.json',state)
        return None,None
    checkpoint=source/'packing.restart'
    if not checkpoint.exists() or len(checkpoint.read_text().splitlines()) != case['chains']:
        raise RuntimeError('Interrupted packing has no complete saved poses; no random reset')
    out=root/f'resume-{len(list(root.glob("resume-*")))+1:04d}'
    out.mkdir()
    write_json(out/'previous-status.json',state)
    shutil.copy2(source/'molecule.pdb',out/'molecule.pdb')
    shutil.copy2(checkpoint,out/'source.restart')
    text=re.sub(r'^restart_from .*\n','',(source/'packmol.inp').read_text(),flags=re.M)
    (out/'packmol.inp').write_text('restart_from source.restart\n'+text)
    resumed=dict(case,input_sha256=sha(out/'packmol.inp'),restart={
        'source':str(checkpoint.relative_to(ROOT)), 'sha256':sha(checkpoint),
        'scope':'saved rigid-body poses; optimizer/RNG state not retained'})
    write_json(out/'case.json',resumed)
    return out,resumed


def run_case(case, queue_fd=None):
    root = BASE/case['name']
    with (root/'case.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = json.loads((root/'status.json').read_text())
        if state['state'] not in ('prepared','running'):
            return
        if state['state']=='running':
            rows=subprocess.run(['ps','-C','packmol','-o','ppid=,stat='],capture_output=True,text=True).stdout.splitlines()
            if any(r.split()[0]=='1' and not r.split()[1].startswith('Z') for r in rows):
                raise RuntimeError('Orphan Packmol process requires inspection before recovery')
            out,case=recover_attempt(root,state,case)
            if out is None:return
        else:
            out=root
            if (out/'packmol.log').exists():raise RuntimeError('Existing output preserved')
        state = {'state': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(), 'MD_started': False,
                 'attempt_directory':str(out.relative_to(root))}
        write_json(root/'status.json', state)
        try:
            with (out/'packmol.inp').open() as inp, (out/'packmol.log').open('x') as log:
                subprocess.run([shutil.which('packmol') or '/workspace/md-tools/root/usr/bin/packmol'],
                               cwd=out, stdin=inp, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=10800,
                               pass_fds=tuple(fd for fd in (lock.fileno(),queue_fd) if fd is not None))
            result = audit(out, case); write_json(out/'geometry-audit.json', result)
            state['state'] = 'geometry_passed'
            state['accepted_directory'] = str(out.relative_to(root))
        except Exception as error:
            state.update(state='blocked_preserved', reason=str(error))
        state['finished_utc'] = datetime.now(timezone.utc).isoformat(); write_json(root/'status.json', state)
        print(case['name'], state['state'], state.get('reason', ''), flush=True)


def run():
    with (BASE/'queue.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan = json.loads(PLAN.read_text())
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(lambda case:run_case(case,lock.fileno()), plan['cases']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run'])
    args = parser.parse_args(); prepare() if args.action == 'prepare' else run()
