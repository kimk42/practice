"""Advance the authorized PC preparation through strict dry numerical controls.

Resume coordinate-checkpointed QM workers when inactive; never duplicate a live
worker. Only audited parameters and packed cases advance. Scientific failures
are recorded and held for review, not converted into success or bypassed.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from pc_length_parameters import ROOT
from run_pc_sqm_checkpointed import BASE as QM_BASE, write_json
from pc_qm_route import route

STATUS=ROOT/'runs/pc-completion-status.json'


def packing_ready_or_launch():
    base=ROOT/'runs/pc-dry-packing-v4'
    with (base/'queue.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return 'running'
        ready={'prepared','running','reviewed_local_polish','prepared_anchor_variant','reviewed_relocation'}
        states=[json.loads(p.read_text())['state'] for p in base.glob('*/status.json')]
        if not any(state in ready for state in states):return 'idle_completed_or_held'
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    with (ROOT/'runs'/f'pc-packing-auto-{stamp}.log').open('x') as log:
        subprocess.Popen([sys.executable,'-u',str(ROOT/'scripts/run_pc_packing_margin.py'),'run'],
                         cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    return 'launched_ready_or_explicitly_reviewed_cases'


def qm_ready_or_launch(dp):
    active_base, worker = route(dp)
    folder=active_base/f'pc{dp}'
    if (folder/'completed.json').exists():
        return 'gradient_passed'
    QM_BASE.mkdir(exist_ok=True)
    with (QM_BASE/f'pc{dp}.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return 'running'
        state=json.loads((folder/'status.json').read_text()) if (folder/'status.json').exists() else {}
        if state.get('state')=='blocked_preserved':return 'blocked_preserved'
        # An orphaned QM process may no longer be covered by its parent lock.
        process_rows=subprocess.run(['ps','-C','sqm','-o','ppid=,stat='],capture_output=True,text=True).stdout.splitlines()
        if any(row.split()[0]=='1' and not row.split()[1].startswith('Z') for row in process_rows):
            return 'orphan_QM_requires_inspection'
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    log=ROOT/f'runs/pc{dp}-auto-continuation-{stamp}.log'
    with log.open('x') as stream:
        subprocess.Popen([sys.executable,'-u',str(ROOT/'scripts'/worker),'--dp',str(dp)],
                         cwd=ROOT,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
    return 'launched_from_saved_coordinates_or_authorized_initial_input'


def main(watch):
    with (ROOT/'runs/pc-completion.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        failures={}
        if STATUS.exists():failures=json.loads(STATUS.read_text()).get('held_failures',{})
        started=time.monotonic()
        cases=json.loads((ROOT/'models/surface-revision-v2/dry-packing-margin-v4.json').read_text())['cases']
        while True:
            packing_queue=packing_ready_or_launch()
            quantum={str(dp):qm_ready_or_launch(dp) for dp in (4,8,16)}
            for dp in (4,8,16):
                key=f'parameters-PC{dp}'
                if quantum[str(dp)]=='gradient_passed' and key not in failures:
                    if not (ROOT/f'runs/pc-length-forcefield-v2/pc{dp}/parameter-audit.json').exists():
                        with (ROOT/f'runs/pc{dp}-automatic-finalization.log').open('a') as log:
                            result=subprocess.run([sys.executable,str(ROOT/'scripts/finalize_pc_checkpointed.py'),'--dp',str(dp)],
                                                  cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                        if result.returncode:failures[key]='Postprocessing failed;inspect preserved log/directory before retry'
            completed=[];packing={}
            for case in cases:
                name=case['name'];key='bulk-'+name
                packed=ROOT/'runs/pc-dry-packing-v4'/name
                packing[name]=json.loads((packed/'status.json').read_text())['state']
                result_path=ROOT/'runs/pc-bulk-controls-v2'/name/'result.json'
                if result_path.exists():completed.append(name);continue
                if key in failures or packing[name]!='geometry_passed':continue
                if not (ROOT/f'runs/pc-length-forcefield-v2/pc{case["DP"]}/parameter-audit.json').exists():continue
                running={'state':'running_numerical_control','active_case':name,'quantum':quantum,
                         'completed_controls':completed,'held_failures':failures,'updated_utc':datetime.now(timezone.utc).isoformat()}
                write_json(STATUS,running)
                with (ROOT/'runs'/('bulk-'+name+'-automatic.log')).open('a') as log:
                    result=subprocess.run([sys.executable,str(ROOT/'scripts/run_pc_bulk_control.py'),name],
                                          cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                if result.returncode:failures[key]='Numerical control failed;inspect checkpoint/log before any continuation'
                elif result_path.exists():completed.append(name)
            state={'state':'completed_numerical_screen' if len(completed)==len(cases) else 'awaiting_audited_prerequisites',
                   'quantum':quantum,'packing':packing,'completed_controls':completed,'target_controls':len(cases),
                   'held_failures':failures,'updated_utc':datetime.now(timezone.utc).isoformat(),
                   'packing_queue':packing_queue,
                   'scope':'PC4/8/16 dry EM/2ps NVT/20ps NPT screen only;not full coating/FXII/cPBSzI or physical approval'}
            write_json(STATUS,state)
            if len(completed)==len(cases):
                reports=[json.loads((ROOT/'runs/pc-bulk-controls-v2'/n/'result.json').read_text()) for n in completed]
                write_json(ROOT/'models/surface-revision-v2/completed-dry-screen.json',{'cases':reports,'physical_approval':False})
                return
            if not watch or time.monotonic()-started>86400:
                return
            time.sleep(30)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--watch',action='store_true');a=p.parse_args();main(a.watch)
