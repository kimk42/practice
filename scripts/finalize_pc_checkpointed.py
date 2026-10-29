"""Convert gradient-accepted AM1 output to audited GAFF2/AM1-BCC candidates.

Uses the accepted SQM output directly: no second quantum optimization. ACPYPE's
documented charge balancing is recorded; candidates permit numerical controls
only, not physical-PC/adsorption inference.
"""
import argparse
import fcntl
import json
from pathlib import Path
import re
import shutil
import subprocess

from pc_length_parameters import ROOT, audit, preflight, sha
from run_pc_sqm_checkpointed import BASE as QM_BASE, amber_environment, parse_completed, write_json
from pc_qm_route import route
from sqm_antechamber_view import write_view

BASE = ROOT/'runs/pc-length-forcefield-v2'


def finalize(dp):
    mol, source = preflight(dp)
    active_base, worker = route(dp)
    completed_path = active_base/f'pc{dp}/completed.json'
    if not completed_path.exists():
        return None
    completed = json.loads(completed_path.read_text())
    sqm = ROOT/completed['final_output']
    if sha(sqm) != completed['output_sha256']:
        raise ValueError('Accepted QM output hash changed')
    if worker == 'pc_sqm_lbfgs.py':
        from pc_sqm_lbfgs import parse_singlepoint, validate_coordinates
        _, _, convergence = parse_singlepoint(sqm, mol.GetNumAtoms())
        validate_coordinates(sqm, completed['coordinates_A'], [a.GetAtomicNum() for a in mol.GetAtoms()])
        convergence['optimizer'] = 'SQM analytic gradient with SciPy L-BFGS-B'
    else:
        _, convergence = parse_completed(sqm, [a.GetAtomicNum() for a in mol.GetAtoms()])
    if not convergence['gradient_gate_passed']:
        raise ValueError('QM output fails actual gradient criterion')
    BASE.mkdir(exist_ok=True)
    out = BASE/f'pc{dp}'
    with (BASE/f'pc{dp}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (out/'parameter-audit.json').exists():
            result = json.loads((out/'parameter-audit.json').read_text())
            for file, digest in result['outputs_sha256'].items():
                if sha(out/f'PC{dp}.amb2gmx'/file) != digest:
                    raise ValueError('Existing finalized forcefield changed')
            return out/f'PC{dp}.amb2gmx'
        if out.exists():
            raise RuntimeError('Partial postprocessing preserved; inspect before recovery')
        out.mkdir(); name=f'PC{dp}'; env=amber_environment()
        shutil.copy2(sqm, out/'accepted-sqm.out')
        write_view(sqm,out/'accepted-sqm-antechamber.out',mol.GetNumAtoms())
        write_json(out/'source.json', {'accepted_QM': str(sqm.relative_to(ROOT)), 'sha256': sha(sqm), 'input': source})
        def run(command, label):
            with (out/f'{label}.log').open('x') as log:
                subprocess.run(command, cwd=out, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1800)
        run(['antechamber','-i','accepted-sqm-antechamber.out','-fi','sqmout','-o',f'{name}_bcc_gaff2.mol2',
             '-fo','mol2','-c','bcc','-at','gaff2','-nc','0','-m','1','-dr','no','-pf','n','-s','2'], 'bcc')
        run(['parmchk2','-i',f'{name}_bcc_gaff2.mol2','-f','mol2','-o',f'{name}_AC.frcmod','-s','gaff2'], 'parmchk')
        if re.search(r'ATTN|need\s+revision', (out/f'{name}_AC.frcmod').read_text(), re.I):
            raise ValueError('Unassigned bonded parameters require review')
        (out/'leap.in').write_text(f'source leaprc.gaff2\nloadAmberParams {name}_AC.frcmod\n'
                                   f'MOL = loadMol2 {name}_bcc_gaff2.mol2\ncheck MOL\n'
                                   f'saveAmberParm MOL {name}_AC.prmtop {name}_AC.inpcrd\nquit\n')
        run(['tleap','-f','leap.in'], 'tleap')
        leap=(out/'tleap.log').read_text()
        counters=re.findall(r'Errors =\s*(\d+); Warnings =\s*(\d+)', leap)
        if not counters or int(counters[-1][0]):
            raise ValueError('TLEaP error or missing completion summary')
        if int(counters[-1][1]):
            raise ValueError('TLEaP warnings preserved for review before topology export')
        run(['acpype','-p',f'{name}_AC.prmtop','-x',f'{name}_AC.inpcrd','-b',name,'-a','gaff2'], 'acpype-convert')
        folder=out/f'{name}.amb2gmx'
        topology=(folder/f'{name}_GMX.top').read_text()
        begin=re.search(r'^\s*\[\s*atomtypes\s*\]',topology,re.M)
        end=re.search(r'^\s*\[\s*system\s*\]',topology,re.M)
        if not begin or not end:
            raise ValueError('Cannot isolate molecule/types from generated topology')
        (folder/f'{name}_GMX.itp').write_text(topology[begin.start():end.start()])
        for filename in [f'{name}_bcc_gaff2.mol2',f'{name}_AC.frcmod']:
            shutil.copy2(out/filename, folder/filename)
        shutil.copy2(sqm, folder/'sqm.out')
        result=audit(dp, folder)
        # Trace the raw residual to the printed QM charges, instead of inventing
        # missing precision or treating automatic balancing as invisible.
        text=sqm.read_text().split('Mulliken Charge',1)[1].split('Total Mulliken Charge',1)[0]
        rows=[line.split() for line in text.splitlines() if re.match(r'^\s*\d+\s+[A-Z][a-z]?\s+[-\d.]',line)]
        if len(rows)!=mol.GetNumAtoms():
            raise ValueError('Incomplete printed QM charge table')
        printed_sum=sum(float(r[-1]) for r in rows)
        if abs(printed_sum-result['mol2_net_charge_e']) > mol.GetNumAtoms()*5e-7:
            raise ValueError('BCC residual not explained by printed QM charge conservation')
        if abs(printed_sum) > mol.GetNumAtoms()*.0005 + 1e-8:
            raise ValueError('Charge residual exceeds original QM rounding envelope')
        result.update(state='numerical-control candidate;physical forcefield not approved',
                      QM_convergence=convergence,
                      charge_review={'printed_QM_sum_e': printed_sum,
                                     'BCC_conserves_printed_QM_sum': True,
                                     'ACPYPE_conversion_reconstructed': True,
                                     'decision': 'retain documented ACPYPE convention for bounded numerical controls;review charge/model sensitivity before physical inference'},
                      topology_folder=str(folder.relative_to(ROOT)))
        write_json(out/'parameter-audit.json', result)
        return folder


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--dp', type=int, choices=[4,8,16], required=True)
    args=p.parse_args(); print(finalize(args.dp))
