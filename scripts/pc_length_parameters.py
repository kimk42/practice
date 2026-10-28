"""Independent whole-chain AM1-BCC assignments; never concatenate PC4 charges.

Run each new chain once in its own directory. Existing partial results require
inspection/recovery, not an automatic rerun. Numerical gates do not validate
bulk PC, aqueous behavior, or the physical accuracy of GAFF2 analog terms.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from surface_validation import section_rows

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'runs/pc-length-parameters-v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_path(dp):
    folder = 'models/pc-graft' if dp == 4 else f'models/surface-revision-v2/pc{dp}'
    return ROOT / folder / f'bpa-pc-dp{dp}-phenolic-ends.mol'


def preflight(dp):
    path = source_path(dp)
    mol = Chem.MolFromMolFile(str(path), removeHs=False)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1 or Chem.GetFormalCharge(mol) != 0:
        raise ValueError('Input must be one neutral connected PC molecule')
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula != f'C{16*dp-1}H{14*dp+2}O{3*dp-1}' or mol.GetNumAtoms() != 33*dp:
        raise ValueError('PC repeat count/formula mismatch')
    ends = mol.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]-c'))
    carbonates = mol.GetSubstructMatches(Chem.MolFromSmarts('[C](=[O])([O]-c)[O]-c'))
    if len(ends) != 2 or len(carbonates) != dp-1:
        raise ValueError('Incorrect phenolic ends or carbonate links')
    xyz = mol.GetConformer().GetPositions()
    lengths = [np.linalg.norm(xyz[b.GetBeginAtomIdx()] - xyz[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    if not np.isfinite(xyz).all() or min(lengths) < .7 or max(lengths) > 2.0:
        raise ValueError('Nonfinite or invalid starting bond geometry')
    return mol, {
        'PC_BPA_units': dp, 'source': str(path.relative_to(ROOT)), 'sha256': sha(path),
        'atoms': mol.GetNumAtoms(), 'formula': formula, 'molecular_mass_g_mol': Descriptors.MolWt(mol),
        'phenolic_OH_ends': len(ends), 'carbonate_links': len(carbonates),
        'bond_lengths_A': [float(min(lengths)), float(max(lengths))],
        'charge_route': 'whole-chain AM1-BCC/GAFF2; no PC4 charge transfer',
        'scope': 'single initial conformation and assumed OH termini; not experimental PC molecular weight/end chemistry',
    }


def audit(dp, folder):
    mol, initial = preflight(dp)
    name = f'PC{dp}'
    itp, gro = folder / f'{name}_GMX.itp', folder / f'{name}_GMX.gro'
    mol2 = folder / f'{name}_bcc_gaff2.mol2'
    sqm = folder / 'sqm.out'
    if 'Calculation Completed' not in sqm.read_text():
        raise ValueError('No SQM completion marker')
    text = itp.read_text()
    atoms = list(section_rows(text, 'atoms'))
    types = {r[0]: r for r in section_rows(text, 'atomtypes')}
    n = mol.GetNumAtoms()
    if len(atoms) != n or int(gro.read_text().splitlines()[1]) != n:
        raise ValueError('Topology/coordinates atom count mismatch')
    if [int(r[0]) for r in atoms] != list(range(1, n+1)):
        raise ValueError('Topology atom numbering changed')
    sections, active = {}, None
    for line in mol2.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            active = line[9:]; sections[active] = []
        elif active and line.strip():
            sections[active].append(line.split())
    rows = sections['ATOM']
    if len(rows) != n or [int(r[0]) for r in rows] != list(range(1, n+1)):
        raise ValueError('MOL2 atom numbering/count changed')
    elements = {'c': 6, 'o': 8, 'h': 1}
    if [elements.get(r[5][0]) for r in rows] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
        raise ValueError('MOL2 atom order/elements changed')
    expected = {tuple(sorted((b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1))):
                ('ar' if b.GetIsAromatic() else str(int(b.GetBondTypeAsDouble())))
                for b in mol.GetBonds()}
    observed = {tuple(sorted(map(int, r[1:3]))): r[3] for r in sections['BOND']}
    if expected != observed:
        raise ValueError('MOL2 connectivity or bond order differs from input')
    if {tuple(sorted(map(int, r[:2]))) for r in section_rows(text, 'bonds')} != set(expected):
        raise ValueError('ITP bond graph differs from input')
    charges = np.array([float(r[6]) for r in atoms])
    qmol2 = np.array([float(r[8]) for r in rows])
    if not np.isfinite(charges).all() or not np.isfinite(qmol2).all():
        raise ValueError('Nonfinite charges')
    if abs(charges.sum()) > 1e-4 or round(float(qmol2.sum())) != 0:
        raise ValueError('Topology not neutral or raw charge rounds to a nonzero integer')
    # ACPYPE 2026.10.15 balances the residual on all atoms tied at the largest
    # absolute charge. Reconstruct this documented conversion, not an invented
    # uniform correction. A nontrivial correction stays a separate review gate.
    extreme = max(qmol2) if abs(max(qmol2)) >= abs(min(qmol2)) else min(qmol2)
    balanced_ids = np.flatnonzero(qmol2 == extreme)
    predicted = qmol2.copy()
    predicted[balanced_ids] -= (qmol2.sum() - round(float(qmol2.sum()))) / len(balanced_ids)
    if not np.allclose(charges, predicted, atol=5.1e-7, rtol=0):
        raise ValueError('ITP charges differ from documented ACPYPE balancing/serialization')
    charge_review = abs(float(qmol2.sum())) > 1e-4
    for atom, row, mrow in zip(mol.GetAtoms(), atoms, rows):
        if row[1] != mrow[5] or row[1] not in types or row[4] != mrow[1]:
            raise ValueError('Atom name/type mismatch or undefined atom type')
        if abs(float(row[7]) - atom.GetMass()) > .2:
            raise ValueError('Atom mass mismatch')
        if not all(np.isfinite(float(v)) for v in row[6:8]):
            raise ValueError('Nonfinite atomic parameters')
    for row in types.values():
        if not all(np.isfinite(float(v)) and float(v) > 0 for v in row[-2:]):
            raise ValueError('Invalid LJ parameters')
    frcmods = list(folder.glob('*.frcmod'))
    if not frcmods:
        raise ValueError('Missing frcmod evidence')
    for path in frcmods:
        if re.search(r'ATTN|need\s+revision', path.read_text(), re.I):
            raise ValueError('Unassigned forcefield terms: ' + path.name)
    # OH charge is reported as O+its bonded H, not as a neutral repeat fragment.
    end_charges = []
    for o, _ in mol.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]-c')):
        ids = [o] + [a.GetIdx() for a in mol.GetAtomWithIdx(o).GetNeighbors() if a.GetAtomicNum() == 1]
        end_charges.append({'atoms_1based': [i+1 for i in ids], 'OH_charge_e': float(charges[ids].sum())})
    return {
        'state': 'assignment audited; charge-balancing review required' if charge_review else
                 'parameter assignment numerical audit passed; MD/physical validation pending',
        'input': initial, 'atoms': n, 'net_charge_e': float(charges.sum()),
        'mol2_net_charge_e': float(qmol2.sum()), 'atom_type_counts': dict(Counter(r[1] for r in atoms)),
        'OH_end_charge_diagnostics': end_charges,
        'charge_conversion': {
            'method': 'ACPYPE balanceCharges: residual distributed equally over atoms tied at largest absolute raw charge',
            'balanced_atoms_1based': (balanced_ids+1).tolist(),
            'max_abs_change_e': float(np.max(np.abs(charges-qmol2))),
            'review_required_before_new_MD': charge_review,
            'no_files_modified_by_audit': True,
            'precision_scope': 'MOL2/ITP arithmetic and converter provenance, not recovered QM digits or physical accuracy',
        },
        'checks': ['SQM completion', 'atom order/elements', 'connectivity/bond order', 'ITP bond graph',
                   'topology neutrality', 'ACPYPE charge-balancing reconstruction', 'masses', 'defined LJ types', 'frcmod revision flags'],
        'outputs_sha256': {p.name: sha(p) for p in [itp, gro, mol2, sqm, *frcmods]},
        'physical_forcefield_approved': False,
        'remaining': ['strict GROMACS preprocessing', 'numerical MD control', 'condensed-phase and conformation assessment'],
    }


def run(dp):
    if dp not in (8, 16):
        raise ValueError('Only new PC8/PC16 assignments may be launched here')
    _, initial = preflight(dp)
    BASE.mkdir(exist_ok=True)
    with (BASE / f'pc{dp}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out = BASE / f'pc{dp}'
        out.mkdir()  # Refuse all existing work, including incomplete quantum outputs.
        source = source_path(dp)
        shutil.copy2(source, out / 'input.mol')
        state = {'state': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
                 'worker_pid': os.getpid(), 'input': initial, 'SQM_limit_seconds': 86400,
                 'BLAS_threads': 1, 'output': str(out.relative_to(ROOT)),
                 'automatic_retry': False, 'automatic_MD': False}
        def save():
            temp = out / 'status.tmp'
            temp.write_text(json.dumps(state, indent=2) + '\n')
            temp.replace(out / 'status.json')
        save()
        try:
            env = dict(os.environ, MD_SQM_MAX_SECONDS='86400')
            with (out / 'worker.log').open('x') as log:
                subprocess.run(['bash', str(ROOT / 'scripts/parameterize_candidate.sh'),
                                str(out / 'input.mol'), f'PC{dp}', str(out)],
                               cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            if sha(source) != initial['sha256'] or sha(out / 'input.mol') != initial['sha256']:
                raise ValueError('Input changed during parameterization')
            result = audit(dp, out / f'PC{dp}.acpype')
            (out / 'parameter-audit.json').write_text(json.dumps(result, indent=2) + '\n')
            state['state'] = 'charge_balancing_review_required' if result['charge_conversion']['review_required_before_new_MD'] else 'parameter_audit_passed'
        except Exception as error:
            state.update(state='blocked_preserved', reason=str(error))
            raise
        finally:
            state['finished_utc'] = datetime.now(timezone.utc).isoformat()
            save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['preflight', 'run', 'audit'])
    parser.add_argument('--dp', type=int, choices=[4, 8, 16], required=True)
    parser.add_argument('--folder', type=Path, help='Existing ACPYPE folder for audit only')
    args = parser.parse_args()
    if args.action == 'run':
        run(args.dp)
    elif args.action == 'audit':
        folder = args.folder or BASE / f'pc{args.dp}' / f'PC{args.dp}.acpype'
        print(json.dumps(audit(args.dp, folder), indent=2))
    else:
        print(json.dumps(preflight(args.dp)[1], indent=2))
