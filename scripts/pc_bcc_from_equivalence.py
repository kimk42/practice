"""Installed BCC correction with prevalidated full-path charge equivalence.

Only the equivalence sorting route is replaced. Raw printed Mulliken charges,
full-path classes, installed BCC parameters and GAFF2 typing remain explicit.
"""
import json
from pathlib import Path
import re
import shutil

import numpy as np
from rdkit import Chem

from pc_length_parameters import ROOT, preflight, sha
from run_pc_sqm_checkpointed import amber_environment, write_json
from sqm_antechamber_view import write_view
from validate_equatom_qsort import ac_rows, average_charges, printed_charges
from managed_subprocess import run_group


def replace_ac_charges(text,charges):
    lines=[];i=0
    for line in text.splitlines():
        if line.startswith('ATOM'):
            match=re.search(r'\s+([-+]?\d+\.\d+)\s+\S+\s*$',line)
            if not match:
                raise ValueError('Unexpected AC charge record')
            # Preserve all original columns, names, types, bonds and coordinates.
            begin,end=match.span(1)
            line=line[:begin]+f'{charges[i]:.6f}'+line[end:];i+=1
        lines.append(line)
    if i!=len(charges):
        raise ValueError('Charge table length differs')
    return '\n'.join(lines)+'\n'


def mol2_sections(path):
    parts={};current=None
    for line in path.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            current=line[9:];parts[current]=[]
        elif line.strip() and current:
            parts[current].append(line.split())
    return parts


def validate_equivalence_input(dp,record,raw,typed):
    """Bind a saved equivalence table to the current molecule and QM charges."""
    mol,info=preflight(dp)
    if info['sha256']!=record['molecule']['sha256']:
        raise ValueError('Equivalence molecule source differs')
    rows=ac_rows(typed)
    if [int(r[1]) for r in rows]!=list(range(1,mol.GetNumAtoms()+1)):
        raise ValueError('Typed atom count/order differs')
    edges=[tuple(sorted(map(int,line.split()[2:4])))
           for line in typed.read_text().splitlines() if line.startswith('BOND')]
    expected={tuple(sorted([b.GetBeginAtomIdx()+1,b.GetEndAtomIdx()+1])) for b in mol.GetBonds()}
    if len(edges)!=len(expected) or set(edges)!=expected:
        raise ValueError('Typed bond graph differs')
    ranks=list(Chem.CanonicalRankAtoms(mol,breakTies=False,includeChirality=False))
    representatives=[ranks.index(r)+1 for r in ranks]
    if representatives!=record['representatives_1based']:
        raise ValueError('Full-path equivalence no longer matches the reviewed graph partition')
    averaged=average_charges(printed_charges(raw,len(rows)),representatives)
    if not np.array_equal(averaged,np.array(record['averaged_charges_e'])):
        raise ValueError('Saved equivalence charges differ from accepted QM charges')
    return averaged


def calculate(dp,out,lock_fd=None):
    validation=ROOT/'models/surface-revision-v2/equatom-qsort-validation-20261030.json'
    report=json.loads(validation.read_text())
    if not report['cases']['4']['reference_and_qsort_identical']:
        raise ValueError('Reference sorting parity not verified')
    record=report['cases'][str(dp)]
    source=ROOT/f'runs/pc-length-forcefield-v2/pc{dp}'
    raw=source/'accepted-sqm.out';typed=source/'ANTECHAMBER_AC.AC'
    expected=json.loads((source/'source.json').read_text())['sha256']
    if sha(raw)!=expected:
        raise ValueError('Accepted QM source changed')
    averaged=validate_equivalence_input(dp,record,raw,typed)
    out.mkdir()
    shutil.copy2(raw,out/'accepted-sqm.out')
    write_view(raw,out/'accepted-sqm-antechamber.out',len(record['averaged_charges_e']))
    pre=replace_ac_charges(typed.read_text(),averaged)
    (out/'full-path-averaged.ac').write_text(pre)
    env=amber_environment()
    def run(command,label):
        with (out/(label+'.log')).open('x') as log:
            run_group(command,cwd=out,env=env,log=log,timeout=600,
                      pass_fds=(lock_fd,) if lock_fd is not None else ())
    run(['am1bcc','-i','full-path-averaged.ac','-o','bcc.ac','-f','ac','-t','bcc',
         '-p',env['AMBERHOME']+'/dat/antechamber/BCCPARM.DAT','-s','2','-j','1'],'bcc')
    rows=ac_rows(out/'bcc.ac');original=ac_rows(typed)
    if [r[1:3] for r in rows]!=[r[1:3] for r in original]:
        raise ValueError('BCC atom names/order changed')
    charges=[float(r[8]) for r in rows]
    if not np.isfinite(charges).all() or abs(sum(charges)-record['raw_QM_sum_e'])>len(charges)*5e-7:
        raise ValueError('BCC charge conservation differs beyond serialized precision')
    (out/'gaff2-bcc.ac').write_text(replace_ac_charges(typed.read_text(),charges))
    # Retain the installed route's AC coordinate serialization. The unchanged
    # raw SQM output remains the optimized-coordinate/convergence authority.
    name=f'PC{dp}_bcc_gaff2.mol2'
    run(['antechamber','-i','gaff2-bcc.ac','-fi','ac','-o',name,'-fo','mol2',
         '-at','gaff2','-j','0','-eq','0','-dr','no','-pf','n','-s','2'],'convert')
    sections=mol2_sections(out/name);atoms=sections['ATOM']
    if [r[1] for r in atoms]!=[r[2] for r in original] or [r[5] for r in atoms]!=[r[9] for r in original]:
        raise ValueError('Conversion changed atom names or GAFF2 types')
    if not np.array_equal(np.array(charges),np.array([float(r[8]) for r in atoms])):
        raise ValueError('Conversion changed BCC charges')
    result={'source_QM_sha256':sha(raw),'source_typed_AC_sha256':sha(typed),
            'equivalence_review_sha256':sha(validation),'molecule':record['molecule'],
            'charge_route':'full-path equivalence with qsort,upstream averaging,installed am1bcc -j1;no path truncation or new QM',
            'raw_printed_QM_sum_e':record['raw_QM_sum_e'],'mol2_charge_sum_e':float(sum(charges)),
            'serialized_preBCC_sum_e':float(sum(float(r[8]) for r in ac_rows(out/'full-path-averaged.ac'))),
            'BCC_parameters_sha256':sha(Path(env['AMBERHOME'])/'dat/antechamber/BCCPARM.DAT'),
            'outputs_sha256':{p.name:sha(p) for p in [out/'full-path-averaged.ac',out/'bcc.ac',out/name]}}
    if dp==4:
        reference=mol2_sections(source/name)
        if sections['BOND']!=reference['BOND']:
            raise ValueError('PC4 reference bond table differs')
        for a,b in zip(atoms,reference['ATOM']):
            if a[:2]!=b[:2] or a[5:]!=b[5:] or not np.array_equal(np.array(a[2:5],float),np.array(b[2:5],float)):
                raise ValueError('Split BCC conversion differs from original PC4 MOL2')
        result['installed_PC4_final_MOL2_fields_identical']=True
    write_json(out/'bcc-route-audit.json',result)
    return result


if __name__=='__main__':
    print(json.dumps(calculate(4,ROOT/'runs/pc4-split-bcc-parity-v2'),indent=2))
