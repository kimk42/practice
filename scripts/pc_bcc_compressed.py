"""Use validated exact compressed path multiplicities before installed BCC."""
import json
from pathlib import Path

import numpy as np

from managed_subprocess import run_group
from pc_bcc_from_equivalence import apply_reviewed_bcc,mol2_sections
from pc_compressed_paths import score_histograms,representatives
from pc_length_parameters import ROOT,preflight,sha
from run_pc_sqm_checkpointed import amber_environment,write_json
from validate_equatom_qsort import average_charges,printed_charges
from sqm_antechamber_view import write_view

VALIDATION=ROOT/'models/surface-revision-v2/pc-compressed-path-validation.json'
PARITY=ROOT/'models/surface-revision-v2/pc-compressed-bcc-parity.json'


def reviewed_classes(dp,require_bcc_parity=True):
    report=json.loads(VALIDATION.read_text())
    if sha(ROOT/'scripts/pc_compressed_paths.py')!=report['compressed_algorithm_sha256']:
        raise ValueError('Compressed-path implementation changed since validation')
    for small in ('4','8'):
        if not report['cases'][small]['every_exact_binary64_score_and_integer_multiplicity_matches_enumeration']:
            raise ValueError('Exact upstream path-score comparison not complete')
    if require_bcc_parity:
        parity=json.loads(PARITY.read_text())
        if parity['path_validation_sha256']!=sha(VALIDATION):
            raise ValueError('BCC parity refers to different path validation')
        if not all(parity['cases'][small]['final_MOL2_fields_and_bonds_identical'] for small in ('4','8')):
            raise ValueError('Installed BCC parity is incomplete')
    mol,info=preflight(dp);record=report['cases'][str(dp)]
    if info['sha256']!=record['molecule']['sha256']:
        raise ValueError('Molecule differs from reviewed graph')
    histograms=score_histograms(mol);ids=representatives(histograms)
    if ids!=record['representatives_1based'] or [sum(h.values()) for h in histograms]!=record['unrestricted_path_counts']:
        raise ValueError('Recomputed full-path signatures differ')
    return mol,info,ids


def calculate(dp,source,out,lock_fd=None,require_bcc_parity=True):
    mol,info,ids=reviewed_classes(dp,require_bcc_parity)
    typed=source/'typed-no-charge';typed.mkdir()
    raw=source/'accepted-sqm.out'
    write_view(raw,typed/'accepted-view.out',mol.GetNumAtoms())
    with (typed/'typing.log').open('x') as log:
        run_group(['antechamber','-i','accepted-view.out','-fi','sqmout','-o','typed.ac','-fo','ac',
                   '-at','gaff2','-j','4','-eq','0','-dr','no','-pf','n','-s','2'],
                  cwd=typed,env=amber_environment(),log=log,timeout=600,
                  pass_fds=(lock_fd,) if lock_fd is not None else ())
    if (typed/'sqm.in').exists() or (typed/'sqm.out').exists():
        raise ValueError('Typing-only conversion unexpectedly produced a QM calculation')
    charges=printed_charges(raw,mol.GetNumAtoms());averaged=average_charges(charges,ids)
    record={'molecule':info,'representatives_1based':ids,'averaged_charges_e':averaged.tolist(),
            'raw_QM_sum_e':float(charges.sum()),'path_validation_sha256':sha(VALIDATION),
            'raw_QM_sha256':sha(raw),'no_paths_truncated':True,'new_QM_performed':False}
    evidence=source/'compressed-equivalence-audit.json';write_json(evidence,record)
    result=apply_reviewed_bcc(dp,out,record,raw,typed/'typed.ac',evidence,lock_fd,
                             'exact compressed unrestricted path-score multiplicities,upstream charge averaging,installed am1bcc -j1;no path truncation or new QM')
    return result


def validate_parity():
    import shutil
    base=ROOT/'runs/pc-compressed-bcc-parity-v1';base.mkdir();results={}
    for dp in (4,8):
        source=base/f'pc{dp}';source.mkdir()
        original=ROOT/f'runs/pc-length-forcefield-v2/pc{dp}'
        shutil.copy2(original/'accepted-sqm.out',source/'accepted-sqm.out')
        calculate(dp,source,source/'candidate',require_bcc_parity=False)
        parameter=json.loads((original/'parameter-audit.json').read_text())
        reference=ROOT/parameter['topology_folder']/f'PC{dp}_bcc_gaff2.mol2'
        generated=source/'candidate'/f'PC{dp}_bcc_gaff2.mol2'
        a=mol2_sections(reference);b=mol2_sections(generated)
        if a['ATOM']!=b['ATOM'] or a['BOND']!=b['BOND']:
            raise ValueError('Compressed route MOL2 differs from validated candidate')
        results[str(dp)]={'final_MOL2_fields_and_bonds_identical':True,
                          'reference_MOL2_sha256':sha(reference),'generated_MOL2_sha256':sha(generated),
                          'typed_without_new_QM':True}
    report={'cases':results,'path_validation_sha256':sha(VALIDATION),
            'decision':'Exact compressed full-path route may prepare PC16 numerical candidates only after its unchanged QM gradient gate passes;physical approval remains false.',
            'physical_approval':False,'PC16_charge_assigned':False}
    write_json(PARITY,report);print(json.dumps(report,indent=2))


if __name__=='__main__':validate_parity()
