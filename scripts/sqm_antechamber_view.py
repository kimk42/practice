"""Remove only verbose derivative records for Amber's older sqmout reader.

Keep the raw SQM output as the convergence authority. This format-only copy
preserves every coordinate, charge, energy and completion line verbatim.
"""
import re
from pc_length_parameters import sha
from run_pc_sqm_checkpointed import write_json


def write_view(source, target, natoms):
    text=source.read_text();lines=text.splitlines(keepends=True)
    headers=[line for line in lines if line.startswith('QMMM: Forces on QM atoms from SCF calculation')]
    derivatives=[line for line in lines if re.match(r'^QMMM: Atm\s+\d+:',line)]
    if headers and (len(headers)!=1 or len(derivatives)!=natoms):
        raise ValueError('Unexpected verbose derivative layout; no conversion view written')
    if derivatives and not headers:raise ValueError('Derivative rows lack a recognized header')
    filtered=[];in_final=False;blank_records=0
    for line in lines:
        if line.strip()=='Final Structure':in_final=True
        if line in headers or line in derivatives:continue
        if headers and in_final and re.fullmatch(r'QMMM:\s*',line):
            blank_records+=1;continue
        filtered.append(line)
    view=''.join(filtered)
    for marker in ('Mulliken Charge','Final Structure','Calculation Completed'):
        if marker not in view:raise ValueError('Missing SQM output block: '+marker)
    target.write_text(view)
    report={'source_sha256':sha(source),'view_sha256':sha(target),'removed_derivative_lines':len(headers)+len(derivatives),
            'removed_empty_diagnostic_records_after_final_coordinates':blank_records,
            'all_other_lines_preserved_verbatim':True,'purpose':'compatibility with Antechamber sqmout parser;not a convergence certificate',
            'raw_convergence_evidence':str(source)}
    write_json(target.with_suffix('.conversion.json'),report)
    return report
