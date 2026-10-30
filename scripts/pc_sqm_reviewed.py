"""Validate an explicitly reviewed QM endpoint; never launch a new calculation."""
import argparse
import json
from pc_length_parameters import ROOT, preflight, sha
from pc_sqm_lbfgs import parse_singlepoint, validate_coordinates

BASE=ROOT/'runs/pc-length-sqm-reviewed-v4'


def validate(dp):
    path=BASE/f'pc{dp}/completed.json'
    if not path.exists():
        raise RuntimeError('No explicitly reviewed endpoint; no automatic fallback')
    completed=json.loads(path.read_text())
    report_path=ROOT/completed['independent_review']
    if sha(report_path)!=completed['independent_review_sha256']:
        raise ValueError('Independent review changed')
    review=json.loads(report_path.read_text())
    raw=ROOT/completed['final_output']
    if not review['numerical_endpoint_robustness_passed'] or sha(raw)!=review['accepted_output_sha256'] or sha(raw)!=completed['output_sha256']:
        raise ValueError('Endpoint not covered by passed independent review')
    mol,info=preflight(dp)
    if info!=review['molecule_source']:
        raise ValueError('Reviewed molecule source changed')
    _,_,gate=parse_singlepoint(raw,mol.GetNumAtoms())
    validate_coordinates(raw,completed['coordinates_A'],[a.GetAtomicNum() for a in mol.GetAtoms()])
    if not gate['gradient_gate_passed']:
        raise ValueError('Actual gradient criterion failed')
    return completed


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dp',type=int,choices=[8],required=True)
    print(validate(parser.parse_args().dp)['state'])
