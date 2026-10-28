"""Preserve original AM1-BCC outputs and derive explicitly normalized candidates."""
from decimal import Decimal, ROUND_HALF_EVEN
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/phosphate-charge-projection-v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project(charges, target, groups, phosphorus_index):
    """Equal-weight least-squares correction under total-charge/equivalence constraints.

    Groups are already equal in these source candidates. Serialization uses 12
    decimals and assigns only its sub-picoelectron residual to the unique P.
    The added decimal places are arithmetic precision, not recovered QM digits.
    """
    assert sorted(i for group in groups for i in group) == list(range(len(charges)))
    if any(len({charges[i] for i in group}) != 1 for group in groups):
        raise ValueError('Unexpected source asymmetry: do not silently symmetrize')
    correction = (Decimal(target) - sum(charges)) / Decimal(len(charges))
    if abs(correction) > Decimal('0.0005'):
        raise ValueError('Correction exceeds original atomic-output rounding envelope')
    result = [(q + correction).quantize(Decimal('0.000000000001'), rounding=ROUND_HALF_EVEN)
              for q in charges]
    result[phosphorus_index] += Decimal(target) - sum(result)
    assert sum(result) == Decimal(target)
    assert all(len({result[i] for i in group}) == 1 for group in groups)
    assert max(abs(a-b) for a, b in zip(charges, result)) <= Decimal('0.0005')
    return result, correction


def main():
    if OUT.exists():
        raise RuntimeError('Existing derived candidates preserved;inspect instead of overwriting')
    base = ROOT / json.loads((ROOT / 'runs/buffer-charge-active.json').read_text())['directory']
    specifications = [('dihydrogen_phosphate-', -1, [[0, 2], [1], [3, 4], [5, 6]], 1),
                      ('hydrogen_phosphate2-', -2, [[0, 2, 3], [1], [4], [5]], 1)]
    prepared = []
    for species, target, groups, p_index in specifications:
        source = base / species / 'candidate.mol2'
        lines = source.read_text().splitlines(keepends=True)
        start = next(i for i, line in enumerate(lines) if line.strip() == '@<TRIPOS>ATOM') + 1
        stop = next(i for i in range(start, len(lines)) if lines[i].startswith('@<TRIPOS>'))
        rows = [line.split() for line in lines[start:stop]]
        charges = [Decimal(row[8]) for row in rows]
        assert rows[p_index][5] == 'p5' and sum(charges) == Decimal(target)-Decimal('0.001')
        for group in groups:
            assert len({rows[i][5] for i in group}) == 1
        normalized, delta = project(charges, target, groups, p_index)
        result = list(lines)
        for i, value in enumerate(normalized):
            spans = list(re.finditer(r'\S+', lines[start+i]))
            token = spans[8]
            result[start+i] = lines[start+i][:token.start()] + f'{value:.12f}' + lines[start+i][token.end():]
            assert result[start+i].split()[:8] == rows[i][:8]
            assert result[start+i].split()[9:] == rows[i][9:]
        assert result[:start] == lines[:start] and result[stop:] == lines[stop:]
        prepared.append((species, ''.join(result), {
            'species': species, 'source_file': str(source.relative_to(ROOT)),
            'source_sha256': sha(source), 'source_QM_sha256': sha(source.parent/'sqm.out'),
            'expected_integer_charge_e': target, 'original_sum_e': str(sum(charges)),
            'projected_sum_e': str(sum(normalized)),
            'equivalence_groups_1based': [[i+1 for i in group] for group in groups],
            'equivalence_basis': 'phosphate nonprotonated O resonance equivalence;OH O and H groups;unique P;source BCC already equal within groups',
            'least_squares_equal_shift_e': str(delta),
            'maximum_absolute_change_e': str(max(abs(a-b) for a,b in zip(charges,normalized))),
            'charges': [{'atom_id': int(row[0]), 'name': row[1], 'type': row[5],
                         'original_e': str(q), 'projected_e': str(new)}
                        for row,q,new in zip(rows,charges,normalized)],
            'geometry_connectivity_atom_types_unchanged': True,
            'aqueous_parameters_validated': False,
        }))
    OUT.mkdir()
    for species, text, record in prepared:
        target = OUT / (species+'.mol2')
        target.write_text(text)
        record['derived_file'] = str(target.relative_to(ROOT))
        record['derived_sha256'] = sha(target)
        assert sha(ROOT/record['source_file']) == record['source_sha256']
    report = {
        'state': 'explicit constrained-charge candidates prepared;original quantum/BCC files preserved',
        'method': 'minimum equal-weight squared change with integer net charge and retained source-equivalence groups',
        'precision_note': '12 decimal serialization is arithmetic;no claim of recovering hidden QM precision',
        'species': [row for _,_,row in prepared],
        'MD_topologies_exported': False, 'physical_force_field_approved': False,
        'remaining': ['review phosphate resonance/type representation and bonded coverage',
                      'numerical solute-water and aqueous ion-pairing checks',
                      'pH populations/ionic-condition model or explicit sensitivity assumptions'],
    }
    (OUT/'projection-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'models/buffer-parameters/phosphate-charge-projection-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{'species':r['species'],'before':r['original_sum_e'],
                      'after':r['projected_sum_e'],'max_change_e':r['maximum_absolute_change_e']}
                     for r in report['species']],indent=2))


if __name__ == '__main__':
    main()
