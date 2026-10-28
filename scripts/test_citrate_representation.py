"""Reject chemical/protonation changes masquerading as resonance serialization."""
from pathlib import Path
import tempfile
import unittest
from audit_citrate_representation import review
from pc_length_parameters import ROOT,sha


class CitrateReviewTests(unittest.TestCase):
    def test_existing_candidates_are_unchanged_and_microstate_is_explicit(self):
        for name in ('citrate3-','hydrogen_citrate2-'):
            path=ROOT/'runs/buffer-charge-candidates-g5tzm8pz'/name/'candidate.mol2'
            before=sha(path);r=review(name)
            self.assertEqual(before,sha(path));self.assertFalse(r['approved_for_full_buffer_MD'])
            if name.startswith('hydrogen'):
                self.assertEqual(r['protonated_carboxyl_sites'][0]['site'],'terminal')

    def test_loss_of_acid_double_bond_is_rejected(self):
        source=ROOT/'runs/buffer-charge-candidates-g5tzm8pz/hydrogen_citrate2-/candidate.mol2'
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'changed.mol2';text=source.read_text()
            before,after=text.split('@<TRIPOS>BOND',1)
            rows=after.splitlines();fields=rows[1].split();self.assertEqual(fields[3],'2')
            fields[3]='1';rows[1]=' '.join(fields)
            path.write_text(before+'@<TRIPOS>BOND'+'\n'.join(rows))
            with self.assertRaisesRegex(ValueError,'outside carboxylate'):
                review('hydrogen_citrate2-',path)


if __name__=='__main__':unittest.main()
