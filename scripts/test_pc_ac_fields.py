"""Four-character AC names must never shift the charge/type columns."""
from pathlib import Path
import tempfile
import unittest
from pc_bcc_from_equivalence import replace_ac_charges
from validate_equatom_qsort import ac_rows


class ACFieldsTests(unittest.TestCase):
    def test_adjacent_atom_and_residue_names(self):
        text=('ATOM      1  O1  MOL     1      -4.117   0.780  15.578 -0.505100        31\n'
              'ATOM    118  C100MOL     1       1.596   5.796  -9.163  0.041100        11\n'
              'ATOM    264  H114MOL     1       4.978  -2.890  -6.288  0.155500        91\n')
        with tempfile.TemporaryDirectory() as name:
            path=Path(name)/'atoms.ac';path.write_text(text)
            rows=ac_rows(path)
            self.assertEqual([r[2] for r in rows],['O1','C100','H114'])
            self.assertEqual([float(r[8]) for r in rows],[-.5051,.0411,.1555])
            self.assertEqual([r[9] for r in rows],['31','11','91'])
            path.write_text(replace_ac_charges(text,[.125,-.25,-.375]))
            changed=ac_rows(path)
            self.assertEqual([float(r[8]) for r in changed],[.125,-.25,-.375])
            for before,after in zip(rows,changed):
                self.assertEqual(before[:8],after[:8]);self.assertEqual(before[9:],after[9:])


if __name__=='__main__':
    unittest.main()
