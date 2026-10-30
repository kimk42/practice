"""Prevent changed Hamiltonians or silent restarts in the dry-PC extension."""
from pathlib import Path
import tempfile
import unittest
from continue_pc_npt100 import check_tpr_comparison,resume_command


class ContinuationTests(unittest.TestCase):
    def test_only_step_count_may_change(self):
        good='header\ncomparing inputrec\ninputrec->nsteps (40000 - 200000)\ncomparing mtop topology\n'
        check_tpr_comparison(good)
        with self.assertRaisesRegex(ValueError,'not limited'):
            check_tpr_comparison(good+'inputrec->ref-t[0] (300 - 310)\n')
        with self.assertRaises(ValueError):check_tpr_comparison('incomplete comparison')

    def test_partial_run_never_falls_back_to_original_checkpoint(self):
        with tempfile.TemporaryDirectory() as name:
            out=Path(name);(out/'source.cpt').touch()
            self.assertEqual(resume_command(out)[-3:],['-cpi','source.cpt','-noappend'])
            (out/'npt.part0002.log').touch()
            with self.assertRaisesRegex(RuntimeError,'without checkpoint'):resume_command(out)
            (out/'npt.cpt').touch()
            self.assertEqual(resume_command(out)[-3:],['-cpi','npt.cpt','-append'])


if __name__=='__main__':unittest.main()
