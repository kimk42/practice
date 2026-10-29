"""Recovery must use saved poses and leave all prior artifacts intact."""
import json
from pathlib import Path
import tempfile
import unittest
from pc_length_parameters import ROOT, sha
from run_pc_packing_margin import recover_attempt


class PackingRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='packing-recovery-test-',dir=ROOT/'runs')
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        for name,text in {'molecule.pdb':'retained molecule\n', 'packing.restart':'1 2 3 4 5 6\n',
                          'packmol.inp':'restart_from old.restart\ntolerance 2.004\n',
                          'packmol.log':'interrupted optimization\n'}.items():
            (self.root/name).write_text(text)

    def test_recovery_preserves_prior_files_and_reuses_latest_poses(self):
        before={p:sha(p) for p in self.root.iterdir()}
        out,case=recover_attempt(self.root,{'state':'running'},{'chains':1})
        self.assertEqual(before,{p:sha(p) for p in before})
        self.assertEqual(sha(out/'source.restart'),sha(self.root/'packing.restart'))
        self.assertTrue((out/'packmol.inp').read_text().startswith('restart_from source.restart\n'))
        self.assertEqual(case['input_sha256'],sha(out/'packmol.inp'))
        self.assertEqual(json.loads((out/'previous-status.json').read_text())['state'],'running')

    def test_incomplete_pose_checkpoint_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError,'no complete saved poses'):
            recover_attempt(self.root,{'state':'running'},{'chains':2})

    def test_reviewed_local_polish_changes_only_optimizer_controls(self):
        (self.root/'packmol.inp').write_text('tolerance 2.004\nprecision 0.01\nmaxit 200\npbc 0 0 0 50 50 50\n')
        before={p:sha(p) for p in self.root.iterdir()}
        out,case=recover_attempt(self.root,{'state':'reviewed_local_polish','recovery_method':'local_polish'},{'chains':1})
        self.assertEqual(before,{p:sha(p) for p in before})
        text=(out/'packmol.inp').read_text()
        for line in ('tolerance 2.004','precision 0.01','pbc 0 0 0 50 50 50','maxit 1000','disable_movebad'):
            self.assertIn(line,text)
        self.assertEqual(case['restart']['method'],'local_polish')

    def test_relocation_recovery_preserves_generator_and_distance(self):
        (self.root/'packmol.inp').write_text('disable_movebad\ntolerance 2.004\nmaxit 1000\n')
        (self.root/'case.json').write_text(json.dumps({'chains':1,'generator_binary':'verified-variant','generator_sha256':'retained-digest'}))
        out,case=recover_attempt(self.root,{'state':'reviewed_relocation','recovery_method':'resume_with_relocation'},{'chains':1})
        self.assertNotIn('disable_movebad',(out/'packmol.inp').read_text())
        self.assertIn('tolerance 2.004',(out/'packmol.inp').read_text())
        self.assertEqual(case['generator_sha256'],'retained-digest')


if __name__=='__main__':unittest.main()
