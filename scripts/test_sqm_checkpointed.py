"""Installed-SQM restart/convergence checks on a small neutral-water probe."""
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from run_pc_sqm_checkpointed import amber_environment, input_text, parse_completed


class SQMCheckpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='sqm-checkpoint-test-')
        cls.addClassCleanup(cls.temp.cleanup)
        cls.root=Path(cls.temp.name); cls.env=amber_environment()
        cls.initial=np.array([[0.,0.,0.],[1.1,0.,0.],[-.2,1.1,0.]])
        cls.short=cls.run_probe('short',cls.initial,2)
        cls.full=cls.run_probe('full',cls.initial,50)
        xyz,_=parse_completed(cls.short,[8,1,1])
        cls.resumed=cls.run_probe('resumed',xyz,50)

    @classmethod
    def run_probe(cls,name,xyz,cycles):
        out=cls.root/name;out.mkdir()
        (out/'sqm.in').write_text(input_text([8,1,1],['O','H1','H2'],xyz,cycles))
        with (out/'launch.log').open('w') as log:
            subprocess.run([cls.env['AMBERHOME']+'/bin/sqm','-O','-i','sqm.in','-o','sqm.out'],
                           cwd=out,env=cls.env,stdout=log,stderr=subprocess.STDOUT,check=True)
        return out/'sqm.out'

    def test_false_convergence_marker_does_not_pass(self):
        self.assertIn('geometry converged',self.short.read_text())
        _,result=parse_completed(self.short,[8,1,1])
        self.assertGreater(result['rms_gradient_kcal_mol_A_printed'],1.)
        self.assertFalse(result['gradient_gate_passed'])

    def test_coordinate_restart_reaches_same_small_molecule_minimum(self):
        a,first=parse_completed(self.full,[8,1,1]);b,second=parse_completed(self.resumed,[8,1,1])
        self.assertTrue(first['gradient_gate_passed'] and second['gradient_gate_passed'])
        self.assertAlmostEqual(first['heat_of_formation_kcal_mol_printed'],second['heat_of_formation_kcal_mol_printed'],places=4)
        np.testing.assert_allclose(np.linalg.norm(a[:,None]-a[None,:],axis=2),
                                   np.linalg.norm(b[:,None]-b[None,:],axis=2),atol=.001,rtol=0)

    def test_changed_final_element_is_rejected(self):
        text=self.full.read_text();before,after=text.rsplit('Final Structure',1)
        after=after.replace('      O ', '      C ',1)
        path=self.root/'bad.out';path.write_text(before+'Final Structure'+after)
        with self.assertRaisesRegex(ValueError,'elements changed'):
            parse_completed(path,[8,1,1])


if __name__=='__main__':unittest.main()
