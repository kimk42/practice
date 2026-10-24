import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from md.run_pe_complex import surface_metrics
from md.build_pe_complex import principal_orientation
class PEMetricsTest(unittest.TestCase):
    def test_periodic_contact_and_distinct_residue_count(self):
        xyz=np.array([[.1,1,1],[.2,1,1],[1.5,1,1],[9.9,1,1]])
        d,res=surface_metrics(xyz,np.array([0,1,2]),np.array([3]),{0:143,1:143,2:144},np.array([10.,10.,10.]))
        self.assertAlmostEqual(d,.2)
        self.assertEqual(res,[143])
    def test_orientation_preserves_pair_distances_and_chirality(self):
        xyz=np.array([[0.,0.,0.],[5.,1.,0.],[2.,0.,1.],[3.,-.5,-1.]])
        rotated=principal_orientation(xyz.tolist())
        np.testing.assert_allclose(np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2),np.linalg.norm(rotated[:,None]-rotated[None,:],axis=2),atol=1e-12)
        before=np.linalg.det(xyz[1:]-xyz[0]);after=np.linalg.det(rotated[1:]-rotated[0])
        self.assertGreater(before*after,0)
if __name__=='__main__':unittest.main()
