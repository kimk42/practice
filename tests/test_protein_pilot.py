import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from md.protein_pilot import aligned_rmsd
class RMSDTest(unittest.TestCase):
    def test_alignment_removes_rigid_motion_but_retains_deformation(self):
        a=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        rot=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
        b=a@rot+np.array([3.,4.,5.])
        self.assertLess(aligned_rmsd(b,a),1e-12)
        b[0]+=np.array([.4,.2,.1])
        self.assertGreater(aligned_rmsd(b,a),.05)
