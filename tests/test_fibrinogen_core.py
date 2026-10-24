import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from md.prepare_fibrinogen_core import internal_atom,methylamide_cap
class CapGeometryTest(unittest.TestCase):
    def test_internal_coordinate_bond_and_angle(self):
        a=np.array([0.,1.,0.]);b=np.zeros(3);c=np.array([1.,0.,0.])
        d=internal_atom(a,b,c,1.329,121.7,45)
        self.assertAlmostEqual(np.linalg.norm(d-c),1.329)
        angle=np.rad2deg(np.arccos(np.dot((b-c)/np.linalg.norm(b-c),(d-c)/np.linalg.norm(d-c))))
        self.assertAlmostEqual(angle,121.7)
    def test_methylamide_preserves_peptide_geometry(self):
        r={'CA':{'xyz':np.array([0.,0.,0.])},'C':{'xyz':np.array([1.5,0.,0.])},'O':{'xyz':np.array([2.1,1.,0.])}}
        cap=methylamide_cap(r)
        self.assertAlmostEqual(np.linalg.norm(cap['N']-r['C']['xyz']),1.335)
        self.assertAlmostEqual(np.linalg.norm(cap['C']-cap['N']),1.458)
        # The new amide nitrogen lies opposite the existing carbonyl oxygen.
        self.assertLess(cap['N'][1],r['C']['xyz'][1])
if __name__=='__main__':unittest.main()
