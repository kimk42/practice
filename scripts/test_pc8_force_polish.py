"""Verify bounded force relaxation on a stiff analytic potential."""
import unittest
import numpy as np
from pc8_force_polish import fire_step


class ForcePolishTests(unittest.TestCase):
    def test_stiff_harmonic_minimum_and_step_bound(self):
        x=np.array([[.1,.01,-.02]])
        v=np.zeros_like(x);dt=.005;alpha=.1;positive=0
        stiffness=np.array([[1.,100.,3000.]])
        for _ in range(5000):
            gradient=stiffness*x
            if np.linalg.norm(gradient)<1e-6:
                break
            dx,v,dt,alpha,positive=fire_step(gradient,v,dt,alpha,positive)
            self.assertLessEqual(np.linalg.norm(dx,axis=1).max(),.00100000001)
            x+=dx
        self.assertLess(np.linalg.norm(stiffness*x),1e-6)


if __name__ == '__main__':
    unittest.main()
