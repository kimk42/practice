import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
from md.equilibrate_small_pe import density_g_cm3
class DensityUnitsTest(unittest.TestCase):
    def test_one_water_per_30_cubic_angstroms(self):
        # 1 water molecule occupies approximately 0.030 nm3 at ambient density.
        self.assertAlmostEqual(density_g_cm3(18.01528,.030),.99713,places=4)
    def test_constant_mass_compression_changes_density(self):
        a=density_g_cm3(1000,10)
        self.assertAlmostEqual(density_g_cm3(1000,8),a*1.25)
if __name__=='__main__':unittest.main()
