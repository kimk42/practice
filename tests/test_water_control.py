import unittest
import numpy as np
from md.water_control import initial_box,summarize

class WaterControlTests(unittest.TestCase):
 def test_box_density_and_rigid_geometry(self):
  length,xyz=initial_box();self.assertEqual(xyz.shape,(1536,3));self.assertAlmostEqual(512*18.016*.00166053906660/length**3,.86)
  mol=xyz.reshape(512,3,3);np.testing.assert_allclose(np.linalg.norm(mol[:,1]-mol[:,0],axis=1),.09572,atol=1e-12);np.testing.assert_allclose(np.linalg.norm(mol[:,2]-mol[:,0],axis=1),.09572,atol=1e-12)
  np.testing.assert_array_equal(xyz,initial_box()[1]);self.assertFalse(np.allclose(xyz,initial_box(seed=709)[1]))
 def test_low_density_control_cannot_trigger_large_continuation(self):
  rows=[{'ensemble':'NPT','additional_NPT_ps':i*.5,'density_g_cm3':.86,'volume_nm3':18+i*.0001,'temperature_K':310} for i in range(1,201)]
  self.assertEqual(summarize(rows)['status'],'needs_review')
  for r in rows:r['density_g_cm3']=.99
  result=summarize(rows);self.assertEqual(result['status'],'technical_control_passed');self.assertFalse(result['equilibrium_established'])
  for r in rows:r['volume_nm3']=18
  self.assertEqual(summarize(rows)['status'],'needs_review')
