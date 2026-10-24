import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rdkit.Chem import rdMolDescriptors,rdMolTransforms
from md.prepare_fluoropolymers import molecule,conformer,MOTIFS,normalize_charge
from workbench.research import place,rotation,pmf_plan,budget,registry

class MolecularCandidateTests(unittest.TestCase):
 def test_formulas_and_topology_shared_across_pvdf_motifs(self):
  for family,formula in [('PTFE','C24F50'),('PVDF','C24H26F24')]:
   m=molecule(family);self.assertEqual(rdMolDescriptors.CalcMolFormula(m),formula);self.assertEqual(m.GetNumAtoms(),74)
  m=molecule('PVDF');bonds=[(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType())) for b in m.GetBonds()]
  for phase in ['pvdf-alpha','pvdf-beta','pvdf-gamma']:
   c=conformer(m,phase);self.assertEqual(bonds,[(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType())) for b in c.GetBonds()]);carbon=[a.GetIdx() for a in c.GetAtoms() if a.GetSymbol()=='C']
   for i in range(len(carbon)-3):
    angle=rdMolTransforms.GetDihedralDeg(c.GetConformer(),*carbon[i:i+4]);target=MOTIFS[phase][i%len(MOTIFS[phase])];self.assertAlmostEqual(((angle-target+180)%360)-180,0,places=7)
 def test_charge_correction_preserves_original_and_rejects_large_residual(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'polymer.mol2';original='@<TRIPOS>ATOM\n1 C1 0 0 0 c3 1 POL 0.104\n2 F1 1 0 0 f 1 POL -0.100\n@<TRIPOS>BOND\n';p.write_text(original);info=normalize_charge(p);self.assertAlmostEqual(info['raw_charge_sum_e'],.004);self.assertEqual(p.with_suffix('.raw.mol2').read_text(),original);self.assertAlmostEqual(sum(float(x.split()[-1]) for x in p.read_text().splitlines()[1:3]),0)
   p.write_text(original.replace('0.104','0.200'))
   with self.assertRaises(ValueError):normalize_charge(p)

class ResearchGateTests(unittest.TestCase):
 def test_placements_are_proper_and_preserve_all_distances(self):
  x=np.array([[0.,0,0],[2,1,0],[.5,0,2],[1,-2,-.5]]);d=np.linalg.norm(x[:,None]-x[None,:],axis=2)
  for roll,pitch in [(0,0),(90,0),(180,0),(270,0),(0,90),(0,-90)]:
   y,r=place(x,np.arange(4),roll,pitch);self.assertGreater(np.linalg.det(r),.999999);np.testing.assert_allclose(np.linalg.norm(y[:,None]-y[None,:],axis=2),d,atol=1e-12);self.assertAlmostEqual(y[:,2].min(),.5)
 def test_resume_budget_excludes_carried_time(self):
  b=budget({'last':{'context_time_ps':6.2,'additional_time_ps':3.2},'elapsed_wall_seconds':100,'current_process_start_context_time_ps':6.1})
  self.assertAlmostEqual(b['observed_current_process_simulated_ps'],.1);self.assertAlmostEqual(b['observed_preparation_ns_per_wall_day'],.0864);self.assertEqual(b['pilot_matrix']['replicas'],90)
 def test_design_does_not_claim_phases_or_free_energy(self):
  self.assertTrue(all(not x['comparative_adsorption_production_allowed'] for x in registry().values()));p=pmf_plan();self.assertEqual(p['status'],'design_only_not_executable');self.assertEqual(p['MD_steps'],0);self.assertIsNone(p['free_energy_values']);self.assertEqual(len(p['displacement_grid_nm']),41)

class TrajectoryMetricTests(unittest.TestCase):
 def test_periodic_contacts_and_fitted_geometry(self):
  from md.research_metrics import diagnostics
  patch=np.array([[1.,1,1],[3,1,1],[1,3,1],[3,3,1],[1,1,1.1],[3,1,1.1],[1,3,1.1],[3,3,1.1]])
  protein=np.array([[1,1,1.4],[2,1,1.4],[3,1,1.4]]);water=np.array([[1,1,1.25]])
  xyz=np.vstack([patch,protein,water]);box=np.array([10.,10,10]);ca=np.array([8,9,10]);ph=ca;poly=np.arange(8);w=np.array([11]);labels=['A:1','A:2','A:3']
  a=diagnostics(xyz,box,xyz,box,ca,ph,poly,w,labels);self.assertEqual(a['contact_residue_count'],2);self.assertEqual(a['water_O_within_0_35nm_of_both_groups'],1);self.assertAlmostEqual(a['protein_axis_angle_to_fitted_patch_plane_deg'],0)
  shifted=xyz.copy();shifted[:8,0]+=10;b=diagnostics(shifted,box,xyz,box,ca,ph,poly,w,labels);self.assertEqual(a['contact_residues'],b['contact_residues']);self.assertAlmostEqual(a['patch_plane_fit_RMS_height_nm'],b['patch_plane_fit_RMS_height_nm'])
 def test_degenerate_plane_is_not_reported_as_a_valid_orientation(self):
  from md.research_metrics import diagnostics
  xyz=np.array([[1.,1,1],[2,1,1],[3,1,1],[1,1,1.3],[2,1,1.3],[3,1,1.3]])
  d=diagnostics(xyz,[10,10,10],xyz,[10,10,10],[3,4,5],[3,4,5],[0,1,2],[],['A:1','A:2','A:3']);self.assertFalse(d['patch_plane_fit_usable']);self.assertIsNone(d['protein_axis_angle_to_fitted_patch_plane_deg'])

class SolventRegionTests(unittest.TestCase):
 def test_domain_volume_and_periodic_translation(self):
  from md.solvent_audit import region_statistics
  xyz=np.array([[5.,5,5],[1,1,1],[9,9,9]]);a=region_statistics(xyz,[10,10,10],[0],[1,2],18.016)
  expected=1000-4*np.pi*.8**3/3
  self.assertLess(abs(a['solvent_region_volume_nm3']-expected),6*a['volume_MC_standard_error_nm3']);self.assertEqual(a['water_O_in_region'],2)
  xyz[1,0]+=10;b=region_statistics(xyz,[10,10,10],[0],[1,2],18.016);self.assertEqual(a,b)
