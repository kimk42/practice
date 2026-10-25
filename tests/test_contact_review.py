import unittest
from md.contact_review import summarize

class ContactReviewTests(unittest.TestCase):
 def test_identity_uses_saved_frames_not_finer_count_series(self):
  analysis={'job_id':'example','saved_frame_context_time_range_ps':[30.4,31.4],'geometry_rows':[{'context_time_ps':30.4,'contact_residues':['E:344']},{'context_time_ps':30.9,'contact_residues':[]},{'context_time_ps':31.4,'contact_residues':['E:344']}]};mapping=[{'source_chain':'E','source_label_seq':344,'source_auth_seq':'312','name':'LYS','output_residue':5,'artificial_cap':False}];summary=summarize(analysis,mapping,[{'context_time_ps':30.05,'contact_native_residue_count':1}]);r=summary['contacts'][0]
  self.assertEqual(r['observed_context_times_ps'],[30.4,31.4]);self.assertAlmostEqual(r['saved_frame_fraction'],2/3);self.assertEqual(r['source_auth_seq'],'312');self.assertEqual(summary['first_any_contact_count_observation_ps'],30.05);self.assertFalse(summary['binding_sites_identified'])
