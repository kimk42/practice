"""Large path counts must be detected without enumerating billions of paths."""
import unittest
from pc_equivalence_scalability import path_lower_bound,require_supported_full_path_count
from pc_length_parameters import preflight


class PathScalabilityTests(unittest.TestCase):
    def test_end_to_end_bound_and_pc16_gate(self):
        for dp in (4,8,16):
            mol,_=preflight(dp);record=path_lower_bound(mol)
            self.assertEqual(record['rings_on_end_to_end_path'],2*dp)
            self.assertEqual(record['end_to_end_simple_paths_lower_bound'],4**dp)
        require_supported_full_path_count(4)
        require_supported_full_path_count(8)
        with self.assertRaisesRegex(ValueError,'unrestricted equivalence held'):
            require_supported_full_path_count(16)


if __name__=='__main__':unittest.main()
