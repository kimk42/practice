"""Regression checks for scientific gates, using disposable PC4 copies only."""
from pathlib import Path
import shutil
import tempfile
import unittest

from pc_length_parameters import ROOT, audit, sha


class ParameterAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='pc-audit-test-')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        source = ROOT / 'runs/parameterization-long/pc4/PC4.acpype'
        for name in ['PC4_GMX.itp', 'PC4_GMX.gro', 'PC4_bcc_gaff2.mol2', 'PC4_AC.frcmod', 'sqm.out']:
            shutil.copy2(source / name, self.folder / name)

    def test_archived_pc4_conversion_is_reported_without_mutation(self):
        before = {p.name: sha(p) for p in self.folder.iterdir()}
        result = audit(4, self.folder)
        self.assertAlmostEqual(result['mol2_net_charge_e'], -.004)
        self.assertTrue(result['charge_conversion']['review_required_before_new_MD'])
        self.assertEqual(result['charge_conversion']['balanced_atoms_1based'], [35, 54, 73])
        self.assertEqual(before, {p.name: sha(p) for p in self.folder.iterdir()})

    def test_charge_corruption_is_rejected_even_if_total_charge_is_unchanged(self):
        path = self.folder / 'PC4_GMX.itp'
        text = path.read_text().replace('-0.497600', '-0.487600', 1).replace('0.122600', '0.112600', 1)
        path.write_text(text)
        with self.assertRaisesRegex(ValueError, 'documented ACPYPE balancing'):
            audit(4, self.folder)

    def test_bond_order_corruption_is_rejected(self):
        path = self.folder / 'PC4_bcc_gaff2.mol2'
        parts = path.read_text().split('@<TRIPOS>BOND\n')
        lines = parts[1].splitlines()
        row = lines[0].split(); row[3] = '2' if row[3] == '1' else '1'
        lines[0] = ' '.join(row)
        path.write_text(parts[0] + '@<TRIPOS>BOND\n' + '\n'.join(lines) + '\n')
        with self.assertRaisesRegex(ValueError, 'bond order differs'):
            audit(4, self.folder)

    def test_missing_term_flag_is_rejected(self):
        with (self.folder / 'PC4_AC.frcmod').open('a') as stream:
            stream.write('\nATTN: need revision\n')
        with self.assertRaisesRegex(ValueError, 'Unassigned forcefield terms'):
            audit(4, self.folder)


if __name__ == '__main__':
    unittest.main()
