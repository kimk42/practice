"""Verify the installed Antechamber can consume a lossless compatibility view."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from pc_sqm_lbfgs import Evaluator
from pc_length_parameters import sha
from run_pc_sqm_checkpointed import amber_environment
from sqm_antechamber_view import write_view


class SQMViewTests(unittest.TestCase):
    def test_charged_carboxylate_bcc_without_altering_raw_qm(self):
        with tempfile.TemporaryDirectory(prefix='sqm-view-test-') as temp:
            root=Path(temp);mol=Chem.AddHs(Chem.MolFromSmiles('CC(=O)[O-]'))
            self.assertEqual(AllChem.EmbedMolecule(mol,randomSeed=2026),0)
            ev=Evaluator(root,[a.GetAtomicNum() for a in mol.GetAtoms()],
                         [a.GetSymbol()+str(a.GetIdx()+1) for a in mol.GetAtoms()],qmcharge=-1)
            ev.evaluate(mol.GetConformer().GetPositions().ravel())
            source=root/'evaluation-000001/sqm.out';original=source.read_text();before=sha(source)
            target=root/'compatible.out';report=write_view(source,target,mol.GetNumAtoms())
            self.assertEqual(before,sha(source))
            self.assertEqual(report['removed_derivative_lines'],mol.GetNumAtoms()+1)
            self.assertEqual(original.split('Mulliken Charge')[1].split('Total Mulliken Charge')[0],
                             target.read_text().split('Mulliken Charge')[1].split('Total Mulliken Charge')[0])
            command=['antechamber','-i','compatible.out','-fi','sqmout','-o','candidate.mol2','-fo','mol2',
                     '-c','bcc','-at','gaff2','-j','4','-nc','-1','-m','1','-eq','1','-dr','no','-pf','n']
            result=subprocess.run(command,cwd=root,env=amber_environment(),capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            atoms=[r.split() for r in (root/'candidate.mol2').read_text().split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0].splitlines() if r.strip()]
            self.assertEqual(len(atoms),mol.GetNumAtoms())
            self.assertLess(abs(sum(float(r[8]) for r in atoms)+1),.004)


if __name__=='__main__':unittest.main()
