"""Validate compressed path scores against direct enumeration on small graphs."""
from collections import Counter
import unittest

from rdkit import Chem

from pc_compressed_paths import score_histograms
from pc_bcc_compressed import reviewed_classes


def brute_force(mol):
    result=[]
    for atom in mol.GetAtoms():
        root=atom.GetIdx();histogram=Counter()
        def visit(path,score):
            histogram[score]+=1
            for neighbor in mol.GetAtomWithIdx(path[-1]).GetNeighbors():
                i=neighbor.GetIdx()
                if i not in path:
                    visit(path+(i,),score+((len(path)+1)*.11+neighbor.GetAtomicNum()*.08))
        visit((root,),.11+atom.GetAtomicNum()*.08);result.append(histogram)
    return result


class CompressedPathTests(unittest.TestCase):
    def test_exact_scores_with_nonuniform_elements_and_multiple_rings(self):
        for smiles in ('CC(O)CO','OC1CCNCC1O','Oc1ccc(Cc2ccc(O)cc2)cc1'):
            mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
            self.assertEqual(score_histograms(mol),brute_force(mol),smiles)

    def test_fused_rings_require_separate_review(self):
        mol=Chem.AddHs(Chem.MolFromSmiles('c1ccc2ccccc2c1'))
        with self.assertRaisesRegex(ValueError,'disjoint'):score_histograms(mol)

    def test_pc16_review_recomputes_large_integer_multiplicities(self):
        mol,_,ids=reviewed_classes(16)
        self.assertEqual(mol.GetNumAtoms(),528)
        self.assertEqual(len(set(ids)),153)


if __name__=='__main__':unittest.main()
