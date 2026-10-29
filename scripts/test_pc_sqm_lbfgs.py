"""Check installed-SQM derivative sign/units and optimizer recovery evidence."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pc_sqm_lbfgs import Evaluator, optimize, validate_coordinates, parse_singlepoint
from run_pc_sqm_checkpointed import amber_environment, input_text, parse_completed
import subprocess


class AnalyticGradientTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='sqm-gradient-test-')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def test_water_and_carbon_gradients_match_finite_difference(self):
        for smiles,atom in [('O',0),('CO',0)]:
            mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
            self.assertEqual(AllChem.EmbedMolecule(mol,randomSeed=2026),0)
            out=self.root/smiles;out.mkdir()
            ev=Evaluator(out,[a.GetAtomicNum() for a in mol.GetAtoms()],
                         [a.GetSymbol()+str(a.GetIdx()+1) for a in mol.GetAtoms()])
            x=mol.GetConformer().GetPositions().ravel();_,gradient=ev.evaluate(x)
            h=1e-4
            for axis in range(3):
                j=3*atom+axis;p=x.copy();m=x.copy();p[j]+=h;m[j]-=h
                finite=(ev.evaluate(p)[0]-ev.evaluate(m)[0])/(2*h)
                self.assertLess(abs(finite-gradient[j]),2e-4)

    def test_water_minimum_matches_native_and_saved_coordinates(self):
        xyz=np.array([[0.,0.,0.],[1.1,0.,0.],[-.2,1.1,0.]])
        native=self.root/'native';native.mkdir()
        (native/'sqm.in').write_text(input_text([8,1,1],['O','H1','H2'],xyz,50))
        env=amber_environment()
        subprocess.run([env['AMBERHOME']+'/bin/sqm','-O','-i','sqm.in','-o','sqm.out'],
                       cwd=native,env=env,check=True,capture_output=True)
        _,reference=parse_completed(native/'sqm.out',[8,1,1])
        out=self.root/'lbfgs';out.mkdir()
        result=optimize(out,[8,1,1],['O','H1','H2'],xyz)
        self.assertTrue(result['gradient_gate_passed'])
        self.assertLess(abs(result['energy_kcal_mol']-reference['heat_of_formation_kcal_mol_printed']),1e-4)
        saved=np.load(out/result['evaluation']/'input-coordinates.npy')
        np.testing.assert_array_equal(saved,result['coordinates_A'])
        with self.assertRaisesRegex(ValueError,'coordinates differ'):
            validate_coordinates(out/result['evaluation']/'sqm.out',saved+.01,[8,1,1])

    def test_scaled_optimizer_reaches_same_physical_gradient_gate(self):
        xyz=np.array([[0.,0.,0.],[1.1,0.,0.],[-.2,1.1,0.]])
        out=self.root/'scaled';out.mkdir()
        result=optimize(out,[8,1,1],['O','H1','H2'],xyz,coordinate_scale=.125)
        self.assertLess(result['rms_gradient_kcal_mol_A'],.0005)
        self.assertLess(abs(result['energy_kcal_mol']-(-59.2288)),1e-4)

    def test_print_overflow_is_rejected_without_inventing_a_derivative(self):
        path=self.root/'overflow.out'
        path.write_text('SCF Converged\nQMMM: SCF Energy = -1.0 KCal/mol\n'
                        'Forces on QM atoms from SCF calculation\n'
                        'QMMM: Atm 1: ******************** 2.000000 -3.000000\nCalculation Completed\n')
        with self.assertRaisesRegex(ValueError,'gradient overflow'):
            parse_singlepoint(path,1)

    def test_charged_carboxylate_gradient_and_total_charge(self):
        mol=Chem.AddHs(Chem.MolFromSmiles('CC(=O)[O-]'))
        self.assertEqual(AllChem.EmbedMolecule(mol,randomSeed=2026),0)
        out=self.root/'anion';out.mkdir()
        ev=Evaluator(out,[a.GetAtomicNum() for a in mol.GetAtoms()],
                     [a.GetSymbol()+str(a.GetIdx()+1) for a in mol.GetAtoms()],qmcharge=-1)
        x=mol.GetConformer().GetPositions().ravel();_,gradient=ev.evaluate(x)
        h=1e-4;p=x.copy();m=x.copy();p[6]+=h;m[6]-=h
        finite=(ev.evaluate(p)[0]-ev.evaluate(m)[0])/(2*h)
        self.assertLess(abs(finite-gradient[6]),2e-4)
        self.assertIn('qmcharge=-1',(out/'evaluation-000001/sqm.in').read_text())


if __name__=='__main__':unittest.main()
