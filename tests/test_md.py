"""Engine integration tests use a synthetic two-particle system, not a polymer model."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from md.run import run, validate
from openmm import System, HarmonicBondForce, XmlSerializer, Vec3, unit
from openmm.app import Topology, PDBFile, element

class MDTest(unittest.TestCase):
    def config(self):
        return {'surface':{'polymer':'PVDF','phase':'beta','face':'(100)','polarization':'plus-normal'},'temperature_K':310,'timestep_fs':1,'steps':10,'seed':42}

    def test_reject_phase_and_polarization(self):
        c=self.config(); c['surface']['phase']='delta'
        with self.assertRaises(ValueError): validate(c)
        c['surface']['phase']='alpha'
        with self.assertRaises(ValueError): validate(c)

    def test_md_outputs_and_model_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); bundle=root/'bundle'; bundle.mkdir()
            system=System(); system.addParticle(12); system.addParticle(12)
            force=HarmonicBondForce(); force.addBond(0,1,0.15,1000); system.addForce(force)
            (bundle/'system.xml').write_text(XmlSerializer.serialize(system))
            top=Topology(); chain=top.addChain(); res=top.addResidue('TST',chain)
            a=top.addAtom('C1',element.carbon,res); b=top.addAtom('C2',element.carbon,res); top.addBond(a,b)
            with (bundle/'coordinates.pdb').open('w') as f:
                PDBFile.writeFile(top,[Vec3(0,0,0),Vec3(.16,0,0)]*unit.nanometer,f)
            c=self.config(); config=root/'config.json'; config.write_text(json.dumps(c))
            manifest={'surface':c['surface'],'sha256':{n:hashlib.sha256((bundle/n).read_bytes()).hexdigest() for n in ('system.xml','coordinates.pdb')}}
            for key in ('protein_source','polymer_forcefield','protein_forcefield','water_model','preparation_notes'):
                manifest[key]='SYNTHETIC ENGINE TEST ONLY; no scientific model validation'
            (bundle/'manifest.json').write_text(json.dumps(manifest))
            result=run(config,bundle,root/'output')
            self.assertEqual(result['status'],'completed')
            for name in ('energy.csv','trajectory.dcd','final.pdb','result.json'):
                self.assertGreater((root/'output'/name).stat().st_size,0)
            c['surface']['phase']='gamma'; config.write_text(json.dumps(c))
            with self.assertRaisesRegex(ValueError,'differs'): run(config,bundle,root/'rejected')
            c['surface']['phase']='beta'; config.write_text(json.dumps(c))
            with (bundle/'system.xml').open('a') as f: f.write('\n')
            with self.assertRaisesRegex(ValueError,'checksum'): run(config,bundle,root/'tampered')

if __name__=='__main__': unittest.main()
