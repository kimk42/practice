import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from md.prepare_protein import fetch, audit
import pdbfixer

class ProteinPreparationTest(unittest.TestCase):
    def test_audit_preserves_source_and_detects_missing_atom(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(pdbfixer.__file__).parent/'templates/ALA.pdb'
            text=source.read_text()
            modified='\n'.join(line for line in text.splitlines() if not (line.startswith('ATOM') and line[12:16].strip()=='CB'))+'\n'
            p=Path(tmp)/'incomplete.pdb';p.write_text(modified)
            original=p.read_bytes()
            result=audit(p,Path(tmp)/'audit.json')
            self.assertEqual(result['modifications_performed'],[])
            self.assertTrue(any('CB' in r['atoms'] for r in result['missing_heavy_atoms']))
            self.assertEqual(p.read_bytes(),original)
            self.assertEqual(result['sha256'],hashlib.sha256(original).hexdigest())
    def test_download_failure_is_recorded_not_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('urllib.request.urlopen',side_effect=OSError('network blocked')):
                report=fetch('3GHG',Path(tmp)/'structure')
            self.assertEqual(report['status'],'blocked')
            self.assertFalse((Path(tmp)/'structure/source.pdb').exists())
            self.assertEqual(json.loads((Path(tmp)/'structure/download.json').read_text())['status'],'blocked')

if __name__=='__main__':unittest.main()
