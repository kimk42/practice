import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from md.project_paths import amber_home, resolve_project_path
from scripts.restore_raw_results import restore


class RepositoryToolsTest(unittest.TestCase):
    def archive(self, path, entries, bad_hash=False):
        manifest = []
        with zipfile.ZipFile(path, 'w') as z:
            for name, data in entries.items():
                z.writestr(name, data)
                manifest.append({'path': name, 'bytes': len(data), 'sha256': '0'*64 if bad_hash else hashlib.sha256(data).hexdigest()})
            z.writestr('FILES.json', json.dumps({'files': manifest}))

    def test_restore_verifies_raw_inputs_and_keeps_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);dest=root/'checkout';dest.mkdir();(dest/'README.md').write_text('current code')
            self.archive(root/'result.zip', {'data/model.bin':b'input','runs/run/state.xml':b'state','README.md':b'old archive code'})
            first=restore(root/'result.zip',dest)
            self.assertEqual(first['restored'],2)
            self.assertEqual((dest/'README.md').read_text(),'current code')
            self.assertEqual(restore(root/'result.zip',dest)['already_present'],2)

    def test_corrupt_data_is_not_installed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.archive(root/'bad.zip',{'data/model.bin':b'input'},bad_hash=True)
            with self.assertRaisesRegex(ValueError,'integrity'):restore(root/'bad.zip',root/'checkout')
            self.assertFalse((root/'checkout/data/model.bin').exists())

    def test_lfs_pointer_has_actionable_restore_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);pointer=root/'results-complete.zip'
            pointer.write_bytes(b'version https://git-lfs.github.com/spec/v1\noid sha256:'+b'0'*64+b'\nsize 100\n')
            with self.assertRaisesRegex(ValueError,'git lfs pull'):restore(pointer,root/'checkout')

    def test_existing_different_data_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);dest=root/'checkout';(dest/'data').mkdir(parents=True);(dest/'data/model.bin').write_bytes(b'user data')
            self.archive(root/'result.zip',{'data/model.bin':b'archived data'})
            with self.assertRaisesRegex(ValueError,'refusing overwrite'):restore(root/'result.zip',dest)
            self.assertEqual((dest/'data/model.bin').read_bytes(),b'user data')

    def test_archive_paths_cannot_escape_checkout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.archive(root/'bad.zip',{'data/../../outside':b'bad'})
            with self.assertRaisesRegex(ValueError,'Unsafe'):restore(root/'bad.zip',root/'checkout')
            self.assertFalse((root/'outside').exists())

    def test_relocated_local_dataset_takes_priority(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);old=root/'old/runs/pilot/state.xml';new=root/'new/runs/pilot/state.xml'
            for p in [old,new]:p.parent.mkdir(parents=True);p.write_text(p.parent.parent.parent.name)
            self.assertEqual(resolve_project_path(old,root/'new'),new)
            self.assertEqual(resolve_project_path('runs/pilot/state.xml',root/'new'),new)

    def test_explicit_amberhome_is_respected(self):
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'AMBERHOME':tmp}):
            self.assertEqual(amber_home(),Path(tmp).resolve())


if __name__=='__main__':unittest.main()
