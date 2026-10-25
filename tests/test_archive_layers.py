import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import zipfile
from scripts.archive_layers import initialize,publish_snapshot,load_index
from scripts.restore_raw_results import restore


def archive(path,files):
 with zipfile.ZipFile(path,'w') as z:
  rows=[]
  for name,data in files.items():z.writestr(name,data);rows.append({'path':name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
  z.writestr('FILES.json',json.dumps({'files':rows}))

class LayerTests(unittest.TestCase):
 def test_latest_versions_restore_without_rewriting_baseline_or_code(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);data=root/'datasets';data.mkdir();base=data/'results-complete.zip';archive(base,{'data/a':b'old','runs/one':b'unchanged','README.md':b'old code'});original=base.read_bytes();initialize(data)
   full=root/'new.zip';archive(full,{'data/a':b'new','runs/one':b'unchanged','runs/two':b'new run','README.md':b'new code'});index=publish_snapshot(full,data)
   self.assertEqual(len(index['archives']),2);self.assertEqual(base.read_bytes(),original);latest={r['path']:r for r in index['files']};self.assertEqual(latest['runs/one']['archive'],base.name);self.assertNotEqual(latest['data/a']['archive'],base.name)
   dest=root/'checkout';dest.mkdir();(dest/'README.md').write_text('user code');result=restore(data/'results-index.json',dest);self.assertEqual(result['restored'],3);self.assertEqual((dest/'data/a').read_bytes(),b'new');self.assertEqual((dest/'README.md').read_text(),'user code')
   self.assertEqual(restore(data/'results-index.json',dest)['already_present'],3)
   (dest/'data/a').write_text('user raw data')
   with self.assertRaisesRegex(ValueError,'refusing overwrite'):restore(data/'results-index.json',dest)
 def test_corrupt_layer_rejected_before_any_restore(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);data=root/'datasets';data.mkdir();archive(data/'results-complete.zip',{'data/a':b'old'});initialize(data);archive(root/'new.zip',{'data/a':b'new'});index=publish_snapshot(root/'new.zip',data);(data/index['archives'][-1]['name']).write_bytes(b'corrupted')
   with self.assertRaisesRegex(ValueError,'integrity'):restore(data/'results-index.json',root/'checkout')
   self.assertFalse((root/'checkout/data/a').exists())
 def test_archive_traversal_and_duplicate_paths_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'index.json';good={'format':'codex-results-layers-v1','archives':[{'name':'results-complete.zip','bytes':1,'sha256':'0'*64}],'files':[{'path':'data/a','archive':'results-complete.zip','bytes':1,'sha256':'0'*64}]}
   for bad in ['../outside.zip','/outside.zip']:
    good['archives'][0]['name']=bad;p.write_text(json.dumps(good))
    with self.assertRaises(ValueError):load_index(p)
   good['archives'][0]['name']='results-complete.zip';good['files'][0]['path']='data/../../outside';p.write_text(json.dumps(good))
   with self.assertRaises(ValueError):load_index(p)
 def test_bounded_layers_split_and_failed_update_keeps_valid_index(self):
  from unittest.mock import patch
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);data=root/'datasets';data.mkdir();archive(data/'results-complete.zip',{'data/a':b'a'});initialize(data);archive(root/'new.zip',{'data/a':b'a','data/b':b'bbb','runs/c':b'ccc'})
   with patch('scripts.archive_layers.MAX_LAYER_BYTES',5):index=publish_snapshot(root/'new.zip',data)
   self.assertEqual(len(index['archives']),3);self.assertEqual(restore(data/'results-index.json',root/'checkout')['restored'],3)
   before=(data/'results-index.json').read_bytes();archive(root/'bad.zip',{'data/large':b'123456'})
   with patch('scripts.archive_layers.MAX_LAYER_BYTES',5),self.assertRaises(ValueError):publish_snapshot(root/'bad.zip',data)
   self.assertEqual((data/'results-index.json').read_bytes(),before)
