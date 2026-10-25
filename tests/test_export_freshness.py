import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from scripts.finalize_workbench import export

class ExportFreshnessTests(unittest.TestCase):
 def test_failed_job_does_not_inherit_previous_geometry_or_solvent(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);job=root/'runs/workbench'/('a'*32);job.mkdir(parents=True);(job/'config.json').write_text(json.dumps({'protocol_id':'continue-npt'}));(job/'state.json').write_text(json.dumps({'status':'failed'}));(job/'report.html').write_text('<body>failed</body>')
   for p in [root/'artifacts/workbench/preparation',root/'artifacts/pe-continuation']:
    (p/'solvent-audit').mkdir(parents=True);(p/'solvent-audit/solvent-audit.json').write_text('old source');(p/'geometry-analysis.json').write_text('old geometry')
   with patch('scripts.finalize_workbench.ROOT',root),patch('workbench.report.publish'):export(job)
   for p in [root/'artifacts/workbench/preparation',root/'artifacts/pe-continuation']:
    self.assertFalse((p/'solvent-audit/solvent-audit.json').exists());self.assertFalse((p/'geometry-analysis.json').exists())
   self.assertEqual((root/'artifacts/pe-continuation/report.html').read_text(),'<body>failed</body>')
