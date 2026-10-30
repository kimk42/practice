"""Recovery must reuse the reviewed folder and reject stale charge tables."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import finalize_pc_checkpointed as finalizer
from pc_bcc_from_equivalence import validate_equivalence_input
from pc_length_parameters import ROOT, sha


class RecoveryInputTests(unittest.TestCase):
    def test_nested_candidate_returned_without_reprocessing(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name); qm=root/'qm'; (qm/'pc8').mkdir(parents=True)
            raw=root/'accepted.out'; raw.write_text('accepted QM fixture')
            (qm/'pc8/completed.json').write_text(json.dumps({
                'final_output':'accepted.out','output_sha256':sha(raw)}))
            base=root/'forcefield'; out=base/'pc8'; nested=out/'recovered/PC8.amb2gmx'
            nested.mkdir(parents=True); (nested/'PC8_GMX.itp').write_text('topology fixture')
            (out/'parameter-audit.json').write_text(json.dumps({
                'topology_folder':str(nested.relative_to(root)),
                'outputs_sha256':{'PC8_GMX.itp':sha(nested/'PC8_GMX.itp')}}))
            with patch.object(finalizer,'ROOT',root), patch.object(finalizer,'BASE',base), \
                 patch.object(finalizer,'route',return_value=(qm,'native')), \
                 patch.object(finalizer,'parse_completed',return_value=(None,{'gradient_gate_passed':True})), \
                 patch.object(finalizer,'run_group',side_effect=AssertionError('must not rerun')):
                self.assertEqual(finalizer.finalize(8),nested)
                (nested/'PC8_GMX.itp').write_text('changed')
                with self.assertRaisesRegex(ValueError,'forcefield changed'):
                    finalizer.finalize(8)

    def test_saved_pc8_charges_are_bound_to_original_qm(self):
        record=json.loads((ROOT/'models/surface-revision-v2/equatom-qsort-validation-20261030.json').read_text())['cases']['8']
        source=ROOT/'runs/pc-length-forcefield-v2/pc8'
        validate_equivalence_input(8,record,source/'accepted-sqm.out',source/'ANTECHAMBER_AC.AC')
        changed=copy.deepcopy(record)
        changed['averaged_charges_e'][0]+=.01
        changed['averaged_charges_e'][1]-=.01
        with self.assertRaisesRegex(ValueError,'accepted QM charges'):
            validate_equivalence_input(8,changed,source/'accepted-sqm.out',source/'ANTECHAMBER_AC.AC')


if __name__=='__main__':
    unittest.main()
