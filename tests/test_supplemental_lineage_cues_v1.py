from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from native_supplemental_lineage_evidence_v1 import GENES, cue_support


class SupplementaryCues(unittest.TestCase):
    def row(self, *genes, missing=()):
        measured = np.array([g not in missing for g in GENES])
        detected = np.array([[g in genes for g in GENES]])
        return cue_support(detected, measured)[0]

    def test_secreted_liver_pair_alone_is_not_supported(self):
        self.assertFalse(self.row('ALB', 'TTR')[0])
        self.assertTrue(self.row('ALB', 'ASGR1')[0])

    def test_generic_AXL_NGFR_is_not_pigment_lineage(self):
        self.assertFalse(self.row('AXL', 'NGFR', 'MITF')[1])
        self.assertTrue(self.row('PMEL', 'MLANA')[1])

    def test_secretory_neural_signal_needs_transcription_factor(self):
        self.assertFalse(self.row('CHGA', 'CHGB', 'SYP')[2])
        self.assertTrue(self.row('CHGA', 'INSM1')[2])

    def test_missing_anchor_not_counted_as_detected(self):
        self.assertFalse(self.row('ALB', 'ASGR1', missing=('ASGR1',))[0])
        self.assertFalse(self.row('CHGA', 'INSM1', missing=('INSM1',))[2])


if __name__ == '__main__':
    unittest.main()
