from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_supplemental_native_context_v1 import anatomy_values, context_counts


class NativeContext(unittest.TestCase):
    def test_cancer_origin_or_mixed_free_text_not_assigned_biopsy_site(self):
        actual = anatomy_values(['BRCA', 'breast liver metastasis', 'liver_subject_key', ' Liver '])
        self.assertEqual(actual.tolist(), ['unreviewed_value'] * 3 + ['liver'])

    def test_unknown_missing_and_known_blood_remain_separate(self):
        actual = anatomy_values(['unknown', '', 'unreviewed_term', 'peripheral blood'])
        self.assertEqual(actual.tolist(), ['missing_or_explicit_unknown'] * 2 + ['unreviewed_value', 'blood'])

    def test_multiple_cues_not_counted_as_exclusive_identity(self):
        result = context_counts(np.array(['liver', 'liver', 'blood']),
            np.array([[1, 1], [0, 1], [0, 0]], bool), ['cue_a', 'cue_b'])
        liver = next(r for r in result if r['context'] == 'liver')
        self.assertEqual(liver['n_archive_rows'], 2)
        self.assertEqual(liver['cue_positive_counts'], {'cue_a': 1, 'cue_b': 2})


if __name__ == '__main__':
    unittest.main()
