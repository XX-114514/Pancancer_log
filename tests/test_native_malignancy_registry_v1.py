from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_native_malignancy_evidence_registry_v1 import pair_counts


class NativeMalignancyFlags(unittest.TestCase):
    def test_false_result_without_evaluation_is_not_negative(self):
        self.assertEqual(pair_counts(['False', 'unknown'], ['False', 'False']),
            {'explicitly_not_evaluable_not_negative': 1, 'evaluation_unknown_result_uninterpreted': 1})

    def test_true_evaluation_separates_positive_negative_and_missing(self):
        self.assertEqual(pair_counts(['True'] * 3, ['True', 'False', 'nan']),
            {'evaluable_positive_flag': 1, 'evaluable_negative_flag': 1, 'evaluable_result_unknown': 1})

    def test_conflict_and_unknown_encoding_are_held(self):
        self.assertEqual(pair_counts(['False', 'mystery'], ['True', 'False']),
            {'conflict_positive_flag_while_not_evaluable': 1, 'unrecognized_boolean_encoding_requires_review': 1})


if __name__ == '__main__':
    unittest.main()
