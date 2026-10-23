import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import numpy as np
from native_local_fit_pool_v1 import boolean_evidence, fit_pool


class FitPoolTests(unittest.TestCase):
    def test_missing_does_not_mean_assessed_negative(self):
        yes, known = boolean_evidence(['True', 'False', '__SOURCE_VALUE_MISSING__', 'unreviewed'], 4)
        np.testing.assert_array_equal(yes, [True, False, False, False])
        np.testing.assert_array_equal(known, [True, True, False, False])

    def test_archive_records_retained_separate_from_fit_pool(self):
        r = fit_pool([1, 0, 1, 1], [1, 1, 0, 1], ['lib1'] * 4, [False, False, False, True])
        np.testing.assert_array_equal(r['fit_pool_candidate'], [True, False, False, False])
        self.assertTrue(all(len(v) == 4 for v in r.values()))

    def test_missing_unit_and_explicit_exclusion(self):
        r = fit_pool([1] * 4, [1] * 4, ['lib1', '', 'nan', 'lib2'], excluded=['true', 'false', 'false', 'unknown'])
        np.testing.assert_array_equal(r['fit_pool_candidate'], [False, False, False, True])
        self.assertFalse(r['carried_primary_exclusion_assessed_value'][-1])

    def test_reject_shape_errors(self):
        with self.assertRaises(ValueError): fit_pool([1, 1], [1], ['a', 'b'])
        with self.assertRaises(ValueError): boolean_evidence([True], 2)


if __name__ == '__main__': unittest.main()
