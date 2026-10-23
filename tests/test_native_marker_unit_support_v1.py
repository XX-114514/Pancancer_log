from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import unittest
import numpy as np
from native_marker_unit_support_v1 import summarize_unit_direction


class Checks(unittest.TestCase):
    def test_missing_and_zero_are_different(self):
        r=summarize_unit_direction([20,20,20],[0,40,40],[20,0,19],[0,0,19])
        self.assertEqual(r['n_comparable_expression_units'],1)
        self.assertEqual(r['median_log2_mean_ratio'],0.)
        self.assertEqual(r['n_units_positive_mean_ratio'],0)

    def test_direction_and_own_group_excluded_comparator(self):
        r=summarize_unit_direction([20,20],[40,20],[20,20],[20,40])
        self.assertEqual(r['n_units_positive_mean_ratio'],1)
        self.assertEqual(r['n_units_negative_mean_ratio'],1)
        self.assertAlmostEqual(r['median_log2_mean_ratio'],0.)
        r=summarize_unit_direction([5],[10],[0],[0])
        self.assertIsNone(r['median_log2_mean_ratio'])

    def test_invalid_axes_and_sums(self):
        for args in [([20],[2,3],[20],[1]),([20],[-1],[20],[1]),([20],[np.nan],[20],[1])]:
            with self.assertRaises(ValueError): summarize_unit_direction(*args)


if __name__ == '__main__':unittest.main()
