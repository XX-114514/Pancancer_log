import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import unittest
import numpy as np
from source_local_graph_audit_rules_v1 import validate_partition, validate_model_axes


class Checks(unittest.TestCase):
    def base(self):
        return [np.array([3, 5, 8, 10]), np.array([0, 0, 1, -1]), np.array([1, 1, 0, 0], bool),
                np.array([1, 1, 1, 0], bool), np.array([-1, -1, 2, 2]), 2]

    def test_archive_and_coarse_retention(self):
        self.assertTrue(validate_partition(*self.base()))
        args = self.base(); args[1][-1] = 1
        with self.assertRaises(ValueError): validate_partition(*args)

    def test_duplicate_axes_and_routes(self):
        args = self.base(); args[0][-1] = 3
        with self.assertRaises(ValueError): validate_partition(*args)
        args = self.base(); args[4][1] = 2
        with self.assertRaises(ValueError): validate_partition(*args)

    def test_group_gap_and_fine_domain(self):
        args = self.base(); args[1][2] = 3
        with self.assertRaises(ValueError): validate_partition(*args)
        args = self.base(); args[2][-1] = True
        with self.assertRaises(ValueError): validate_partition(*args)

    def test_model_swap_and_variance(self):
        p = np.array([0, 1, 2]); rows = np.array([3, 5, 8]); g = np.array([0, 0, 1]); c = np.array([5, 9])
        model = dict(source_row_positions=p, native_rows=rows, source_local_group_code=g, native_columns=c,
                     components=np.eye(2), mean=np.zeros(2), fit_embedding=np.array([[1., 0.], [-1., 0.], [0., 2.]]),
                     explained_variance=np.array([1., 2.]))
        self.assertTrue(validate_model_axes(model, p, rows, g, c))
        with self.assertRaises(ValueError): validate_model_axes(model, p, rows[::-1], g, c)
        model['explained_variance'] = np.array([1., 3.])
        with self.assertRaises(AssertionError): validate_model_axes(model, p, rows, g, c)


if __name__ == '__main__': unittest.main()
