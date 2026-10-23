"""Counterexamples for unit balance, rare strata and deterministic sampling."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import numpy as np
    from annotation_pilot_sampling_v1 import balanced_allocations, balanced_cell_selection, stable_seed
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'NumPy unavailable')
class SamplingTests(unittest.TestCase):
    def test_waterfill_preserves_exact_budget_without_oversampling(self):
        values = balanced_allocations([1, 100, 100], 9, np.random.default_rng(1))
        np.testing.assert_array_equal(values, [1, 4, 4])
        np.testing.assert_array_equal(balanced_allocations([1, 2], 100, np.random.default_rng(1)), [1, 2])
        for total in range(12):
            x = balanced_allocations([0, 2, 4, 5], total, np.random.default_rng(2))
            self.assertEqual(x.sum(), min(total, 11))
            self.assertTrue((x <= [0, 2, 4, 5]).all())

    def test_deterministic_eligibility_unit_balance_and_rare_reserve(self):
        units = np.array(['small'] * 2 + ['large'] * 100 + ['excluded'] * 3)
        fine = np.array(['rare', 'common'] + ['rare'] + ['common'] * 99 + ['other'] * 3)
        eligible = np.arange(105) < 102
        seed = stable_seed(20261023, 'source', 'T_NK')
        a = balanced_cell_selection(units, fine, eligible, 12, seed)
        b = balanced_cell_selection(units, fine, eligible, 12, seed)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(len(a), 12)
        self.assertEqual(len(set(a)), 12)
        self.assertTrue(set([0, 1, 2]).issubset(set(a)))
        self.assertFalse((a >= 102).any())

    def test_small_budget_and_empty_pool_are_not_fabricated_cells(self):
        a = balanced_cell_selection(['a', 'b', 'c'], ['x', 'y', 'z'], [True] * 3, 1, 4)
        self.assertEqual(len(a), 1)
        self.assertEqual(len(balanced_cell_selection(['a'], ['x'], [False], 10, 4)), 0)
        all_rows = balanced_cell_selection(['a', 'b'], ['x', 'x'], [True, True], 100, 4)
        np.testing.assert_array_equal(all_rows, [0, 1])

    def test_missing_units_and_mismatched_axes_fail(self):
        with self.assertRaises(ValueError):
            balanced_cell_selection([''], ['x'], [True], 1, 4)
        with self.assertRaises(ValueError):
            balanced_cell_selection(['a'], ['x', 'y'], [True], 1, 4)
        self.assertNotEqual(stable_seed(1, 's', 'a'), stable_seed(1, 's', 'b'))


if __name__ == '__main__':
    unittest.main()
