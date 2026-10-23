import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import numpy as np
from representative_local_mapping_v1 import weighted_vote, without_self, support_mask


class MappingTests(unittest.TestCase):
    def test_exact_matches_exclude_nonzero_neighbors(self):
        labels, votes = weighted_vote([[0, 0, .1, .1, .1]], [[4, 4, 1, 1, 1]])
        np.testing.assert_array_equal(labels, [4])
        np.testing.assert_array_equal(votes, [1])

    def test_distance_weight_and_exact_tie(self):
        labels, votes = weighted_vote([[.1, 1, 1], [1, 1, 1]], [[5, 2, 2], [8, 3, 5]])
        np.testing.assert_array_equal(labels, [5, 3])
        self.assertGreater(votes[0], .8)
        self.assertAlmostEqual(votes[1], 1 / 3)

    def test_remove_self_by_index_not_zero_distance(self):
        out = without_self([[0, 0, 2], [0, 0, 3]], [[7, 9, 8], [7, 9, 8]], [9, 7], 2)
        np.testing.assert_array_equal(out, [[0, 2], [0, 3]])

    def test_support_retains_separate_distance_and_vote_requirements(self):
        np.testing.assert_array_equal(support_mask([1, 2, 1], [.9, .9, .4], 1, .6), [True, False, False])

    def test_reject_bad_axes_and_invalid_values(self):
        for d, labels in [([[1]], [[-1]]), ([[np.nan]], [[1]]), ([[1]], [[1, 2]])]:
            with self.assertRaises(ValueError): weighted_vote(d, labels)
        with self.assertRaises(ValueError): without_self([[0]], [[0]], [0], 1)
        with self.assertRaises(ValueError): support_mask([1], [1.1], 1, .6)

    def test_empty_query_preserves_typed_empty_results(self):
        labels, votes = weighted_vote(np.empty((0, 3)), np.empty((0, 3), np.int32))
        self.assertEqual(labels.shape, (0,))
        self.assertEqual(votes.shape, (0,))


if __name__ == '__main__': unittest.main()
