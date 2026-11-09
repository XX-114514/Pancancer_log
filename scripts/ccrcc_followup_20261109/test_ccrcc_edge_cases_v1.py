"""Nontrivial leakage, multiplicity and geometry counterexamples; no source data."""
import unittest
import numpy as np
import shapely
from shapely import STRtree
from statsmodels.stats.multitest import multipletests
from ccrcc_followup_v1 import identities, fdr


class EdgeCases(unittest.TestCase):
    def test_exposure_cannot_supply_missing_EC_identity_and_conflict_excluded(self):
        a = {'PECAM1': np.array([1, 1, 0, 1]), 'VWF': np.array([0, 1, 0, 1]),
             'ACKR1': np.array([100, 0, 100, 100]),
             'CD3D': np.array([0, 0, 0, 1]), 'CD3E': np.array([0, 0, 0, 1])}
        result = identities(a, np.repeat(100, 4), np.repeat(30, 4))
        self.assertEqual(list(result[-1]), ['Other', 'Blood_EC', 'Other', 'Other'])
        self.assertEqual(list(result[3]), [False, False, False, True])
        a['ACKR1'] = np.zeros(4)
        changed = identities(a, np.repeat(100, 4), np.repeat(30, 4))
        np.testing.assert_array_equal(result[-1], changed[-1])
        np.testing.assert_array_equal(result[5], changed[5])

    def test_full_family_BH_BY_with_ties_and_missing_slots_as_one(self):
        p = np.array([.001, .1, 1, .001, .7] * 20)[:96]
        for by, method in [(False, 'fdr_bh'), (True, 'fdr_by')]:
            np.testing.assert_allclose(fdr(p, by), multipletests(p, method=method)[1], atol=1e-14)

    def test_nonunique_and_interior_polygon_queries_are_not_arbitrarily_assigned(self):
        polygons = np.array([shapely.box(0, 0, 1, 1), shapely.box(2, 0, 3, 1)])
        points = shapely.points([[1.5, .5], [.5, .5], [4, .5]])
        pairs, distances = STRtree(polygons).query_nearest(points, all_matches=True, return_distance=True)
        hits = np.bincount(pairs[0], minlength=3)
        accepted = (hits[pairs[0]] == 1) & (distances > 1e-8)
        self.assertEqual(list(pairs[0, accepted]), [2])
        self.assertEqual(list(distances[accepted]), [1.])


if __name__ == '__main__':
    unittest.main(verbosity=2)
