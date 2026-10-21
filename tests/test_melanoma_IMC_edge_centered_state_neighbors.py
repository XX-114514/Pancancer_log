"""Reuse full workflow fixtures plus an independent marginal-affinity case."""
import importlib.util
import itertools
from pathlib import Path
import unittest
import numpy as np

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('raw_organization_fixtures',
    root/'tests'/'test_melanoma_IMC_continuous_state_neighbors.py')
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)


class EdgeCenteredOrganization(fixtures.StateOrganization):
    @classmethod
    def setUpClass(cls):
        path = root/'scripts'/'analyze_melanoma_IMC_continuous_state_neighbors_v2.py'
        spec = importlib.util.spec_from_file_location('edge_centered_under_test', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def test_analytic_null_equals_complete_enumeration(self):
        left = np.array([[-2., 3.], [2., 5.]])
        right = np.array([[1., -3.], [5., 1.]])
        a, b = np.array([0, 0, 1]), np.array([0, 1, 1])
        all_statistics = []
        for p, q in itertools.product(itertools.permutations(range(2)), repeat=2):
            l = left[np.array(p)[a]]
            r = right[np.array(q)[b]]
            independent = sum(np.outer(x-l.mean(0), y-r.mean(0)) for x, y in zip(l, r))/len(a)
            all_statistics.append(independent)
        observed, expected, _ = self.m.roi_permutations(left, right, a, b,
            np.zeros(2), np.zeros(2), 7, np.random.default_rng(5))
        np.testing.assert_allclose(expected, np.mean(all_statistics, axis=0), atol=1e-12)
        l, r = left[a], right[b]
        independent = sum(np.outer(x-l.mean(0), y-r.mean(0)) for x, y in zip(l, r))/len(a)
        np.testing.assert_allclose(observed, independent, atol=1e-12)

    def test_separate_interface_preferences_are_not_pair_covariance(self):
        values, _ = self.m.rank_z(np.arange(20)[:, None])
        edges = np.array(list(itertools.product(range(15, 20), repeat=2)))
        a, b = edges[:, 0], edges[:, 1]
        raw_observed = float((values[a].T @ values[b]/len(a))[0, 0])
        self.assertGreater(raw_observed, 1.)
        observed, expected, null = self.m.roi_permutations(values, values, a, b,
            np.zeros(20), np.zeros(20), 2000, np.random.default_rng(3))
        # Every chosen endpoint is independently paired with every other one:
        # high-state interface preference exists but no pair covariance exists.
        np.testing.assert_allclose(observed, 0, atol=1e-12)
        np.testing.assert_allclose(expected, 0, atol=1e-12)
        np.testing.assert_allclose(null, 0, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
