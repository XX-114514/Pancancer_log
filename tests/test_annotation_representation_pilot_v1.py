"""Counterexamples and real scientific-library API qualification for pilot."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import numpy as np
    import scipy.sparse as sp
    from annotation_representation_pilot_v1 import (normalized_log_counts,
        feature_columns, select_features, pca_fit, harmonize, cluster,
        conditional_readout_metrics)
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'Run in qualified scientific environment')
class RepresentationTests(unittest.TestCase):
    def test_full_native_denominator_zero_and_invalid_count(self):
        x, ok = normalized_log_counts(sp.csr_matrix([[2, 0], [0, 0]]), [10, 0])
        self.assertAlmostEqual(float(x[0, 0]), float(np.log1p(2000)), places=5)
        np.testing.assert_array_equal(ok, [True, False])
        for x, total in [([[2.1]], [4]), ([[4]], [3]), ([[-1]], [4])]:
            with self.assertRaises(ValueError):normalized_log_counts(sp.csr_matrix(x), total)

    def test_physical_readout_alias_leakage_rejected(self):
        lookup = {'MODEL': {0}, 'HELDOUT_ALIAS': {0}, 'M2': {1}, 'OTHER': {2}}
        with self.assertRaises(ValueError):feature_columns(lookup, ['MODEL', 'M2'], ['HELDOUT_ALIAS'])
        np.testing.assert_array_equal(feature_columns(lookup, ['MODEL', 'M2'], ['OTHER']), [0, 1])

    def test_feature_selection_takes_training_only_and_ties_stable(self):
        x = sp.csr_matrix(np.tile([[0, 1, 0], [10, 1, 0]], (10, 1)), dtype=np.float32)
        chosen = select_features(x, np.array(['TRAIN'] * 20), 1)
        np.testing.assert_array_equal(chosen, [0])
        with self.assertRaises(ValueError):select_features(x, np.array(['TRAIN']), 1)

    def test_conditional_gradient_and_constant_signal(self):
        n = 100
        z = np.arange(n, dtype=np.float32)[:, None]
        out = conditional_readout_metrics(z, {'gradient': z.ravel(), 'constant': np.ones(n)},
                                         np.ones(n) * 10, np.zeros(n, bool), np.arange(n), 5)
        self.assertGreater(out['states']['gradient']['depth_adjusted_neighbor_Pearson'], .99)
        self.assertEqual(out['states']['constant']['status'], 'no_variable_signal')
        self.assertIsNone(out['gamma_delta_RNA_candidate']['conditional_neighbor_candidate_fraction'])

    def test_actual_PCA_Harmony_and_graph_APIs_no_accuracy_claim(self):
        rng = np.random.default_rng(5)
        x = sp.csr_matrix(rng.poisson(2, (160, 60)), dtype=np.float32)
        source = np.array(['A'] * 80 + ['B'] * 80)
        raw, model = pca_fit(x, 1, 10, 2)
        cfg = dict(theta=1., **{'lambda': 1.}, max_iter_harmony=2, threads=4,
                   neighbors=10, resolution=1.)
        corrected, info = harmonize(raw, source, 1, cfg)
        self.assertEqual(corrected.shape, (160, 10))
        labels, count = cluster(corrected, 1, cfg)
        self.assertEqual(len(labels), 160)
        self.assertGreater(count, 0)
        self.assertTrue(np.isfinite(model.components_).all())
        self.assertEqual(info['status'], 'discovery_all_sources_jointly_fit')


if __name__ == '__main__':unittest.main(verbosity=2)
