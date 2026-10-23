"""Independent score replay and changed-model-axis counterexamples."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import numpy as np
    from audit_annotation_source_holdout_v1 import prototype_prediction, independent_scores, restore_reference, check_score_bounds
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'Run in actual scientific environment')
class HoldoutAuditTests(unittest.TestCase):
    def test_inverse_distance_exact_match_and_unequal_weight_prediction(self):
        centers = np.array([[0.], [2.]])
        outcomes = np.array([[1., 4.], [3., 2.]])
        np.testing.assert_allclose(prototype_prediction(centers, outcomes, np.array([[0.]]), 1), [[1., 4.]])
        np.testing.assert_allclose(prototype_prediction(centers, outcomes, np.array([[.5]]), 2), [[1.5, 3.5]])
        with self.assertRaises(ValueError):prototype_prediction(centers, outcomes[:1], np.array([[0.]]), 1)

    def test_depth_residual_scoring_and_missing_denominator(self):
        n = 120;depth = np.arange(1, n + 1, dtype=float)
        y = np.sin(depth) + .7 * np.log1p(depth)
        result = independent_scores(y, y + 2., depth, 100)
        self.assertAlmostEqual(result['MAE_native_logCP10K'], 2.)
        self.assertAlmostEqual(result['depth_adjusted_Pearson'], 1.)
        self.assertAlmostEqual(result['depth_adjusted_Spearman'], 1.)
        missing = independent_scores(y, np.full(n, np.nan), depth, 100)
        self.assertEqual(missing, {'status': 'insufficient_evaluable_query_cells', 'n_cells': 0})
        with self.assertRaises(ValueError):check_score_bounds({'n_cells': 121, 'n_query_route_cells': 120, 'status': 'estimated'})

    def test_minimal_reference_gene_order_or_index_mutation_fails(self):
        path = Path(tempfile.mkdtemp(prefix='heldout-audit-counterexample-')) / 'reference.npz'
        np.savez_compressed(path, genes=np.array(['a', 'b', 'c']), mean=np.zeros(3), std=np.ones(3),
                            loadings=np.eye(3), C=np.ones((2, 3)), Nr=np.ones(2) * 10,
                            sigma=np.ones(2) * .1, feature_candidate_indices=np.array([4, 8, 10]),
                            metadata_json=np.asarray(json.dumps({'synthetic': True})))
        ref, meta = restore_reference(path, ['a', 'b', 'c'], [4, 8, 10])
        self.assertEqual(ref.shape, (0, 3));self.assertTrue(meta['synthetic'])
        with self.assertRaises(ValueError):restore_reference(path, ['b', 'a', 'c'], [4, 8, 10])
        with self.assertRaises(ValueError):restore_reference(path, ['a', 'b', 'c'], [4, 8, 11])

    def test_transposed_nonsquare_Harmony_centers_rejected(self):
        path = Path(tempfile.mkdtemp(prefix='heldout-center-axis-counterexample-')) / 'reference.npz'
        np.savez_compressed(path, genes=np.array(['a', 'b', 'c']), mean=np.zeros(3), std=np.ones(3),
                            loadings=np.eye(3), C=np.ones((3, 2)), Nr=np.ones(2) * 10,
                            sigma=np.ones(2) * .1, feature_candidate_indices=np.array([4, 8, 10]),
                            metadata_json=np.asarray('{}'))
        with self.assertRaises(ValueError):restore_reference(path, ['a', 'b', 'c'], [4, 8, 10])


if __name__ == '__main__':unittest.main()
