"""Counterexamples for fold leakage, feature scope and paired missingness."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import numpy as np
import scipy.sparse as sp
from lineage_holdout_diagnostics_v1 import fold_members, training_feature_sets, paired_profiles, VARIANTS
from annotation_representation_pilot_v1 import select_features


class HoldoutTests(unittest.TestCase):
    def test_query_perturbation_cannot_change_training_features(self):
        rng = np.random.default_rng(14)
        x = sp.csr_matrix(rng.poisson(2, (90, 18)).astype(float))
        source = np.repeat(['a', 'b', 'q'], 30)
        route = np.tile(np.repeat(['T', 'M'], 15), 3)
        train, query = fold_members(source, source, np.ones(90, bool), 'q')
        first = training_feature_sets(x, source, route, train, 'T', 7)
        changed = x.toarray(); changed[query] = rng.poisson(100, (30, 18))
        second = training_feature_sets(sp.csr_matrix(changed), source, route, train, 'T', 7)
        for key in VARIANTS:
            np.testing.assert_array_equal(first[key], second[key])
        local = train[route[train] == 'T']
        np.testing.assert_array_equal(first[VARIANTS[1]], select_features(x[local], source[local], 7))

    def test_possible_overlap_sources_share_query_fold(self):
        source = np.array(['a', 'b', 'c', 'd'])
        fold = np.array(['a', 'paired', 'paired', 'd'])
        train, query = fold_members(source, fold, np.ones(4, bool), 'paired')
        self.assertEqual(query.tolist(), [1, 2]); self.assertEqual(train.tolist(), [0, 3])

    def test_split_source_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Source split'):
            fold_members(['a', 'a', 'b'], ['q', 't', 't'], [True] * 3, 'q')

    def test_zero_RNA_not_fitted_or_scored(self):
        train, query = fold_members(['a', 'a', 'b', 'b'], ['a', 'a', 'b', 'b'], [1, 0, 0, 1], 'b')
        self.assertEqual(train.tolist(), [0]); self.assertEqual(query.tolist(), [3])

    def test_paired_missing_predictions_share_denominator(self):
        n = 8; truth = np.arange(n, dtype=float)[:, None]
        pred = {v: truth.copy() for v in VARIANTS}; pred[VARIANTS[1]][0] = np.nan
        support = {v: np.ones(n, bool) for v in VARIANTS}; support[VARIANTS[0]][1] = False
        rows = paired_profiles(np.repeat('s', n), np.repeat('unresolved', n), truth, np.ones(n),
            {'all_route_cells_not_pure_type': pred}, support, ['readout'], 3)
        self.assertEqual(rows[0]['n_common_query_cells'], 6)
        self.assertEqual({m['n_cells'] for m in rows[0]['metrics'].values()}, {6})
        self.assertAlmostEqual(rows[0]['delta_lineage_minus_shared_depth_adjusted_Pearson'], 0)

    def test_missing_direct_candidate_is_not_a_fine_identity(self):
        n = 6; truth = np.arange(n, dtype=float)[:, None]
        rows = paired_profiles(np.repeat('s', n), np.repeat('unresolved', n), truth, np.ones(n),
            {'all_route_cells_not_pure_type': {v: truth for v in VARIANTS}},
            {v: np.ones(n, bool) for v in VARIANTS}, ['readout'], 10)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['status'], 'insufficient_coverage')


if __name__ == '__main__':
    unittest.main()
