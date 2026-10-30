"""Actual standard API/domain and conditional-coordinate counterexamples."""
import sys
import unittest
from pathlib import Path
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from standard_annotation_order_primitives_v1 import (
    select_features, pca_fit, conditional_readout_metrics, coarse_routes_from_local,
    technical, GATE_ORDER)


class StandardOrderTests(unittest.TestCase):
    def test_official_count_HVG_and_logged_input_refusal(self):
        rng = np.random.default_rng(2026103001)
        x = sp.csr_matrix(rng.poisson(np.geomspace(.2, 8, 400), size=(160, 400)).astype(np.float32))
        genes = ['HBB', 'HBEGF'] + ['G'+str(i) for i in range(398)]
        source = np.repeat(['a', 'b'], 80)
        selected, table = select_features(x, source, 50, genes, 100)
        self.assertEqual(len(selected), 50)
        self.assertEqual(int(table.selected_HVG.sum()), 50)
        self.assertNotIn(0, selected)
        self.assertFalse(technical('HBEGF'))
        logged = x.copy()
        logged.data = np.log1p(logged.data)
        with self.assertRaises(ValueError):
            select_features(logged, source, 50, genes, 100)

    def test_standard_centered_PCA_projection(self):
        rng = np.random.default_rng(2026103002)
        x = sp.csr_matrix(np.log1p(rng.poisson(2, size=(100, 60))).astype(np.float32))
        z, model = pca_fit(x, 42, 5, 3)
        expected = (x.toarray() - model.mean_) @ model.components_.T
        np.testing.assert_allclose(z, expected, rtol=3e-4, atol=3e-3)
        np.testing.assert_allclose(model.components_ @ model.components_.T, np.eye(5), atol=3e-4)

    def test_official_Morans_I_matches_direct_definition_and_missing_stays_missing(self):
        n, k = 60, 4
        t = np.arange(n, dtype=float)
        y = np.sin(t / 7)
        z = np.column_stack([t, np.zeros(n)])
        result = conditional_readout_metrics(z, {'signal': y}, np.full(n, 1000),
            np.full(n, np.nan), np.arange(n), k)
        W = np.zeros((n, n))
        for i in range(n):
            order = np.argsort(np.abs(t - t[i]), kind='stable')
            W[i, order[order != i][:k]] = 1
        centered = y - y.mean()
        direct = n * np.sum(W * centered[:, None] * centered[None, :]) / (W.sum() * (centered**2).sum())
        self.assertAlmostEqual(result['states']['signal']['depth_adjusted_official_Scanpy_Morans_I'], direct, places=10)
        self.assertEqual(result['gamma_delta_RNA_candidate']['status'], 'unmeasured_not_negative')

    def test_independent_coordinate_spaces_are_not_concatenated_for_neighbors(self):
        n = 60
        with self.assertRaises(ValueError):
            conditional_readout_metrics(np.arange(n)[:, None], {}, np.ones(n),
                np.zeros(n), np.arange(n), 4, space_ids=np.tile(['a', 'b'], 30))

    def test_source_local_RNA_routing_keeps_mixed_group_unresolved(self):
        flags = np.zeros((120, len(GATE_ORDER)), bool)
        flags[:40, GATE_ORDER.index('T')] = True
        flags[40:80, GATE_ORDER.index('B')] = True
        flags[80:, GATE_ORDER.index('T')] = True
        flags[80:, GATE_ORDER.index('B')] = True
        routes, rows = coarse_routes_from_local(np.repeat('source', 120),
            np.repeat([0, 1, 2], 40), flags, {'source': dict.fromkeys(GATE_ORDER, True)})
        self.assertTrue((routes[:40] == 'T_NK').all())
        self.assertTrue((routes[40:80] == 'B').all())
        self.assertTrue((routes[80:] == 'unresolved_context').all())
        self.assertEqual(sum(r['n_cells'] for r in rows), 120)


if __name__ == '__main__':
    unittest.main()
