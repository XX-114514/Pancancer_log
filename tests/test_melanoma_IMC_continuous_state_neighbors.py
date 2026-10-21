"""Synthetic mathematical and full-wiring checks; no real source records."""
import importlib.util
import itertools
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd


class StateOrganization(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1]/'scripts'/'analyze_melanoma_IMC_continuous_state_neighbors_v1.py'
        spec = importlib.util.spec_from_file_location('organization_under_test', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def test_rank_ties_and_constant_column(self):
        scores, valid = self.m.rank_z([[0, 1], [0, 1], [2, 1], [4, 1]])
        self.assertEqual(valid.tolist(), [True, False])
        self.assertAlmostEqual(scores[:, 0].mean(), 0)
        self.assertAlmostEqual(np.mean(scores[:, 0]**2), 1)
        self.assertEqual(scores[0, 0], scores[1, 0])
        self.assertTrue((scores[:, 1] == 0).all())

    def test_whole_marker_vectors_and_strata_preserved(self):
        values = np.column_stack([np.arange(12), 3*np.arange(12)+2])
        groups = np.repeat([0, 1, 2], 4)
        shuffled = self.m.shuffled(values, groups, np.random.default_rng(14))
        for group in range(3):
            take = groups == group
            self.assertEqual(sorted(shuffled[take, 0]), sorted(values[take, 0]))
        np.testing.assert_array_equal(shuffled[:, 1], shuffled[:, 0]*3+2)

    def test_analytic_null_equals_complete_enumeration(self):
        left = np.array([[-2., 3.], [2., 5.]])
        right = np.array([[1., -3.], [5., 1.]])
        a, b = np.array([0, 0, 1]), np.array([0, 1, 1])
        all_statistics = []
        for p, q in itertools.product(itertools.permutations(range(2)), repeat=2):
            statistic = np.zeros((2, 2))
            for i, j in zip(a, b):
                statistic += np.outer(left[p[i]], right[q[j]])/len(a)
            all_statistics.append(statistic)
        observed, expected, _ = self.m.roi_permutations(left, right, a, b,
            np.zeros(2), np.zeros(2), 7, np.random.default_rng(5))
        np.testing.assert_allclose(expected, np.mean(all_statistics, axis=0))
        manual = sum(np.outer(left[i], right[j]) for i, j in zip(a, b))/len(a)
        np.testing.assert_allclose(observed, manual)

    def test_regional_gradient_removed_by_stratified_null(self):
        values = np.array([[-1.], [-1.], [1.], [1.]])
        a = b = np.arange(4)
        group = np.array([0, 0, 1, 1])
        raw, expected_raw, _ = self.m.roi_permutations(values, values, a, b,
            np.zeros(4), np.zeros(4), 30, np.random.default_rng(1))
        observed, expected, null = self.m.roi_permutations(values, values, a, b,
            group, group, 30, np.random.default_rng(1))
        self.assertEqual(float((raw-expected_raw)[0, 0]), 1.)
        self.assertEqual(float((observed-expected)[0, 0]), 0.)
        np.testing.assert_allclose(null, np.broadcast_to(expected, null.shape))

    def test_known_paired_arrangement_and_permutation_resolution(self):
        values, _ = self.m.rank_z(np.arange(10)[:, None])
        a = b = np.arange(10)
        observed, expected, null = self.m.roi_permutations(values, values, a, b,
            np.zeros(10), np.zeros(10), 2000, np.random.default_rng(71))
        self.assertAlmostEqual(float(observed[0, 0]), 1.)
        self.assertAlmostEqual(float(expected[0, 0]), 0.)
        p = (1+np.count_nonzero(np.abs(null[:, 0, 0]) >= 1-1e-12))/2001
        self.assertGreaterEqual(p, 1/2001)
        self.assertLess(p, .01)

    def test_patient_ROI_aggregation_and_missing_family_slots(self):
        n = 10
        observed = np.ones((n, 2))
        observed[0, 0] = 20
        expected = np.zeros_like(observed)
        draws = np.zeros((n, 2000, 2))
        roi_counts = np.zeros((n, 2), dtype=int)
        roi_counts[:, 0] = 1
        roi_counts[0, 0] = 20
        rows = self.m.finish_family(observed, expected, draws, roi_counts,
            [dict(name='eligible'), dict(name='missing')], 10)
        self.assertEqual(rows[0]['source_case_units'], 10)
        self.assertEqual(rows[0]['contributing_ROIs'], 29)
        self.assertEqual(rows[0]['effect'], 1.)
        self.assertEqual(rows[0]['case_bootstrap_95_CI'], [1., 1.])
        self.assertAlmostEqual(rows[0]['BH'], 2/2001)
        self.assertIsNone(rows[1]['P'])
        self.assertIsNone(rows[1]['effect'])
        self.assertEqual(rows[1]['status'], 'INSUFFICIENT_SOURCE_CASES')

    def test_common_grid_translation_scale_and_constant_axis(self):
        coords = np.array([[0., 0.], [1., 2.], [2., 4.], [4., 6.]])
        groups = self.m.tile_groups(coords, (coords.min(0), coords.max(0)))
        changed = coords*np.array([2., 3.])+np.array([8., -7.])
        np.testing.assert_array_equal(groups, self.m.tile_groups(changed, (changed.min(0), changed.max(0))))
        constant = np.array([[1., 0.], [1., 2.]])
        np.testing.assert_array_equal(self.m.tile_groups(constant, (constant.min(0), constant.max(0))), [0, 3])

    def test_oriented_edges_never_include_third_identity(self):
        edges = np.array([[0, 1], [0, 2], [1, 2], [1, 3], [2, 3]])
        a, b = self.m.oriented_edges(edges, np.array(['R', 'E', 'C', 'R']), 'R', 'E')
        self.assertEqual(list(zip(a, b)), [(0, 1), (3, 1)])

    def test_full_synthetic_wiring_2000_draws_and_18_opportunities(self):
        rows = []
        mapping = {}
        for case in range(5):
            roi = 'synthetic_roi_'+str(case)
            mapping[roi] = 'synthetic_case_'+str(case)
            for typ in ['R', 'C', 'E']:
                for cell in range(20):
                    value = 19-cell if typ == 'C' else cell
                    rows.append(dict(sample_id=roi, celltype20=typ,
                        **{marker: float(value) for marker in set(self.m.R_MARKERS+self.m.E_MARKERS)}))
        frame = pd.DataFrame(rows)
        coords = np.column_stack([np.tile(np.arange(60), 5), np.zeros(300)])
        def fixed_edges(coordinates, k):
            self.assertEqual(len(coordinates), 60)
            self.assertEqual(k, 10)
            return np.array([(i, 40+i) for i in range(20)]+[(20+i, 40+i) for i in range(20)])
        bundle = dict(frame=frame, native_coordinates=coords, ROI_to_source_patient=mapping,
            native_knn_edges=fixed_edges)
        original_file = self.m.__file__
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root/'scripts').mkdir()
            (root/'private_output').mkdir()
            module_file = root/'scripts'/'synthetic.py'
            module_file.write_text('synthetic')
            bridge = dict(status='COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW',
                independent_CSV_Pandas_node_replay=True, complete_node_member_CRC_checked=True)
            (root/'private_output'/'private_frozen_step_author_node_native_geometry_v1_source_summary.json').write_text(json.dumps(bridge))
            self.m.__file__ = str(module_file)
            try:
                result = self.m.analyze(bundle, dict(regulatory_author_label='R',
                    conventional_CD4_author_label='C', effector_author_label='E'))
            finally:
                self.m.__file__ = original_file
        self.assertEqual(result['family_size'], 18)
        primary = result['family']['ROI_parent_shuffle']
        self.assertEqual(len(primary), 18)
        self.assertTrue(all(row['source_case_units'] == 5 for row in primary))
        self.assertTrue(all(row['effect'] > 0 for row in primary[:9]))
        self.assertTrue(all(row['effect'] < 0 for row in primary[9:]))
        self.assertTrue(all(row['P'] < .01 for row in primary))
        self.assertTrue(all(row['case_bootstrap_95_CI'] is None for row in primary))
        self.assertTrue(all(not row['clinical_twenty_case_coverage_necessary_only'] for row in primary))
        encoded = json.dumps(result, allow_nan=False)
        self.assertNotIn('synthetic_case_', encoded)
        self.assertNotIn('synthetic_roi_', encoded)


if __name__ == '__main__':
    unittest.main()
