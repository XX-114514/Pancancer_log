"""Synthetic complete-key, label, ambiguity and geometry checks."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd


class AuthorGraphIdentity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1]/'scripts'/'qualify_melanoma_IMC_native_graph_key_bridge_v1.py'
        spec = importlib.util.spec_from_file_location('bridge_under_test', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def setUp(self):
        self.frame = pd.DataFrame(dict(sample_id=['roiA', 'roiA', 'roiB'],
            X=['2', '1', '1'], celltype20=['Treg', 'CD8', 'Treg'],
            **{'Unnamed: 0': ['global0', 'global1', 'global2']}))
        self.nodes = pd.DataFrame({'': ['0', '1', '2'],
            'Imageid': ['rawA', 'rawA', 'roiB'],
            'cell_id': ['1', '2', '1'], 'type': ['CD8', 'Treg', 'Treg']})
        self.aliases = {'rawA': 'roiA'}

    def test_complete_bijection_and_labels(self):
        _, locations, fields, checks = self.m.bind_nodes(self.frame, self.nodes, self.aliases)
        self.assertEqual(locations.tolist(), [1, 0, 2])
        self.assertEqual(fields, ['X'])
        self.assertFalse(checks[1]['complete_literal_key_bijection'])

    def test_equivalent_duplicate_identity_fields_allowed(self):
        self.frame['Unnamed: 0'] = self.frame['X']
        self.assertEqual(self.m.bind_nodes(self.frame, self.nodes, self.aliases)[2], ['X', 'Unnamed: 0'])

    def test_different_complete_maps_are_ambiguous(self):
        self.frame['Unnamed: 0'] = ['1', '2', '1']
        self.frame['celltype20'] = 'Treg'
        self.nodes['type'] = 'Treg'
        with self.assertRaisesRegex(AssertionError, 'ambiguous'):
            self.m.bind_nodes(self.frame, self.nodes, self.aliases)

    def test_label_conflict_rejects_identity(self):
        self.nodes.loc[0, 'type'] = 'wrong'
        with self.assertRaisesRegex(AssertionError, 'no_complete'):
            self.m.bind_nodes(self.frame, self.nodes, self.aliases)

    def test_duplicate_author_keys_rejected(self):
        self.nodes.loc[1, 'cell_id'] = '1'
        with self.assertRaisesRegex(AssertionError, 'duplicate_author'):
            self.m.bind_nodes(self.frame, self.nodes, self.aliases)

    def test_unknown_ROI_not_inferred(self):
        self.nodes.loc[0, 'Imageid'] = 'unknown'
        with self.assertRaisesRegex(AssertionError, 'unmapped'):
            self.m.bind_nodes(self.frame, self.nodes, self.aliases)

    def test_different_node_universe_rejected(self):
        with self.assertRaisesRegex(AssertionError, 'different_node'):
            self.m.bind_nodes(self.frame, self.nodes.iloc[:2], self.aliases)

    def test_edge_geometry_uses_mapped_main_rows_and_ROI(self):
        target, locations, _, _ = self.m.bind_nodes(self.frame, self.nodes, self.aliases)
        result = self.m.edge_prefix_summary([['0', 'rawA', '1', '2']], target, locations,
            np.array([[0, 0], [3, 4], [100, 100]]), self.aliases, {'roiA', 'roiB'})
        self.assertEqual(result['native_distance_quantiles'], [5.0]*5)
        self.assertEqual(result['self_edges'], 0)
        self.assertTrue(result['physical_micron_scale_not_inferred'])

    def test_missing_edge_endpoint_rejected(self):
        target, locations, _, _ = self.m.bind_nodes(self.frame, self.nodes, self.aliases)
        with self.assertRaisesRegex(AssertionError, 'edge_endpoint'):
            self.m.edge_prefix_summary([['0', 'rawA', '1', '99']], target, locations,
                np.array([[0, 0], [3, 4], [100, 100]]), self.aliases, {'roiA', 'roiB'})


if __name__ == '__main__':
    unittest.main()
