"""Full author adjacency routing and complete edge identity checks."""
import csv
import hashlib
import importlib.util
import io
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fixtures = load(ROOT/'tests/test_melanoma_IMC_edge_centered_state_neighbors.py','existing_centered_fixtures')


class AuthorAdjacencyMath(fixtures.EdgeCenteredOrganization):
    @classmethod
    def setUpClass(cls):
        cls.m = load(ROOT/'scripts/analyze_melanoma_IMC_author_graph_states_v1.py','actual_author_core')
        ordinary = cls.m.analyze
        def author_fixture(bundle, parameters):
            # Existing independent full-wiring fixture supplies the intended
            # author edges. Native kNN must never be called by this core.
            bundle['author_neighbor_edges'] = {
                roi:bundle['native_knn_edges'](bundle['native_coordinates'][indices],10)
                for roi,indices in bundle['frame'].groupby('sample_id').indices.items()}
            def forbidden(*args):raise AssertionError('Native nearest neighbors must not replace author adjacency')
            bundle['native_knn_edges'] = forbidden
            return ordinary(bundle,parameters)
        cls.m.analyze = author_fixture


class CompleteAuthorEdges(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = load(ROOT/'scripts/measure_melanoma_IMC_complete_author_graph_states_v1.py','actual_author_stream')

    def fixture(self):
        frame = pd.DataFrame(dict(sample_id=['a']*4+['b']*2,celltype20=['R','C','E','E','R','E']))
        target = pd.MultiIndex.from_arrays([frame['sample_id'],[str(i) for i in range(6)]])
        locations = np.arange(6)
        lookup = {key:index for index,key in enumerate(target)}
        return frame,target,locations,lookup

    def test_all_endpoints_independent_key_resolution_and_alias(self):
        frame,target,locations,lookup = self.fixture()
        images,a,b=self.m.resolve_rows([['','raw','0','2']],target,locations,lookup,{'raw':'a'},{'a','b'})
        self.assertEqual((images.tolist(),a.tolist(),b.tolist()),(['a'],[0],[2]))

    def test_unknown_ROI_and_missing_endpoint_fail(self):
        frame,target,locations,lookup=self.fixture()
        for rows in [[['','missing','0','2']],[['','a','0','999']]]:
            with self.assertRaises(AssertionError):self.m.resolve_rows(rows,target,locations,lookup,{}, {'a','b'})

    def test_duplicate_directions_self_edges_and_unmapped_margin_not_inflated(self):
        frame,target,locations,lookup=self.fixture()
        rows=[['','a','0','2'],['','a','2','0'],['','a','1','3'],['','a','0','0'],['','b','4','5']]
        edges,q=self.m.selected_edges(iter(rows),frame,target,locations,{}, {'a':'case'},['R','C','E'])
        self.assertEqual(edges['a'].tolist(),[[0,2],[1,3]])
        self.assertEqual(q['all_input_edges_resolved_twice'],5)
        self.assertEqual(q['selected_raw_cross_role_edges'],3)
        self.assertEqual(q['selected_unique_undirected_cross_role_edges'],2)

    def test_last_batch_is_resolved_and_missing_edges_remain_empty(self):
        frame,target,locations,lookup=self.fixture()
        edges,q=self.m.selected_edges(iter([['','b','4','5']]),frame,target,locations,{}, {'a':'case'},['R','C','E'])
        self.assertEqual(edges['a'].shape,(0,2));self.assertEqual(q['all_input_edges_resolved_twice'],1)

    def test_full_byte_digest_and_CSV_replay(self):
        blob=b',image,first,second\n,a,0,2\n,a,1,3\n'
        stream=self.m.DigestReader(io.BytesIO(blob))
        rows=list(csv.reader(io.TextIOWrapper(io.BufferedReader(stream),encoding='utf-8-sig',newline='')))
        self.assertEqual(len(rows),3)
        self.assertEqual(stream.count,len(blob))
        self.assertEqual(stream.digest.hexdigest(),hashlib.sha256(blob).hexdigest())


if __name__ == '__main__':unittest.main()
