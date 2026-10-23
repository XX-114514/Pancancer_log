import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import unittest
import numpy as np
import scipy.sparse as sp
from annotation_graph_diagnostics_v1 import (unit_out_graph_state,direct_unit_out_graph_state,
    identity_neighbor_profiles,rare_neighbor_profiles,checked_graph)
from paired_granularity_review_v1 import paired_scores
from compare_lineage_graph_engines_v2 import qualify_expression_units


class GraphDiagnostics(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(31);self.n=160
        self.unit=np.repeat(np.arange(8),20);self.depth=rng.normal(size=self.n)
        self.y=.3*self.depth+rng.normal(size=self.n)
        a=rng.uniform(size=(self.n,self.n));a[a<.9]=0;self.graph=sp.csr_matrix(a)

    def test_independent_weighted_prediction_and_missing(self):
        self.y[1]=np.nan;self.depth[3]=np.nan
        p=unit_out_graph_state(self.graph,self.unit,self.depth,self.y,return_predictions=True)
        direct,base=direct_unit_out_graph_state(self.graph,self.unit,self.depth,self.y)
        np.testing.assert_allclose(p['prediction'],direct,atol=1e-12,equal_nan=True)
        np.testing.assert_allclose(p['baseline'],base,atol=1e-12,equal_nan=True)

    def test_entire_query_unit_excluded_and_self_excluded(self):
        a=unit_out_graph_state(self.graph,self.unit,self.depth,self.y,return_predictions=True)
        altered=self.y.copy();altered[self.unit==0]+=10000
        b=unit_out_graph_state(self.graph,self.unit,self.depth,altered,return_predictions=True)
        np.testing.assert_allclose(a['prediction'][:20],b['prediction'][:20],equal_nan=True)
        c=self.graph.copy();c.setdiag(1e20)
        d=unit_out_graph_state(c,self.unit,self.depth,self.y,return_predictions=True)
        np.testing.assert_allclose(a['prediction'],d['prediction'],equal_nan=True)

    def test_no_other_unit_no_estimate(self):
        x=unit_out_graph_state(self.graph,np.zeros(self.n),self.depth,self.y)
        self.assertEqual(x['n_evaluable'],0)

    def test_pair_uses_common_queries(self):
        a=unit_out_graph_state(self.graph,self.unit,self.depth,self.y,return_predictions=True)
        b={k:v.copy() if isinstance(v,np.ndarray) else v for k,v in a.items()}
        b['prediction'][:12]=np.nan
        p=paired_scores(self.y,self.unit,a,b)
        self.assertEqual(p['second_minus_first'],0.)
        self.assertLess(p['n_common_query_cells'],self.n)

    def test_unknown_neighbors_and_shared_source_gate(self):
        graph=sp.csr_matrix(np.ones((4,4)))
        p=identity_neighbor_profiles(graph,np.array(['a','a','b','b']),np.array(['T','T','T','unresolved']),minimum_source_cells=1)
        self.assertAlmostEqual(p[0]['unknown_neighbor_weight_fraction'],1/3)
        self.assertAlmostEqual(p[0]['same_candidate_fraction_among_known_neighbor_weight'],1.)
        q=identity_neighbor_profiles(graph,np.array(['a','a','b','b']),np.array(['T','T','T','unresolved']),minimum_source_cells=2)
        self.assertTrue(all(r['same_candidate_neighbor_effective_sources'] is None for r in q))

    def test_rare_self_and_missing(self):
        p=rare_neighbor_profiles(sp.csr_matrix(np.ones((3,3))),np.array(['a']*3),np.array([1,0,np.nan]))[0]
        self.assertEqual(p['n_measured'],2);self.assertEqual(p['candidate_neighbor_weight_fraction'],0.)

    def test_qualified_units_preserve_unicode_and_leading_zeros(self):
        values=np.asarray(['001','A-2','study::sample'],dtype='U30')
        a=qualify_expression_units('source_a',values)
        self.assertEqual(a.tolist(),['source_a::001','source_a::A-2','source_a::study::sample'])
        self.assertFalse(np.isin(a,qualify_expression_units('source_b',values)).any())

    def test_invalid_graph_and_axes(self):
        with self.assertRaises(ValueError):checked_graph(sp.csr_matrix([[-1,0],[0,0]]))
        with self.assertRaises(ValueError):unit_out_graph_state(self.graph,[],[],[])


if __name__=='__main__':unittest.main()
