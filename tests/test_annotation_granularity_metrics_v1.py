from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import numpy as np
    from annotation_granularity_metrics_v1 import granularity,unit_out_state,gamma_group_retention
    AVAILABLE=True
except ImportError:AVAILABLE=False


@unittest.skipUnless(AVAILABLE,'Actual scientific environment required')
class GranularityTests(unittest.TestCase):
    def test_group_numbers_are_namespace_invariant_and_small_groups_explicit(self):
        a=np.array([0]*80+[1]*20+[-1]*10)
        b=np.array([400]*80+[100]*20+[-1]*10)
        self.assertEqual(granularity(a),granularity(b))
        self.assertAlmostEqual(granularity(a)['effective_group_count'],np.exp(-.8*np.log(.8)-.2*np.log(.2)))
        self.assertEqual(granularity([-1,-1])['status'],'unpartitioned')

    def test_shared_groups_predict_other_unit_state(self):
        labels=np.tile(np.repeat([0,1],80),2);units=np.repeat([0,1],160)
        y=labels.astype(float)*2
        out=unit_out_state(labels,units,np.ones(len(y)),y)
        self.assertEqual(out['n_evaluable'],320)
        self.assertAlmostEqual(out['RNA_prediction_improvement_over_depth'],1.)

    def test_unit_exclusive_groups_cannot_leak_their_own_RNA(self):
        labels=np.repeat([0,1],160);units=labels.copy();y=labels.astype(float)
        out=unit_out_state(labels,units,np.ones(320),y)
        self.assertEqual(out['n_evaluable'],0)

    def test_heldout_unit_outcome_does_not_change_its_predictions(self):
        rng=np.random.default_rng(42)
        labels=np.tile(np.repeat([0,1],80),3);units=np.repeat([0,1,2],160)
        depth=rng.normal(8,1,len(labels));y=.4*depth+labels+rng.normal(0,.1,len(labels))
        a=unit_out_state(labels,units,depth,y,return_predictions=True)
        y[units==1]+=rng.normal(0,4,160)
        b=unit_out_state(labels,units,depth,y,return_predictions=True)
        np.testing.assert_allclose(a['prediction'][units==1],b['prediction'][units==1],atol=1e-12)
        np.testing.assert_allclose(a['baseline'][units==1],b['baseline'][units==1],atol=1e-12)

    def test_gamma_singletons_do_not_get_perfect_retention(self):
        out=gamma_group_retention(np.array([0,1,1]),np.array([1,1,0]))
        self.assertEqual(out['n_singleton_candidates'],1)
        self.assertEqual(out['candidate_same_group_fraction_excluding_self'],0.)

    def test_graph_reuse_at_resolution_one_matches_original_API(self):
        from compare_annotation_granularity_v1 import graph_grid
        from annotation_representation_pilot_v1 import cluster
        from sklearn.metrics import adjusted_rand_score
        rng=np.random.default_rng(91);x=rng.normal(size=(160,6)).astype(np.float32)
        x[:80,0]+=5
        cfg={'threads':4,'neighbors':20,'resolution':1.}
        original,_=cluster(x,91,cfg)
        grid=graph_grid(x,91,cfg,[.5,1.,2.])
        self.assertEqual(adjusted_rand_score(original,grid['1.0']),1.)


if __name__=='__main__':unittest.main()
