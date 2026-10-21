"""Independent synthetic positive, negative, confounded and missing cases."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
from CM02_marked_spatial_core_v1 import (rankz,residualize,graph_edges,weighted_corr,
    permutation_groups,permutations,permuted_corr,cohort_summary,bh_fixed_family,paired_roi_contrast)


class Tests(unittest.TestCase):
    def test_average_rank_ties_and_constant(self):
        x=rankz([1,1,4,7]);expected=np.array([1.5,1.5,3,4]);expected=(expected-expected.mean())/expected.std()
        np.testing.assert_allclose(x,expected);self.assertTrue(np.isnan(rankz([2,2])).all())
    def test_fixed_control_confounded_signal_removed(self):
        rng=np.random.default_rng(8);x=rng.normal(size=(40,4));noise=rng.normal(size=40)
        y=x@np.array([20,5,3,-8])+noise;r,reason=residualize(y,x)
        self.assertEqual(reason,'ESTIMABLE');self.assertLess(abs(np.dot(r,x[:,0])),1e-10)
        self.assertGreater(np.corrcoef(r,noise)[0,1],.9)
    def test_collinearity_and_df_never_drop_controls(self):
        x=np.column_stack([np.arange(12),np.arange(12)])
        self.assertEqual(residualize(np.arange(12),x)[1],'FIXED_DESIGN_NOT_FULL_RANK')
        rng=np.random.default_rng(5);x=rng.normal(size=(10,4))
        self.assertEqual(residualize(rng.normal(size=10),x)[1],'RESIDUAL_DF_BELOW_FIXED_MINIMUM')
    def test_positive_negative_edge_mark_and_brute_force(self):
        rng=np.random.default_rng(12);s=rng.normal(size=(40,3));xy=np.column_stack([np.arange(40)*30,np.arange(40)%2*10])
        edges,reason=graph_edges(xy,xy+1,5);self.assertEqual(reason,'ESTIMABLE')
        t=np.column_stack([s[:,0],-s[:,0]]);r=weighted_corr(s,t,edges)
        self.assertAlmostEqual(r[0,0],1);self.assertAlmostEqual(r[0,1],-1)
        np.testing.assert_allclose(r,np.corrcoef(np.column_stack([s,t]),rowvar=False)[:3,3:],atol=1e-12)
        self.assertIsNone(graph_edges(xy,xy+100000,5)[0])
    def test_joint_stratified_permutations_and_vectorized_replay(self):
        rng=np.random.default_rng(2);a=rng.normal(size=(40,3));b=rng.normal(size=(40,2))
        xy=np.column_stack([np.arange(40)*30,np.arange(40)%2*10]);g,_=graph_edges(xy,xy+1,5)
        groups=[np.arange(20),np.arange(20,40)]
        p=permutations(groups,40,20,np.random.default_rng(3));q=permutations(groups,40,20,np.random.default_rng(4))
        for row in p:self.assertEqual(set(row[:20]),set(range(20)))
        np.testing.assert_allclose(np.cov(a[p[0]],rowvar=False),np.cov(a,rowvar=False))
        got=permuted_corr(a,b,g,p,q)
        expected=np.array([weighted_corr(a[pi],b[qi],g) for pi,qi in zip(p,q)])
        np.testing.assert_allclose(got,expected,atol=1e-12)
    def test_grid_and_unknown_interfaces_remain_distinct(self):
        frame=pd.DataFrame(dict(Location_Center_X=[0,1,2,3],Location_Center_Y=[0,1,2,3],
            is_interface=[0,0,-1,-1],is_perivascular=[0,0,0,0]))
        groups,n=permutation_groups(frame,(np.array([0,0]),np.array([10,10])),2,True)
        self.assertEqual(n,4);self.assertEqual(sorted(map(len,groups)),[2,2])
    def test_source_case_inference_and_insufficient_case_slot(self):
        rng=np.random.default_rng(5);null=rng.normal(0,.1,size=(25,2000));obs=np.full(25,.7)
        got=cohort_summary(obs,null,list(range(25)))
        self.assertLess(got['conditional_randomization_P'],.01);self.assertGreater(got['CI95'][0],0)
        small=cohort_summary(obs[:7],null[:7],list(range(7)))
        self.assertNotIn('conditional_randomization_P',small);self.assertIsNone(small['CI95'])
    def test_fixed_BH_denominator_includes_missing(self):
        results={'a':{'conditional_randomization_P':.01},'b':{'conditional_randomization_P':.04},'missing':{'n_cases':0}}
        bh_fixed_family(results,12);self.assertAlmostEqual(results['a']['BH_fixed_family'],.12)
        self.assertAlmostEqual(results['b']['BH_fixed_family'],.24);self.assertNotIn('BH_fixed_family',results['missing'])
    def test_paired_contrast_requires_identical_ROI(self):
        a=[('shared',.5,np.array([0.,.1])),('only_a',.99,np.array([0.,0.]))]
        b=[('shared',.2,np.array([0.,.2])),('only_b',-.99,np.array([0.,0.]))]
        delta,null=paired_roi_contrast(a,b);self.assertAlmostEqual(delta,.3)
        np.testing.assert_allclose(null,[0,-.1])
        self.assertIsNone(paired_roi_contrast(a[1:],b[1:]))


if __name__=='__main__':unittest.main()
