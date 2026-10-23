import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from annotation_granularity_metrics_v1 import unit_out_state
from paired_granularity_review_v1 import direct_unit_out,paired_scores


class PairedReviewTests(unittest.TestCase):
    def test_independent_predictor_and_constant_depth(self):
        rng=np.random.default_rng(41);lab=np.tile([0,1,2],240);units=np.repeat(np.arange(6),120)
        for x in [rng.normal(size=len(lab)),np.ones(len(lab))*4]:
            y=lab+2*x+rng.normal(size=len(lab))*.2
            a=unit_out_state(lab,units,x,y,return_predictions=True);p,b=direct_unit_out(lab,units,x,y)
            np.testing.assert_allclose(a['prediction'],p,atol=1e-10,equal_nan=True)
            np.testing.assert_allclose(a['baseline'],b,atol=1e-10,equal_nan=True)

    def test_common_mask_does_not_reward_low_coverage(self):
        y=np.arange(200,dtype=float);u=np.repeat([0,1],100)
        a={'prediction':y+1,'baseline':np.zeros(200)}
        b={'prediction':y+1,'baseline':np.zeros(200)}
        b['prediction'][100:]=np.nan;a['prediction'][100:]=999
        r=paired_scores(y,u,a,b)
        self.assertEqual(r['second_minus_first'],0);self.assertEqual(r['n_common_query_cells'],100)
        self.assertFalse(r['high_common_coverage_at_least90pct'])

    def test_incompatible_baselines_rejected(self):
        y=np.arange(200);u=np.repeat([0,1],100)
        with self.assertRaises(AssertionError):
            paired_scores(y,u,{'prediction':y,'baseline':np.zeros(200)},
                          {'prediction':y,'baseline':np.ones(200)})

    def test_unit_exclusive_groups_remain_missing(self):
        u=np.repeat([0,1],100);y=np.arange(200,dtype=float)
        a=unit_out_state(u,u,np.ones(200),y,return_predictions=True)
        p,b=direct_unit_out(u,u,np.ones(200),y)
        self.assertTrue(np.isnan(p).all());self.assertTrue(np.isnan(b).all())
        self.assertEqual(paired_scores(y,u,a,a)['status'],'insufficient_common_query_cells')


if __name__=='__main__':unittest.main()
