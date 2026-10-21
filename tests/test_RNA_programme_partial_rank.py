import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

spec=importlib.util.spec_from_file_location('partial_rank',Path(__file__).resolve().parents[1]/'scripts/RNA_programme_partial_rank_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class PartialRank(unittest.TestCase):
    def fixture(self):
        rng=np.random.default_rng(811);n=80;c=rng.normal(size=(n,8));cat=np.array([['D' if i%2 else 'D+T','positive' if i%3 else 'negative'] for i in range(n)])
        x=rng.normal(size=n)+c[:,0];y=.8*x+rng.normal(size=n)+2*c[:,1]
        return x,y,c,cat

    def test_fixed32_order_and_tied_average_ranks(self):
        self.assertEqual(len(m.fixed_slots()),32)
        self.assertEqual(len({tuple(sorted(r.items())) for r in m.fixed_slots()}),32)
        x=np.array([4,2,2,0,0,7.])
        np.testing.assert_array_equal(m.average_rank(x),pd.Series(x).rank(method='average'))

    def test_independent_lstsq_replays_positive_negative_and_permutation(self):
        x,y,c,cat=self.fixture();levels=m.categorical_spec(cat)
        for sign in [-1,1]:
            a=m.fit(x,sign*y,c,cat,levels)
            z=a['z'];rx=pd.Series(x).rank().to_numpy();ry=pd.Series(sign*y).rank().to_numpy()
            ex=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ey=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0]
            self.assertAlmostEqual(a['rho'],np.corrcoef(ex,ey)[0,1],places=12)
            self.assertEqual(int(np.sign(a['rho'])),sign)
            order=np.random.default_rng(33).permutation(len(x));pseudo=z@np.linalg.lstsq(z,ry,rcond=None)[0]+ey[order]
            residual=pseudo-z@np.linalg.lstsq(z,pseudo,rcond=None)[0]
            self.assertAlmostEqual(m.permuted_statistic(a,order),np.corrcoef(ex,residual)[0,1],places=12)

    def test_no_reranking_of_pseudo_outcome(self):
        x,y,c,cat=self.fixture();a=m.fit(x,y,c,cat,m.categorical_spec(cat));order=np.random.default_rng(9).permutation(len(x))
        pseudo=a['yhat']+a['ey'][order];ranked=m.average_rank(pseudo);wrong=ranked-a['q']@(a['q'].T@ranked)
        self.assertGreater(abs(m.permuted_statistic(a,order)-m.correlation(a['ex'],wrong)),1e-5)

    def test_controls_constant_singular_and_residual_df_are_not_dropped(self):
        x,y,c,cat=self.fixture()
        for how in ['constant','singular','df']:
            cc=c.copy();xx=x;yy=y;categories=cat
            if how=='constant':cc[:,2]=1
            if how=='singular':cc[:,2]=cc[:,1]
            if how=='df':xx=x[:20];yy=y[:20];cc=cc[:20];categories=cat[:20]
            result=m.analyze(xx,yy,cc,categories,0,draws=5,min_df=15 if how=='df' else 10,required_bootstrap=4)
            self.assertEqual(result['status'],'NOT_ESTIMABLE');self.assertIsNone(result['P'])
        constantcat=np.repeat([['D','positive']],len(x),axis=0)
        fit=m.fit(x,y,c,constantcat,m.categorical_spec(constantcat));self.assertEqual(fit['constant_factors'],[0,1])

    def test_seed_deterministic_complete_resampling_and_BH_missing(self):
        x,y,c,cat=self.fixture()
        a=m.analyze(x,y,c,cat,4,draws=100,required_bootstrap=90)
        b=m.analyze(x,y,c,cat,4,draws=100,required_bootstrap=90)
        self.assertEqual(a,b);self.assertEqual(a['bootstrap_valid'],100)
        self.assertIsNotNone(a['CI95']);self.assertEqual(a['permutation_draws'],100)
        r=[{'P':.01},{'P':None},{'P':.04},{'P':.5}];m.bh_fixed(r)
        np.testing.assert_allclose([v['BH_fixed_family'] for v in r],[.04,1,.08,2/3])


if __name__=='__main__':unittest.main()
