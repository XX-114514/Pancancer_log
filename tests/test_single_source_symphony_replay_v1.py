from pathlib import Path
import os,sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from frozen_symphony_reference_v1 import load_symphonypy
from single_source_symphony_replay_v1 import fit_single_source_offsets,apply_single_source_offsets


class QueryReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        load_symphonypy(os.environ['ATLAS_SYMPHONY_WHEEL'],os.environ['ATLAS_SYMPHONY_SHA256'])
        from symphonypy._utils import _assign_clusters,_correct_query
        cls.assign=staticmethod(_assign_clusters);cls.correct=staticmethod(_correct_query)
    def check_case(self,k,dtype):
        rng=np.random.default_rng(51);x=rng.normal(size=(211,7)).astype(np.float32);x[:,-1]=np.abs(x[:,-1])+1;centers=rng.normal(size=(k,7)).astype(dtype);nr=rng.uniform(30,70,k).astype(dtype);sigma=np.full(k,.1,dtype=dtype)
        r=self.assign(x,sigma,centers/np.linalg.norm(centers,axis=1,keepdims=True),k);expected=self.correct(x,np.ones((2,len(x))),r,nr,centers,np.diag([0.,1.]))
        a=fit_single_source_offsets(x,centers,nr,sigma,block_rows=17);b=fit_single_source_offsets(x,centers,nr,sigma,block_rows=len(x));np.testing.assert_allclose(a,b,rtol=1e-9,atol=1e-9)
        replay=apply_single_source_offsets(x,centers,sigma,a,block_rows=23);np.testing.assert_allclose(replay,expected,rtol=2e-5,atol=2e-5)
    def test_multiple_clusters_both_center_dtypes(self):
        for dtype in [np.float32,np.float64]:self.check_case(6,dtype)
    def test_one_cluster_and_chunk_invariant_whole_domain(self):self.check_case(1,np.float32)


if __name__=='__main__':unittest.main()
