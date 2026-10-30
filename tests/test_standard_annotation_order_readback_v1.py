"""Independent readback examples, including float deltas and missing evidence."""
import sys,unittest
from pathlib import Path
import numpy as np
import scipy.sparse as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from standard_annotation_order_readback_v1 import gather_windows,independent_conditional_metrics,compare_record
from standard_annotation_order_primitives_v1 import conditional_readout_metrics
class Group(dict):
    pass
class ReadbackTests(unittest.TestCase):
    def test_sparse_windows_preserve_float_zero_order_and_duplicate_rows(self):
        x=sp.csr_matrix(np.array([[0,.25,1],[0,0,0],[.5,0,2],[0,3,0]],float))
        g=Group(data=x.data,indices=x.indices,indptr=x.indptr);g.attrs={'shape':x.shape}
        rows=np.array([3,0,3,1,2]);np.testing.assert_array_equal(gather_windows(g,rows,width=2).toarray(),x[rows].toarray())
        with self.assertRaises(ValueError):gather_windows(g,np.array([4]))
    def test_independent_moran_correlations_constant_and_missing_cases(self):
        n=100;rng=np.random.default_rng(2026103004);z=rng.normal(size=(n,5));depth=rng.integers(500,5000,n)
        states={'positive':z[:,0]+rng.normal(0,.1,n),'opposite':-z[:,1],
                'constant':np.ones(n),'missing':np.full(n,np.nan)};rare=(z[:,0]>1).astype(float)
        args=(z,states,depth,rare,np.arange(n),7)
        compare_record(independent_conditional_metrics(*args),conditional_readout_metrics(*args),atol=1e-10)
        rare[0]=np.nan;self.assertEqual(independent_conditional_metrics(*args)['gamma_delta_RNA_candidate']['status'],'unmeasured_not_negative')
    def test_tampered_numeric_and_coverage_record_rejected(self):
        with self.assertRaises(ValueError):compare_record({'effect':.2},{'effect':.3})
        with self.assertRaises(ValueError):compare_record({'n':20},{'n':21})
        with self.assertRaises(ValueError):compare_record({'status':'estimated'},{'status':'unmeasured'})
if __name__=='__main__':unittest.main()
