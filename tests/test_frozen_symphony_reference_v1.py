"""Real-library tests for frozen-reference and missing-feature invariants."""
from pathlib import Path
import os,sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import numpy as np
    import scipy.sparse as sp
    from annotation_representation_pilot_v1 import normalized_log_counts
    from frozen_symphony_reference_v1 import load_symphonypy,build_reference,map_query,reference_fingerprint
    AVAILABLE=True
except ImportError:
    AVAILABLE=False


@unittest.skipUnless(AVAILABLE and 'ATLAS_SYMPHONY_WHEEL' in os.environ, 'Run in qualified reference runtime')
class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.symphony=load_symphonypy(os.environ['ATLAS_SYMPHONY_WHEEL'],os.environ['ATLAS_SYMPHONY_SHA256'])
        rng=np.random.default_rng(11)
        counts=sp.csr_matrix(rng.poisson(2,(180,60)),dtype=np.int32)
        cls.x,_=normalized_log_counts(counts,np.asarray(counts.sum(axis=1)).ravel())
        cls.genes=['G'+str(i) for i in range(60)]
        cls.params={'components':10,'pca_iterations':2,'theta':1.,'lambda':1.,'max_iter_harmony':3,'threads':4}
        cls.ref,cls.raw,cls.corrected,cls.pca=build_reference(cls.x[:160],np.array(['a']*80+['b']*80),cls.genes,1,cls.params,cls.symphony)

    def test_query_projection_and_reference_immutability(self):
        before=reference_fingerprint(self.ref)
        raw,mapped,valid=map_query(self.ref,self.x[160:],self.genes,np.array(['new']*20),self.symphony)
        self.assertTrue(valid.all());self.assertTrue(np.isfinite(mapped).all())
        np.testing.assert_allclose(raw,self.pca.transform(self.x[160:].toarray()),rtol=2e-5,atol=2e-4)
        self.assertEqual(before,reference_fingerprint(self.ref));self.assertEqual(self.ref.n_obs,0)

    def test_permuted_queries_preserve_correspondence(self):
        order=np.arange(20)[::-1]
        _,a,va=map_query(self.ref,self.x[160:],self.genes,np.array(['new']*20),self.symphony)
        _,b,vb=map_query(self.ref,self.x[160:][order],self.genes,np.array(['new']*20),self.symphony)
        np.testing.assert_array_equal(va,vb[order]);np.testing.assert_allclose(a,b[order],rtol=2e-5,atol=2e-4)

    def test_missing_or_reordered_query_features_fail(self):
        with self.assertRaises(ValueError):map_query(self.ref,self.x[160:,:-1],self.genes[:-1],['new']*20,self.symphony)
        with self.assertRaises(ValueError):map_query(self.ref,self.x[160:],self.genes[::-1],['new']*20,self.symphony)

    def test_single_soft_cluster_small_reference(self):
        ref,_,_,_=build_reference(self.x[:40],np.array(['a']*20+['b']*20),self.genes,1,self.params,self.symphony)
        self.assertEqual(ref.uns['harmony']['K'],1)
        _,mapped,valid=map_query(ref,self.x[160:],self.genes,['new']*20,self.symphony)
        self.assertTrue(valid.any());self.assertTrue(np.isfinite(mapped[valid]).all())


if __name__=='__main__':unittest.main(verbosity=2)
