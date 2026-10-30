"""Actual standard-PCA frozen-mapping checks using the existing pinned wheel."""
import sys,json,unittest
from pathlib import Path
import numpy as np
import scipy.sparse as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from standard_annotation_order_primitives_v1 import select_features
from frozen_symphony_reference_standard_v2 import load_symphonypy,build_reference,map_query,reference_fingerprint
class StandardReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cfg=json.loads((Path(__file__).resolve().parents[1]/'config/source_holdout_spec_v1.json').read_text())
        cls.symphony=load_symphonypy(cfg['wheel_path'],cfg['wheel_SHA256'])
        rng=np.random.default_rng(2026103005);counts=sp.csr_matrix(rng.poisson(np.geomspace(.3,6,300),size=(140,300)).astype(np.float32));genes=['G'+str(i) for i in range(300)]
        source=np.repeat(['a','b'],60);selected,_=select_features(counts[:120],source,50,genes,120)
        totals=np.asarray(counts.sum(axis=1)).ravel();x=counts[:,selected].multiply((10000/totals)[:,None]).tocsr();x.data=np.log1p(x.data)
        cls.genes=[genes[i] for i in selected];cls.query=x[120:];cls.qsource=np.repeat('query',20)
        cls.ref,cls.raw,cls.corrected,cls.pca=build_reference(x[:120],source,cls.genes,2026102301,dict(components=8,pca_iterations=3,theta=2.,**{'lambda':1.},max_iter_harmony=20,threads=2),cls.symphony)
    def test_actual_standard_centering_and_frozen_query_mapping(self):
        before=reference_fingerprint(self.ref)
        raw,mapped,ok=map_query(self.ref,self.query,self.genes,self.qsource,self.symphony)
        expected=(self.query.toarray()-self.ref.var['mean'].to_numpy())@self.ref.varm['PCs']
        np.testing.assert_allclose(raw,expected,rtol=2e-5,atol=2e-4)
        self.assertTrue(ok.any());self.assertTrue(np.isfinite(mapped[ok]).all())
        self.assertEqual(before,reference_fingerprint(self.ref))
    def test_changed_or_unmeasured_query_feature_axis_is_rejected(self):
        before=reference_fingerprint(self.ref)
        with self.assertRaises(ValueError):map_query(self.ref,self.query,self.genes[::-1],self.qsource,self.symphony)
        with self.assertRaises(ValueError):map_query(self.ref,self.query[:,:-1],self.genes[:-1],self.qsource,self.symphony)
        self.assertEqual(before,reference_fingerprint(self.ref))
if __name__=='__main__':unittest.main()
