from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from project_frozen_query_domain_v1 import project_domain,project_raw_context


class DomainTests(unittest.TestCase):
    def setUp(self):
        self.ref=ad.AnnData(sp.csr_matrix((0,3)),var=pd.DataFrame({'mean':np.zeros(3),'std':np.ones(3)},index=['A','B','C']));self.ref.varm['PCs']=np.eye(3,dtype=np.float32)
        self.ref.uns['harmony']={'C':np.eye(3),'Nr':np.ones(3),'sigma':np.ones(3)};self.calls=[]
    def mapper(self,ref,x,genes,source,symphony):
        self.calls.append(x.toarray());raw=x.toarray();valid=(raw.max(axis=1)>0)&(np.linalg.norm(raw,axis=1)>1e-10);mapped=np.full(raw.shape,np.nan);mapped[valid]=raw[valid]-raw[valid].mean(axis=0);return raw,mapped,valid
    def test_excluded_and_no_signal_never_adapt_and_domain_not_routed(self):
        x=sp.csr_matrix((np.array([1,2,3,0,1,1,2,1],np.float32),np.array([0,1,0,0,0,2,1,2]),np.array([0,2,3,4,6,8])),shape=(5,3))
        r=project_domain(self.ref,x,['A','B','C'],'s',[1,1,1,0,1],[1,0,1,1,1],[1,1,1,1,1],mapper=self.mapper)
        self.assertEqual(len(self.calls),1);np.testing.assert_array_equal(self.calls[0],x[[0,4]].toarray());self.assertTrue(np.isfinite(r['raw_PCA'][1]).all());self.assertTrue(np.isnan(r['raw_PCA'][3]).all());self.assertTrue(np.isnan(r['candidate_Harmony'][[1,2,3]]).all());np.testing.assert_array_equal(r['adaptation_included'],[1,0,0,0,1])
    def test_raw_context_never_adapts_and_retains_positive_RNA_projection(self):
        x=sp.csr_matrix([[1,2,0],[0,0,0],[0,0,0]],dtype=np.float32);r=project_raw_context(x,np.array([.5,.5,0]),np.eye(3),[1,1,0])
        np.testing.assert_allclose(r['raw_PCA'][:2],[[.5,1.5,0],[-.5,-.5,0]]);self.assertTrue(np.isnan(r['raw_PCA'][2]).all());self.assertTrue(np.isnan(r['candidate_Harmony']).all());self.assertFalse(r['adaptation_included'].any());np.testing.assert_array_equal(r['model_gene_signal'],[1,0,0])
    def test_complete_batch_permutation_consistency(self):
        x=sp.csr_matrix([[1,2,0],[3,0,1],[0,1,2],[1,0,1]],dtype=np.float32);mask=np.ones(4,bool);order=np.array([2,0,3,1]);inverse=np.argsort(order)
        a=project_domain(self.ref,x,['A','B','C'],'s',mask,mask,mask,mapper=self.mapper);b=project_domain(self.ref,x[order],['A','B','C'],'s',mask,mask,mask,mapper=self.mapper)
        np.testing.assert_allclose(a['candidate_Harmony'],b['candidate_Harmony'][inverse]);self.assertEqual(len(self.calls),2)


if __name__=='__main__':unittest.main()
