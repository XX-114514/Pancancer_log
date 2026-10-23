import sys
from pathlib import Path
import unittest
import numpy as np
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from source_lineage_discovery_v1 import fit_branches


class BranchTests(unittest.TestCase):
    def setUp(self):
        self.source=np.repeat(['s1','s2'],12)
        self.route=np.tile(np.repeat(['T_NK','Myeloid','unresolved_context'],4),2)
        self.good=np.ones(24,bool);self.good[-1]=False
        self.x=sp.csr_matrix(np.eye(24,dtype=np.float32));self.calls=[]
    def fit(self,rows,correct):
        self.calls.append((rows.copy(),correct))
        return {'1.0':np.arange(len(rows))%2}
    def run_branches(self,minimum=4):
        return fit_branches(self.x,self.source,self.route,self.good,41,{'minimum_fit_group':minimum},[1.],fit=self.fit)
    def test_source_and_lineage_boundaries_and_zero_RNA(self):
        result,_,_=self.run_branches()
        for m in ['A','C']:
            labels=result[m]['1.0'];self.assertTrue(np.array_equal(labels>=0,self.good))
            for code in np.unique(labels[self.good]):
                keep=labels==code;self.assertEqual(len(np.unique(self.route[keep])),1)
                if m=='A':self.assertEqual(len(np.unique(self.source[keep])),1)
    def test_unknown_pool_keeps_local_partition(self):
        result,_,_=self.run_branches()
        for source in ['s1','s2']:
            m=(self.source==source)&(self.route=='unresolved_context')&self.good
            self.assertEqual(adjusted_rand_score(result['A']['1.0'][m],result['C']['1.0'][m]),1.)
        self.assertFalse(any(correct and 'unresolved_context' in self.route[rows] for rows,correct in self.calls))
    def test_small_branches_preserved_as_coarse(self):
        result,fine,records=self.run_branches(minimum=5)
        self.assertFalse(fine['A'].any());self.assertTrue((result['A']['1.0'][self.good]>=0).all())
        self.assertTrue(any(r['status']=='small_branch_coarse_only' for r in records))
    def test_axis_mismatch_rejected_before_fitting(self):
        with self.assertRaises(ValueError):
            fit_branches(self.x,self.source,self.route[:-1],self.good,41,{'minimum_fit_group':4},[1.],fit=self.fit)
        self.assertFalse(self.calls)


if __name__=='__main__':unittest.main()
