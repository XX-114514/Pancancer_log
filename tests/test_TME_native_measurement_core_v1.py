from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
 import numpy as np
 import scipy.sparse as sp
 from TME_native_measurement_core_v1 import canonical_native,context_values,update_group_stats,expression_candidates
 AVAILABLE=True
except ImportError:
 AVAILABLE=False

@unittest.skipUnless(AVAILABLE,'Use the existing scientific environment; production rejects skips')
class NativeCoreTests(unittest.TestCase):
 def test_full_RNA_denominator_and_unfitted_zero_row(self):
  x=canonical_native(sp.csr_matrix([[2,0],[0,4],[0,0],[1,2]]));s={'n':np.zeros(2,np.int64),**{k:np.zeros((2,2)) for k in ['cp_sum','log_sum','detected']}}
  update_group_stats(s,x,np.array([10.,20.,0.,10.]),np.array([0,0,-1,1]));np.testing.assert_array_equal(s['n'],[2,1]);np.testing.assert_array_equal(s['cp_sum'],[[2000,2000],[1000,2000]])
  with self.assertRaises(ValueError):update_group_stats(s,x,np.array([10.,20.,0.,10.]),np.array([0,0,1,1]))
 def test_duplicate_slots_coalesced_before_logarithm(self):
  x=sp.csr_matrix((np.array([2,3,0]),np.array([0,0,1]),np.array([0,3])),shape=(1,2));x=canonical_native(x);self.assertEqual(x.nnz,1)
  s={'n':np.zeros(1,np.int64),**{k:np.zeros((1,2)) for k in ['cp_sum','log_sum','detected']}}
  update_group_stats(s,x,np.array([10.]),np.array([0]));self.assertAlmostEqual(s['log_sum'][0,0],np.log1p(5000));self.assertEqual(s['detected'][0,0],1)
 def test_missing_feature_is_not_zero_and_alias_collision_rejected(self):
  x=canonical_native(sp.csr_matrix([[0,2],[1,0]]));v,m=context_values(x,[0,None,1]);self.assertTrue(np.isnan(v[:,1]).all());self.assertEqual(v[0,0],0);np.testing.assert_array_equal(m,[True,False,True])
  with self.assertRaises(ValueError):context_values(x,[1,1])
 def test_descriptive_fallback_and_same_marker_filter(self):
  r=expression_candidates(['RPL3','A','B'],[20,5,1],[1,.7,.4],[1,.5,2],[.5,.2,.8]);self.assertEqual([v['gene'] for v in r],['A']);self.assertIn('not_inferential_DEG',r[0]['marker_evidence_type'])
  r=expression_candidates(['RPL3','A'],[100,5],[1,.7],None,None);self.assertEqual(r[0]['gene'],'A');self.assertIn('without_comparator',r[0]['marker_evidence_type'])

if __name__=='__main__':unittest.main()
