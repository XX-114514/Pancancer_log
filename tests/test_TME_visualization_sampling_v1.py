import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_visualization_sampling_v1 import balanced_indices
class DisplaySamplingTests(unittest.TestCase):
 def test_rare_cluster_and_projection_mask_preserved(self):
  codes=np.array([0]*100+[1]*3+[2]);eligible=np.ones(len(codes),bool);eligible[-1]=False;p=balanced_indices(codes,eligible,20,3,29);self.assertEqual(len(p),20);self.assertTrue(eligible[p].all());self.assertEqual(sum(codes[p]==1),3);self.assertEqual(len(set(p)),len(p))
 def test_reproducible_and_full_scope_when_small(self):
  codes=np.repeat(np.arange(3),4);a=balanced_indices(codes,np.ones(12,bool),30,2,9);np.testing.assert_array_equal(a,np.arange(12));np.testing.assert_array_equal(a,balanced_indices(codes,np.ones(12,bool),30,2,9))
 def test_rare_coverage_may_explicitly_exceed_nominal_cap(self):
  codes=np.arange(10);p=balanced_indices(codes,np.ones(10,bool),3,3,7);np.testing.assert_array_equal(p,np.arange(10))
if __name__=='__main__':unittest.main()
