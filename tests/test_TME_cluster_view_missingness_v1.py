"""Metadata invariants require missing-aware equality, not relaxed equality."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from deliver_TME_hierarchical_cluster_view_v1 import assert_same_metadata
class MissingnessTests(unittest.TestCase):
 def test_equal_object_nan_and_none_and_text_pass(self):assert_same_metadata(np.array([np.nan,None,'unknown',2.],object),np.array([np.nan,None,'unknown',2.],object))
 def test_nan_and_none_remain_distinct(self):
  with self.assertRaises(AssertionError):assert_same_metadata(np.array([np.nan],object),np.array([None],object))
 def test_moved_missing_position_fails(self):
  with self.assertRaises(AssertionError):assert_same_metadata(np.array([np.nan,1.],object),np.array([1.,np.nan],object))
 def test_changed_nonmissing_value_and_dtype_fail(self):
  for a,b in [(np.array([1.,2.]),np.array([1.,3.])),(np.array([1],dtype=np.int32),np.array([1],dtype=np.int64))]:
   with self.assertRaises(AssertionError):assert_same_metadata(a,b)
if __name__=='__main__':unittest.main()
