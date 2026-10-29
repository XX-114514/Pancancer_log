"""Regression cases for legacy object payloads in a primitive NPZ audit."""
import pathlib,sys,tempfile,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from audit_TME_UMAP_marker_review_v1 import load_required_display_axes

class ReadbackKeys(unittest.TestCase):
 def payload(self):
  return dict(UMAP_Harmony=np.ones((2,2)),major=np.array(['T','B']),fine=np.array(['Tcell','B']),source=np.array(['a','b']),raw=np.ones((2,3)),corrected=np.ones((2,3)),genes=np.ones((2,4)))
 def test_unused_legacy_object_is_not_deserialized(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d)/'axes.npz';x=self.payload();x['baseline_major']=np.array(['T',None],object);np.savez(p,**x)
   actual=load_required_display_axes(p);self.assertNotIn('baseline_major',actual);np.testing.assert_array_equal(actual['major'],['T','B'])
 def test_required_object_array_still_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=pathlib.Path(d)/'axes.npz';x=self.payload();x['major']=np.array(['T','B'],object);np.savez(p,**x)
   with self.assertRaises(ValueError):load_required_display_axes(p)

if __name__=='__main__':unittest.main()
