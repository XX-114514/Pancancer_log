import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_assay_aware_T_anchor_policy_v1 import assay_aware_flags
from TME_cluster_annotation_policy_v1 import GATE_ORDER,coherent_flags
class AssayAnchorTests(unittest.TestCase):
 def test_chain_missing_positive_complete_alternative(self):
  genes=['CD3D','CD3E','CD2','LCK','LAT','CD8A','CD8B'];f,a,m=assay_aware_flags(np.ones((8,len(genes))),genes);self.assertTrue(a['T']);self.assertTrue(f[:,0].all());self.assertTrue(f[:,GATE_ORDER.index('CD8_T')].all());self.assertIn('alternate',m)
 def test_measured_negative_chain_never_rescued(self):
  genes=['CD3D','CD3E','CD2','LCK','LAT','TRAC'];x=np.ones((5,len(genes)));x[:,-1]=0;f,a,m=assay_aware_flags(x,genes);self.assertTrue(a['T']);self.assertFalse(f[:,0].any());self.assertNotIn('alternate',m)
 def test_one_signal_not_enough(self):
  genes=['CD3D','CD3E','CD2','LCK'];f,a,m=assay_aware_flags(np.ones((5,len(genes))),genes);self.assertFalse(a['T']);self.assertFalse(f[:,0].any())
 def test_different_cells_cannot_supply_modules(self):
  genes=['CD3D','CD3E','CD2','LCK','LAT'];f,a,m=assay_aware_flags(np.eye(5),genes);self.assertTrue(a['T']);self.assertFalse(f[:,0].any())
 def test_gamma_chain_not_created(self):
  genes=['CD3D','CD3E','CD7','LCK','ZAP70'];f,a,m=assay_aware_flags(np.ones((5,len(genes))),genes);self.assertFalse(a['gamma_delta_T_RNA_candidate']);self.assertFalse(f[:,GATE_ORDER.index('gamma_delta_T_RNA_candidate')].any())
 def test_original_assessable_arrays_identical(self):
  genes=['CD3D','CD3E','TRAC','LCK','LAT','CD8A','CD8B','FOXP3','IL2RA'];x=np.random.default_rng(33).integers(0,3,(30,len(genes)));old,a=coherent_flags(x,genes);new,b,m=assay_aware_flags(x,genes);np.testing.assert_array_equal(old,new);self.assertEqual(a,b)
if __name__=='__main__':unittest.main()
