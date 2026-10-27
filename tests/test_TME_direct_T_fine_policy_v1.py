from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
 import numpy as np
 from TME_direct_T_fine_policy_v1 import direct_fine,TYPES
 AVAILABLE=True
except ImportError:AVAILABLE=False

@unittest.skipUnless(AVAILABLE,'Scientific production requires all tests with zero skips')
class DirectFineTests(unittest.TestCase):
 def case(self,own,ref,major='T',fine='T',native=1,review=0,positive=True):
  return direct_fine([major],[fine],np.array([own],np.uint16),np.array([ref],np.uint16),np.array([native],np.uint16),np.array([review],np.uint8),[positive])
 def test_Treg_child_does_not_become_two_disagreeing_types(self):self.assertEqual(self.case(3,3)['proposed_fine_code'][0],TYPES.index('Treg'))
 def test_CD4_and_CD8_direct_cues(self):
  self.assertEqual(self.case(2,2)['proposed_fine_code'][0],TYPES.index('CD4_T'));self.assertEqual(self.case(4,4)['proposed_fine_code'][0],TYPES.index('CD8_T'))
 def test_regulatory_RNA_alone_does_not_override_CD8_reference(self):self.assertEqual(self.case(1,4)['proposed_fine_code'][0],0)
 def test_reference_without_native_cue_is_coarse(self):self.assertEqual(self.case(0,1)['proposed_fine_code'][0],0)
 def test_competing_direct_fine_cues_are_not_forced(self):self.assertEqual(self.case(5,5)['fine_review_code'][0],4)
 def test_gamma_stays_RNA_candidate(self):self.assertEqual(self.case(8,8)['fine_review_code'][0],2)
 def test_existing_fine_identity_is_not_replaced(self):self.assertEqual(self.case(4,4,fine='Treg')['proposed_fine_code'][0],0)
 def test_NK_and_T_NK_are_not_reclassified(self):
  self.assertEqual(self.case(4,4,major='NK')['proposed_fine_code'][0],0);self.assertEqual(self.case(4,4,major='T_NK')['proposed_fine_code'][0],0)
 def test_native_dropout_is_not_unstated_new_native_support(self):self.assertEqual(self.case(1,1,native=0)['fine_review_code'][0],5)
 def test_accepted_TCR_coarse_review_retains_its_qualified_native_context(self):self.assertEqual(self.case(1,1,native=5,review=6)['proposed_fine_code'][0],1)
 def test_unresolved_and_zero_RNA_remain(self):
  self.assertEqual(self.case(1,1,major='unresolved',fine='unresolved')['proposed_fine_code'][0],0);self.assertEqual(self.case(1,1,positive=False)['fine_review_code'][0],6)
 def test_shape_or_negative_masks_are_rejected(self):
  with self.assertRaises(ValueError):direct_fine(['T'],['T'],np.array([-1]),np.array([1]),np.array([1]),np.array([0]),[True])

if __name__=='__main__':unittest.main()
