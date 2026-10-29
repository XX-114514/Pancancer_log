"""Do not certify a missing-assay read using an original-mode T case."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from deliver_TME_hierarchical_cluster_view_v1 import eligible_alternate_T_case

class DeliveryRoleTests(unittest.TestCase):
 def node(self,mode=None,major='T'):
  n={'cluster_RNA_proposal':{'major':major}}
  if mode is not None:n['assay_aware_T_evidence_mode']=mode
  return n
 def test_nonempty_original_mode_cannot_stand_in_for_alternative(self):
  self.assertFalse(eligible_alternate_T_case(self.node('original_native_module_assay_coverage_in_proposal')))
 def test_exact_alternate_mode_and_T_is_the_required_case(self):
  self.assertTrue(eligible_alternate_T_case(self.node('T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence')))
 def test_missing_grade_or_alternate_non_T_is_not_the_T_case(self):
  self.assertFalse(eligible_alternate_T_case(self.node()))
  self.assertFalse(eligible_alternate_T_case(self.node('T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence','NK')))
if __name__=='__main__':unittest.main()
