"""Separate populations are not the same observation as coexpressed cores."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from measure_TME_unresolved_core_overlap_v1 import pair_pattern
class CoreOverlapTests(unittest.TestCase):
 def test_equal_margins_can_have_opposite_same_cell_patterns(self):
  self.assertEqual(pair_pattern(100,50,50,0),'separated_core_populations');self.assertEqual(pair_pattern(100,50,50,50),'same_cell_core_coexpression')
 def test_small_disjoint_populations_do_not_meet_frozen_size(self):self.assertEqual(pair_pattern(20,10,10,0),'separated_small_populations')
 def test_intermediate_overlap_is_not_forced_into_either_extreme(self):self.assertEqual(pair_pattern(100,60,30,10),'intermediate_core_overlap')
 def test_impossible_joint_counts_are_rejected(self):
  with self.assertRaises(ValueError):pair_pattern(100,80,80,0)
if __name__=='__main__':unittest.main()
