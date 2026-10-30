import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_state_readout_candidate_policy_v1 import candidate_status

class CandidateCases(unittest.TestCase):
    def setUp(self):self.row=dict(status='estimable',n_units=10,n_finite_splits=20,split_half_Spearman_median=.99,positive_split_fraction=1.)
    def test_one_lucky_finite_split_cannot_pass(self):
        self.row['n_finite_splits']=1
        self.assertFalse(candidate_status(self.row,20,10)[0])
    def test_one_positive_material_cannot_pass(self):self.assertFalse(candidate_status(self.row,20,1)[0])
    def test_boundary_measurement_candidate(self):
        self.row.update(n_finite_splits=16,split_half_Spearman_median=.6,positive_split_fraction=.8)
        self.assertTrue(candidate_status(self.row,20,4)[0])
    def test_constant_or_unavailable_is_not_negative_state(self):
        self.row.update(status='no_between_unit_rank_variation',n_finite_splits=0,split_half_Spearman_median=float('nan'))
        passed,reason=candidate_status(self.row,20,0)
        self.assertFalse(passed);self.assertEqual(reason,'no_between_unit_rank_variation')

if __name__=='__main__':unittest.main()
