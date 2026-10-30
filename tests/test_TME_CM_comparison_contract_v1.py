import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from TME_CM_comparison_contract_v1 import frequency_with_measurement_mask,partition_transport,effective_basis,map_module_mass,match_modules,weighted_jaccard

class CMContractTests(unittest.TestCase):
    def test_missing_parent_is_not_observed_zero_and_removed_nodes_keep_denominator(self):
        f,m=frequency_with_measurement_mask(np.array([[0,0,2]]),np.array([0]),np.array([[0,20,40]]),20)
        self.assertTrue(np.isnan(f[0,0]));self.assertEqual(f[0,1],0);self.assertEqual(f[0,2],.05)
        self.assertEqual(m.tolist(),[[False,True,True]])

    def test_split_node_uses_members_not_matching_names_and_test_cells_do_not_train_transport(self):
        a=np.array([0,0,1,1,1]);b=np.array([0,1,2,2,0]);train=np.array([True,True,True,True,False]);T,available,coverage=partition_transport(a,b,2,3,train)
        np.testing.assert_array_equal(T,np.array([[1,1,0],[0,0,1]]));self.assertEqual(coverage['n_training_atoms'],4)
        mass=map_module_mass(np.array([[.3],[.3],[.4]]),T);np.testing.assert_allclose(mass[:,0],[.6,.4])
        self.assertEqual(weighted_jaccard(mass[:,0],np.array([.6,.4])),1)

    def test_unmapped_loaded_node_cannot_vanish_silently(self):
        with self.assertRaises(ValueError):map_module_mass(np.array([[1.],[1.]]),np.array([[1.,0.]]))

    def test_assignment_does_not_reuse_a_single_new_module_and_leaves_orthogonal_module_unmatched(self):
        A=np.eye(3);B=np.array([[0,1],[1,0],[0,0]],float);r=match_modules(A,B,.5)
        self.assertEqual([(i,j) for i,j,s in r['matches']],[(0,1),(1,0)]);self.assertEqual(r['unmatched_A'],[2]);self.assertEqual(len({j for i,j,s in r['matches']}),2)

    def test_nsNMF_uses_effective_WS_and_negative_state_residuals_are_rejected(self):
        W=np.eye(2);S=np.array([[.75,.25],[.25,.75]]);np.testing.assert_array_equal(effective_basis(W,S),S)
        with self.assertRaises(ValueError):effective_basis(np.array([[1.,-1.]]))

    def test_overlapping_state_labels_are_not_mutually_exclusive_subgroup_counts(self):
        with self.assertRaises(ValueError):frequency_with_measurement_mask(np.array([[14],[14]]),np.array([0,0]),np.array([[20]]),20)

    def test_display_name_only_control_keeps_numeric_matrix_identical(self):
        counts=np.array([[2,0],[18,20]]);parent=np.array([0,0]);denominator=np.array([[20,20]])
        old=['old_marker_1','old_marker_2'];new=['new_name_1','new_name_2'];self.assertNotEqual(old,new)
        np.testing.assert_array_equal(frequency_with_measurement_mask(counts,parent,denominator,20)[0],np.array([[.1,0],[.9,1]]))

if __name__=='__main__':unittest.main()
