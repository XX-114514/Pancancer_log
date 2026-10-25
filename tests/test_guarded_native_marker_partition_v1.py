from pathlib import Path
import sys, unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from guarded_native_marker_partition_v1 import partition_guarded_subsets


class PartitionChecks(unittest.TestCase):
    def test_split_parent_and_same_type_different_parents(self):
        code,pairs=partition_guarded_subsets(np.array([0,0,0,1,1,-1]),np.array([1,2,1,1,1,2]),np.array([True]*5+[False]))
        np.testing.assert_array_equal(pairs,[[0,1],[0,2],[1,1]])
        np.testing.assert_array_equal(code,[0,1,0,2,2,-1])
        self.assertNotEqual(code[0],code[3])

    def test_permutation_and_all_excluded(self):
        parent=np.array([4,0,4,-1]);types=np.array([2,3,1,9]);fit=parent>=0
        code,pairs=partition_guarded_subsets(parent,types,fit)
        order=np.array([2,0,3,1]);permuted,other=partition_guarded_subsets(parent[order],types[order],fit[order])
        np.testing.assert_array_equal(other,pairs)
        np.testing.assert_array_equal(permuted,code[order])
        code,pairs=partition_guarded_subsets(np.array([-1]),np.array([3]),np.array([False]))
        np.testing.assert_array_equal(code,[-1]);self.assertEqual(pairs.shape,(0,2))

    def test_exclusion_mismatch_and_invalid_axes_rejected(self):
        for parent,types,fit in [(np.array([0]),np.array([2]),np.array([False])),(np.array([-1]),np.array([2]),np.array([True])),(np.array([0]),np.array([-1]),np.array([True])),(np.array([0]),np.array([1,2]),np.array([True]))]:
            with self.assertRaises(ValueError):partition_guarded_subsets(parent,types,fit)


if __name__=='__main__':unittest.main()
