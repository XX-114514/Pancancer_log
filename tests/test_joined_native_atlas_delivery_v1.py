import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from joined_native_atlas_delivery_v1 import domain_lookup,positions_checked

class Joins(unittest.TestCase):
    def test_sparse_domain_duplicates_order_and_missing(self):
        i,a=domain_lookup(np.array([1,4,9]),np.array([9,0,4,9]),10)
        np.testing.assert_array_equal(a,[True,False,True,True])
        np.testing.assert_array_equal(i[a],[2,1,2])
    def test_empty_domain_not_positive(self):
        i,a=domain_lookup(np.array([],dtype=int),np.array([0,2]),3)
        self.assertFalse(a.any())
    def test_domain_must_have_one_explicit_axis(self):
        for q in [np.array([1,1]),np.array([2,1]),np.array([-1]),np.array([3])]:
            with self.assertRaises(ValueError):domain_lookup(q,np.array([0]),3)
    def test_input_indices_are_not_silently_truncated(self):
        for p in [np.array([1.2]),np.array([True]),np.zeros((1,1)),np.arange(4097)]:
            with self.assertRaises(ValueError):positions_checked(p,10000)
        for p in [[-1],[3]]:
            with self.assertRaises(IndexError):positions_checked(p,3)

if __name__=='__main__':unittest.main()
