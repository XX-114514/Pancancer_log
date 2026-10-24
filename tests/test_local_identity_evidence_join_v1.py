from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import unittest
import numpy as np
from assemble_full_source_local_identity_review_v1 import category_counts,grouped_sums


class Checks(unittest.TestCase):
    def test_missing_label_counts_and_empty_reference_scope(self):
        out=category_counts(['T','__MISSING__','T','B'],np.array([0,0,1,1]),3)
        self.assertEqual(out,[{'T':1,'__MISSING__':1},{'B':1,'T':1},{}])
        self.assertEqual(category_counts([],np.array([],int),3),[{},{},{}])

    def test_group_evidence_sum_is_not_majority_label(self):
        g=np.array([0,0,1,1,1]);a=np.array([1,1,1,0,0],bool);b=np.array([0,1,1,1,0],bool)
        np.testing.assert_array_equal(grouped_sums(g,a,2),[2,1])
        np.testing.assert_array_equal(grouped_sums(g,b,2),[1,2])

    def test_invalid_group_axes_rejected(self):
        for values,g in [(['T'],np.array([-1])),(['T'],np.array([2])),(['T','B'],np.array([0]))]:
            with self.assertRaises(ValueError):category_counts(values,g,2)
        with self.assertRaises(ValueError):grouped_sums(np.array([0]),np.array([1,2]),2)


if __name__=='__main__':unittest.main()
