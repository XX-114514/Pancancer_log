from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import unittest
import numpy as np
from lineage_support_contingency_v1 import contingency


class Checks(unittest.TestCase):
    def test_same_marginals_different_cell_overlap(self):
        first=np.array([1,1,0,0],bool)
        overlapping=contingency(first,first);segregated=contingency(first,~first)
        self.assertEqual(overlapping['n_first_supported'],segregated['n_first_supported'])
        self.assertEqual(overlapping['n_second_supported'],segregated['n_second_supported'])
        self.assertEqual(overlapping['n_both'],2);self.assertEqual(segregated['n_both'],0)
        self.assertEqual(segregated['n_first_only'],2);self.assertEqual(segregated['n_second_only'],2)

    def test_empty_and_absent_support_have_missing_conditional_denominator(self):
        out=contingency(np.array([],bool),np.array([],bool));self.assertIsNone(out['fraction_both_in_group'])
        out=contingency(np.zeros(3,bool),np.array([1,0,1],bool));self.assertIsNone(out['fraction_first_with_second'])
        self.assertEqual(out['fraction_second_with_first'],0.)

    def test_non_boolean_and_mismatched_rows_rejected(self):
        for a,b in [(np.array([1]),np.array([True])),(np.ones(2,bool),np.ones(3,bool))]:
            with self.assertRaises(ValueError):contingency(a,b)


if __name__=='__main__':unittest.main()
