"""Exact group quality and zero/absent parent semantics, no human data."""
import unittest
import numpy as np
from NSCLC_author_state_counts_RAM_v1 import aggregate


class Tests(unittest.TestCase):
    def test_exact_quality_not_log_of_mean(self):
        mask=np.array([True,True,False,True]);case=np.array([0,0,0,2]);RNA=np.array([9,99,999,49]);genes=np.array([4,24,99,9])
        n,q=aggregate(mask,case,RNA,genes,np.array([4,0,2]),3)
        np.testing.assert_array_equal(n,[2,0,1])
        self.assertAlmostEqual(q[0][0],(np.log(10)+np.log(100))/2)
        self.assertNotAlmostEqual(q[0][0],np.log1p(RNA[:2].mean()))
        self.assertEqual(q[0][3],.5);self.assertIsNone(q[1][0]);self.assertIsNone(q[1][3])

    def test_observed_parent_zero_state_is_not_missing_quality(self):
        n,q=aggregate(np.array([True,False]),np.array([0,0]),np.array([9,19]),np.array([3,4]),np.array([2]),1)
        self.assertEqual(n[0],1);self.assertTrue(np.isfinite(q[0][0]));self.assertEqual(q[0][3],.5)


if __name__=='__main__':unittest.main()
