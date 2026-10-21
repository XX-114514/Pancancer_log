"""Critical assay gates/zero-parent semantics and independently parsed raw prefix."""
from pathlib import Path
import gzip
import tempfile
import unittest
import numpy as np
from NSCLC_observed_T_state_RAM_reader_v1 import native_eligibility,row_counts,group_quality,prefix_audit


class NativeInputChecks(unittest.TestCase):
    def test_QC_boundaries(self):
        total=np.array([500,499,500,500,0]);det=np.array([250,300,249,250,1000]);mt=np.array([100,0,0,101,0])
        np.testing.assert_array_equal(native_eligibility(total,det,mt),[True,False,False,False,False])

    def test_case_counts_keep_observed_zero(self):
        flags=np.array([[1,0],[1,1],[0,1]],bool);materials=np.array([0,0,2])
        np.testing.assert_array_equal(row_counts(flags,materials,4),[[2,1],[0,0],[0,1],[0,0]])

    def test_absent_parent_not_filled_with_zero_quality(self):
        count,q=group_quality(np.array([True,False]),np.array([0,1]),np.array([500,800]),np.array([250,400]),np.array([5,3,0]),3)
        np.testing.assert_array_equal(count,[1,0,0]);self.assertEqual(q[0][3],.2)
        self.assertIsNone(q[1][0]);self.assertIsNone(q[2][3])

    def test_independent_raw_prefix_detects_count_mismatch(self):
        text='%%MatrixMarket matrix coordinate real general\n3 4 5\n1 1 2e+00\n1 3 4.0\n2 2 7\n2 4 3\n3 1 1\n'
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'synthetic_counts.mtx.gz'
            with gzip.open(path,'wt') as f:f.write(text)
            expected=np.array([[2,4],[0,0],[1,0]])
            p=prefix_audit(path,np.arange(4),[0,2],expected,np.array([6,10,1]),np.array([2,2,1]),rows=2)
            self.assertTrue(p['all_prefix_markers_and_total_and_detection_exact'])
            wrong=expected.copy();wrong[0,0]=3
            with self.assertRaises(AssertionError):prefix_audit(path,np.arange(4),[0,2],wrong,np.array([6,10,1]),np.array([2,2,1]),rows=2)


if __name__=='__main__':unittest.main()
