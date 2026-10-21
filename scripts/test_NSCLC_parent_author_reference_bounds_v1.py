"""Truth checks use actual Boolean cells, not formulas duplicated in tests."""
import itertools
import unittest
import numpy as np
from NSCLC_parent_author_reference_bounds_v1 import regulatory_bounds,reference_compatibility


class Tests(unittest.TestCase):
    def test_all_single_cell_boolean_truth_tables(self):
        a=np.array(list(itertools.product([0,1],repeat=6)),int)
        truth=(a[:,0]&a[:,1]&a[:,2:5].any(axis=1)&~a[:,5].astype(bool)).astype(int)
        lo,hi=regulatory_bounds(np.ones(len(a),int),*a.T)
        np.testing.assert_array_equal(lo,truth);np.testing.assert_array_equal(hi,truth)

    def test_aggregated_units_cover_actual_unknown_joint(self):
        rng=np.random.default_rng(91673)
        for probability in [.1,.4,.8,.95]:
            a=rng.random((30,40,6))<probability
            truth=a[:,:,0]&a[:,:,1]&a[:,:,2:5].any(axis=2)&~a[:,:,5]
            lo,hi=regulatory_bounds(np.full(30,40),*a.sum(axis=1).T)
            self.assertTrue((lo<=truth.sum(axis=1)).all());self.assertTrue((hi>=truth.sum(axis=1)).all())
            reference=np.arange(30)%2==0
            result=reference_compatibility(lo,hi,reference,int(truth.sum()))
            actual=int(truth[reference].sum())
            self.assertLessEqual(result['compatible_reference_count_lower'],actual)
            self.assertGreaterEqual(result['compatible_reference_count_upper'],actual)

    def test_incoherent_marginals_rejected(self):
        with self.assertRaises(AssertionError):regulatory_bounds([1],[2],[0],[0],[0],[0],[0])
        with self.assertRaises(AssertionError):reference_compatibility([1],[1],[True],0)


if __name__=='__main__':unittest.main()
