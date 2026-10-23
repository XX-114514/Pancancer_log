"""Mutation counterexamples for the independent partition audit."""
from pathlib import Path
import sys,collections,unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import numpy as np
    from audit_annotation_discovery_v1 import reconcile_partition
    AVAILABLE=True
except ImportError:
    AVAILABLE=False


@unittest.skipUnless(AVAILABLE, 'Run in qualified scientific environment')
class AuditTests(unittest.TestCase):
    def setUp(self):
        self.source=np.array(['a','a','b','b','c'])
        self.route=np.array(['T','T','T','T','unknown'])
        self.fine=np.array(['x','y','x','y','z'])
        self.modeled=np.array([True,True,True,True,False])

    def tables(self, labels):
        return [{'partition_code':int(k),'n_cells':int((labels==k).sum()),
                 **{name:dict(collections.Counter(values[labels==k])) for name,values in [
                     ('source_counts',self.source),('inherited_route_counts',self.route),
                     ('inherited_fine_counts_not_truth',self.fine)]}} for k in np.unique(labels)]

    def check(self, labels, method='A', minimum=2, table=None):
        reconcile_partition(labels,self.source,self.route,self.fine,self.tables(labels) if table is None else table,
                            method,self.modeled,minimum)

    def test_cross_source_merge_and_old_fine_axis_mutation_fail(self):
        correct=np.array([0,0,1,1,-1]);self.check(correct)
        with self.assertRaises(ValueError):self.check(np.array([0,0,0,0,-1]))
        table=self.tables(correct);table[1]['inherited_fine_counts_not_truth']={'wrong':2}
        with self.assertRaises(ValueError):self.check(correct,table=table)

    def test_zero_RNA_and_small_group_unestimability_not_fabricated(self):
        self.check(np.array([-1,-1,-1,-1,-1]),minimum=3)
        with self.assertRaises(ValueError):self.check(np.array([0,0,1,1,2]))
        # C may pool across sources within a route; A cannot.
        self.check(np.array([0,0,0,0,-1]),method='C',minimum=3)


if __name__=='__main__':unittest.main(verbosity=2)
