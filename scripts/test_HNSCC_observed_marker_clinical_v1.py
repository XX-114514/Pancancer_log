"""Check clinical endpoint direction, missing gates and conjunction semantics."""
import os
import unittest
import numpy as np
import GSE288199_state_pair_math_v1 as math
import GSE288199_state_pair_R_audit_v1 as reference
from run_HNSCC_observed_marker_clinical_v1 import arrays,select_rows,conjunction,clean


class ClinicalChecks(unittest.TestCase):
    def test_no_small_parent_or_phase_pooling(self):
        rows=[dict(phase=p,QC='human_RNA_QC500_250_MT20',metrics=[0,a,b])
              for p,a,b in [('pre',20,20),('pre',19,100),('post',100,100),('post2',100,100)]]
        allrows,use=select_rows(dict(rows=rows),'pre')
        self.assertEqual((len(allrows),len(use)),(2,1))

    def test_outcome_and_own_parent_quality(self):
        r=dict(metrics=[100,25,50,0,0,1,2,3,4,5,6,0,0],regression=60.,
               clinical=[0,1],quality=[[1,2,3,4],[5,6,7,8]])
        x,y,c=arrays([r],'E_memory')
        self.assertEqual(x[0],6/50)
        self.assertEqual(y[0],60)
        np.testing.assert_array_equal(c,[[0,1,5,6,7,8]])

    def test_conjunction_cannot_be_rescued_by_one_endpoint(self):
        self.assertEqual(conjunction(.001,.8),.8)
        self.assertEqual(conjunction(.001,1.),1.)

    def test_realistic_positive_negative_constant_and_R_replay(self):
        rng=np.random.default_rng(711)
        x=rng.normal(size=28);c=rng.normal(size=(28,6))
        blocks=[]
        for label,y in [('positive',x+.5*rng.normal(size=28)),('negative',-x+.5*rng.normal(size=28)),('constant',np.ones(28))]:
            m,f,status=math.association(x,y,c)
            if label=='constant': self.assertIsNone(f)
            else:
                self.assertEqual(np.sign(f['rho']),1 if label=='positive' else -1)
            signs=np.where(rng.random((99,28))<.5,-1.,1.) if f else None
            wild=math.wild_statistics(m,f,signs) if f else None
            indices=rng.integers(28,size=(30,28))
            boot,_=math.bootstrap(x,y,c,indices)
            blocks.append(clean(dict(label=label,x=x,y=y,C=c,expected=[f[k] for k in ['beta','se','t','rho']] if f else None,
                signs=signs,wild=wild,indices=indices,bootstrap=boot)))
        self.assertTrue(reference.audit(blocks)['status'].startswith('PASS_ALL'))


if __name__=='__main__': unittest.main()
