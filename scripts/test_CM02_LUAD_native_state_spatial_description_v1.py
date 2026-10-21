"""Finite descriptive edge opportunity and missingness checks, no human data."""
import unittest
import numpy as np
from CM02_LUAD_native_state_spatial_description_v1 import edge_effect, summarize, phenotype
import observed_T_marker_evidence_v1 as evidence


class Tests(unittest.TestCase):
    def fixture(self):
        parent=np.arange(40)<20;state=np.arange(40)<10
        receiver=~parent;positive=receiver & (np.arange(40)<30)
        source=np.repeat(np.arange(20),10)
        target=np.r_[np.tile(np.arange(20,30),10),np.tile(np.arange(30,40),10)]
        return parent,state,receiver,positive,source,target

    def test_positive_and_negative_edge_coupling(self):
        p,s,r,y,a,b=self.fixture()
        value,status=edge_effect(p,s,r,y,a,b)
        negative,_=edge_effect(p,s,r,r & ~y,a,b)
        self.assertGreater(value,0);self.assertAlmostEqual(value,-negative)
        self.assertIn('NO_CELL_LEVEL_P',status)

    def test_missing_capture_not_negative(self):
        p,s,r,y,a,b=self.fixture();y[:]=False
        value,status=edge_effect(p,s,r,y,a,b)
        self.assertIsNone(value);self.assertEqual(status,'STATE_OR_CONTRAST_CAPTURE_BELOW10')

    def test_zero_coupling_and_constant_images(self):
        p,s,r,y,a,b=self.fixture();b=np.tile(np.arange(20,40),10);a=np.repeat(np.arange(20),10)
        value,_=edge_effect(p,s,r,y,a,b)
        self.assertEqual(value,0)
        self.assertEqual(summarize([0,0])['median'],0)
        self.assertIsNone(summarize([])['median'])

    def test_full_scalar_parent_and_state_table(self):
        a=np.zeros((4,len(evidence.GENES)),np.int32);idx={g:i for i,g in enumerate(evidence.GENES)}
        for g in ['CD3E','FOXP3','IL2RA','LAG3']:a[0,idx[g]]=1
        for g in ['CD3E','CD8A','MKI67','TOP2A']:a[1,idx[g]]=1
        for g in ['CD3E','FOXP3','IL2RA','CD8A']:a[2,idx[g]]=1
        f,o,z=phenotype(a,evidence)
        self.assertTrue(f[0,6]);self.assertTrue(f[1,9]);self.assertTrue(f[2,3])
        self.assertFalse(o[2]);self.assertEqual(z.sum(),0)


if __name__=='__main__':unittest.main()
