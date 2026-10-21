"""Synthetic cell-level oracle for subtype marginals and quality."""
import unittest
import numpy as np
import pandas as pd
from NSCLC_subtype_marker_measurement_RAM_v1 import summarize_units


class Measurement(unittest.TestCase):
    def test_actual_cell_oracle_and_absent_subtypes(self):
        rng=np.random.default_rng(34873);labels=['r0','r1','r2','o0','o1']
        definitions={'Treg':labels[:3],'CD4T_nonTreg':labels[3:]};names=['a','b','c']
        registry=[];positives=[];expected={}
        for case in names[:2]:
            for label in labels:
                n=int(rng.integers(1,20));x=rng.poisson(.7,(n,3));rna=rng.uniform(500,9000,n);genes=rng.uniform(250,2500,n)
                registry.append(dict(sampleID=case,sub_cell_type=label,n_cells=n,mean_logRNA=np.log1p(rna).mean(),mean_loggenes=np.log1p(genes).mean()))
                positives.append((x>0).sum(axis=0));expected[case,label]=(n,(x>0).sum(axis=0),[np.log1p(rna).mean(),np.log1p(genes).mean()])
        out,_=summarize_units(pd.DataFrame(registry),np.asarray(positives),range(3),names,definitions)
        for label in labels:
            for i,case in enumerate(names[:2]):
                n,m,q=expected[case,label];self.assertEqual(out[label]['count'][i],n)
                np.testing.assert_array_equal(out[label]['markers'][i],m);np.testing.assert_allclose(out[label]['means'][i],q)
            self.assertEqual(out[label]['count'][2],0);np.testing.assert_array_equal(out[label]['markers'][2],0)

    def test_duplicate_source_units_rejected(self):
        r=dict(sampleID='a',sub_cell_type='r',n_cells=1,mean_logRNA=1.,mean_loggenes=1.)
        with self.assertRaises(AssertionError):summarize_units(pd.DataFrame([r,r]),np.ones((2,3),int),range(3),['a'],{'Treg':['r'],'CD4T_nonTreg':['o']})

    def test_impossible_positive_counts_rejected(self):
        r=dict(sampleID='a',sub_cell_type='r',n_cells=1,mean_logRNA=1.,mean_loggenes=1.)
        with self.assertRaises(AssertionError):summarize_units(pd.DataFrame([r]),np.full((1,3),2),range(3),['a'],{'Treg':['r'],'CD4T_nonTreg':['o']})


if __name__=='__main__':unittest.main()
