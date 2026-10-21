"""Exact source key and geometry gates; no source patient records."""
import importlib.util
from pathlib import Path
import unittest
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mold_bridge',ROOT/'scripts/qualify_moldoveanu_ICI_cells_in_RAM_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class ExactCellBridge(unittest.TestCase):
    def fixture(self):
        left=pd.DataFrame({'obj.id':['a','b','outside'],'sample.id':['synthetic','synthetic','other'],
            'coord.x':[1.,2.,0.],'coord.y':[3.,4.,0.],'Cluster':['Treg','Tc.ae','Th.ae'],
            'Ki67.prob':[.2,.8,np.nan]})
        raw=pd.DataFrame({'Cell_ID':['a','b','excluded'],'Sample_ID':['synthetic']*3,
            'Location_Center_X':[1.,2.,5.],'Location_Center_Y':[3.,4.,6.],'Ki67':[.1,1.,0.]})
        return left,raw
    def test_exact_join_keeps_author_QC_and_nonICI_missing_probability(self):
        joined,q=m.bind(*self.fixture());self.assertEqual(len(joined),2)
        self.assertEqual(q['raw_unselected_cells_remain_excluded'],1)
        self.assertEqual(q['coordinate_maximum_absolute_difference'],0.)
    def test_duplicate_or_unmatched_source_key_rejected(self):
        left,raw=self.fixture()
        with self.assertRaises(AssertionError):m.bind(pd.concat([left,left.iloc[:1]]),raw)
        raw.loc[0,'Cell_ID']='wrong'
        with self.assertRaises(AssertionError):m.bind(left,raw)
    def test_coordinate_or_probability_corruption_rejected(self):
        left,raw=self.fixture();raw.loc[0,'Location_Center_X']+=.01
        with self.assertRaises(AssertionError):m.bind(left,raw)
        left,raw=self.fixture();left.loc[0,'Ki67.prob']=1.1
        with self.assertRaises(AssertionError):m.bind(left,raw)


if __name__=='__main__':unittest.main()
