import importlib.util
from pathlib import Path
import unittest
import pandas as pd


class CounterBridge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=Path(__file__).resolve().parents[1]/'scripts/melanoma_IMC_explicit_counter_bridge_v1.py'
        spec=importlib.util.spec_from_file_location('counter_bridge_test',p);cls.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.m)
    def fixture(self):
        f=pd.DataFrame({'Unnamed: 0':['3','1','2'],'sample_id':['b','a','a'],'celltype20':['C','R','E']})
        n=pd.DataFrame({'':['1','2','0'],'Imageid':['raw','b','raw'],'type':['E','C','R'],'cell_id':['cell1','cell2','cell0']})
        return f,n,{'raw':'a'}
    def test_shuffled_rows_and_original_author_endpoints(self):
        f,n,a=self.fixture();t,p,q=self.m.bind(f,n,a)
        self.assertEqual(p.tolist(),[2,0,1]);self.assertEqual(t[0],('a','cell1'));self.assertTrue(q['all_regions_and_types_exact'])
    def test_no_arbitrary_base(self):
        with self.assertRaises(AssertionError):self.m.counter(['4','5','6'])
    def test_no_missing_duplicate_or_formatted_tokens(self):
        for values in [['1','3'],['1','1'],['01','2'],['1.0','2.0']]:
            with self.assertRaises((AssertionError,ValueError)):self.m.counter(values)
    def test_label_conflict_fails(self):
        f,n,a=self.fixture();n.loc[0,'type']='R'
        with self.assertRaises(AssertionError):self.m.bind(f,n,a)
    def test_unknown_region_fails(self):
        f,n,a=self.fixture();n.loc[0,'Imageid']='unknown'
        with self.assertRaises(AssertionError):self.m.bind(f,n,a)
    def test_cell_token_must_be_exact_exported_counter(self):
        f,n,a=self.fixture();n.loc[0,'cell_id']='cell2'
        with self.assertRaises(AssertionError):self.m.bind(f,n,a)
    def test_whole_counter_still_requires_expected_bases(self):
        f,n,a=self.fixture();f['Unnamed: 0']=['2','0','1']
        with self.assertRaises(AssertionError):self.m.bind(f,n,a)


if __name__=='__main__':unittest.main()
