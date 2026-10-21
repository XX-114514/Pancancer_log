import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd


class ReversibleKeys(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=Path(__file__).resolve().parents[1]/'scripts/normalize_melanoma_IMC_ROI_embedded_cell_keys_v1.py'
        spec=importlib.util.spec_from_file_location('reversible_key_test',p)
        cls.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.m)

    def fixture(self):
        frame=pd.DataFrame(dict(sample_id=['a','a','b'],X=['1','2','1'],**{'Unnamed: 0':['11','12','13']},celltype20=['R','E','C']))
        nodes=pd.DataFrame({'':['0','1','2'],'Imageid':['b','a','a'],'type':['C','R','E'],'cell_id':['b_1','raw_1','raw_2']})
        return frame,nodes,{'raw':'a'}

    def test_complete_same_ROI_anchor_and_exact_labels(self):
        f,n,a=self.fixture();target,locations,formats=self.m.bind_embedded_nodes(f,n,a)
        self.assertEqual(locations.tolist(),[2,0,1]);self.assertTrue(formats)
        self.assertEqual(target[1],('a','raw_1'))

    def test_mismatched_label_rejected(self):
        f,n,a=self.fixture();n.loc[0,'type']='R'
        with self.assertRaises(AssertionError):self.m.bind_embedded_nodes(f,n,a)

    def test_unknown_ROI_anchor_and_number_offset_rejected(self):
        for replacement in ['unrelated_1','b_2']:
            f,n,a=self.fixture();n.loc[0,'cell_id']=replacement
            with self.assertRaises(AssertionError):self.m.bind_embedded_nodes(f,n,a)

    def test_duplicate_and_partial_identity_rejected(self):
        f,n,a=self.fixture();n.loc[2,'cell_id']='raw_1'
        with self.assertRaises(AssertionError):self.m.bind_embedded_nodes(f,n,a)
        f,n,a=self.fixture();n.loc[2,'cell_id']='raw_999'
        with self.assertRaises(AssertionError):self.m.bind_embedded_nodes(f,n,a)

    def test_suffix_format_is_reversible(self):
        f,n,a=self.fixture();n['cell_id']=['1:b','1:raw','2:raw']
        _,locations,_=self.m.bind_embedded_nodes(f,n,a)
        self.assertEqual(locations.tolist(),[2,0,1])

    def test_leading_zero_token_is_preserved(self):
        result=self.m.embedded_variants(['a_001'],['a'],{})
        self.assertIn(['001'],result.values());self.assertNotIn(['1'],result.values())

    def test_all_rows_must_follow_the_format_not_just_prefix_probes(self):
        values=['a_1']*1500;values[700]='wrong_1'
        result=self.m.embedded_variants(values,['a']*1500,{})
        self.assertEqual(list(result),['literal'])


if __name__=='__main__':unittest.main()
