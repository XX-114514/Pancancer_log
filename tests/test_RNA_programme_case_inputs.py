import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

spec=importlib.util.spec_from_file_location('case_inputs',Path(__file__).resolve().parents[1]/'scripts/RNA_programme_case_inputs_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class CaseInputs(unittest.TestCase):
    def fixture(self):
        tab=pd.DataFrame([dict(library_id='synthetic_library',patient_registry='synthetic_case',timepoint_registry='pre',neo_type='D',**{'HPV/p16':'positive'})])
        cur=pd.DataFrame(dict(position=np.arange(100),library_id='synthetic_library',timepoint='pre',major='T',node_id=['member']*60+['outside']*40))
        unique=cur.iloc[:40].copy();unique['label']=['LAG3+ Treg']*20+['Tdys']*20
        score=np.column_stack([np.arange(40)/100,np.arange(40)/50]);total=np.arange(40)+10;detected=np.arange(40)+2
        return tab,cur,unique,score,total,detected

    def test_all_T_denominator_includes_unmatched_reference_and_own_depth_means(self):
        tab,cur,unique,scores,total,detected=self.fixture()
        units,coverage=m.aggregate(tab,cur,unique,{'member'},scores,total,detected)
        whole=units[m.CONTEXTS[0],'pre','synthetic_case',m.PARENTS[0]]
        member=units[m.CONTEXTS[1],'pre','synthetic_case',m.PARENTS[0]]
        self.assertEqual(whole['fraction'],.2);self.assertAlmostEqual(member['fraction'],1/3)
        self.assertAlmostEqual(whole['mean_log1p_UMI'],np.log1p(total[:20]).mean())
        self.assertNotAlmostEqual(whole['mean_log1p_UMI'],np.log1p(total[:20].mean()))
        self.assertEqual(len(coverage),12)
        slot=dict(context=m.CONTEXTS[0],phase='pre',left_parent=m.PARENTS[0],right_parent=m.PARENTS[2],left_programme='s.genes',right_programme='g2m.genes')
        x,y,c,cat,cases=m.pair_arrays(units,slot)
        self.assertAlmostEqual(x[0],scores[:20,0].mean());self.assertAlmostEqual(y[0],scores[20:,1].mean())
        self.assertEqual(c.shape,(1,8));self.assertEqual(len(cases),1)

    def test_absent_reference_unknown_background_phase_and_duplicate_not_reassigned(self):
        tab,cur,unique,scores,total,detected=self.fixture();unique.loc[unique.index[0],'label']='undeclared'
        with self.assertRaises(ValueError):m.aggregate(tab,cur,unique,{'member'},scores,total,detected)
        tab,cur,unique,scores,total,detected=self.fixture();tab['timepoint_registry']='post'
        with self.assertRaises(AssertionError):m.aggregate(tab,cur,unique,{'member'},scores,total,detected)
        tab,cur,unique,scores,total,detected=self.fixture()
        units,coverage=m.aggregate(tab,cur,unique,{'outside'},scores,total,detected)
        self.assertFalse(any(k[0]==m.CONTEXTS[1] for k in units))
        tab['neo_type']='unknown';units,coverage=m.aggregate(tab,cur,unique,{'member'},scores,total,detected)
        slot=dict(context=m.CONTEXTS[0],phase='pre',left_parent=m.PARENTS[0],right_parent=m.PARENTS[2],left_programme='s.genes',right_programme='s.genes')
        self.assertEqual(len(m.pair_arrays(units,slot)[0]),0)
        unique['position']=0
        with self.assertRaises(AssertionError):m.aggregate(tab,cur,unique,{'member'},scores,total,detected)


if __name__=='__main__':unittest.main()
