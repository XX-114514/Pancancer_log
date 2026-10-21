import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def load(name,filename):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/filename);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


class NumericReplay(unittest.TestCase):
    def cores(self):
        a=load('plain','analyze_melanoma_IMC_author_graph_states_v3.py');b=load('audited','analyze_melanoma_IMC_author_graph_states_v3.py')
        auditor=load('audit','melanoma_IMC_numeric_replay_audit_v1.py');counts=auditor.attach(b);return a,b,counts
    def family(self):
        rng=np.random.default_rng(63);observed=rng.normal(size=(12,18));observed[:,0]+=2;observed[:,1]-=2
        expected=rng.normal(scale=.1,size=(12,18));draws=rng.normal(size=(12,2000,18));counts=np.ones((12,18),dtype=int);counts[:,17]=0;counts[5:,16]=0
        return dict(observed=observed,expected=expected,draws=draws,roi_counts=counts,family=[dict(slot=i) for i in range(18)],seed=77)
    def test_2000_covariance_stream_unchanged(self):
        a,b,c=self.cores();rng=np.random.default_rng(8);left=rng.normal(size=(20,3));right=rng.normal(size=(25,3));i=rng.integers(0,20,100);j=rng.integers(0,25,100)
        r1=np.random.default_rng(17);r2=np.random.default_rng(17)
        x=a.roi_permutations(left,right,i,j,np.zeros(20),np.zeros(25),2000,r1)
        y=b.roi_permutations(left,right,i,j,np.zeros(20),np.zeros(25),2000,r2)
        for u,v in zip(x,y):np.testing.assert_array_equal(u,v)
        self.assertEqual(r1.random(),r2.random());self.assertEqual(c['covariance_calls'],2002)
    def test_all_family_effects_P_BH_bootstrap_unchanged(self):
        a,b,c=self.cores();data=self.family();x=a.finish_family(**data);y=b.finish_family(**data)
        self.assertEqual(x,y);self.assertEqual(c['family_positions'],18);self.assertEqual(c['permutation_source_means'],34000);self.assertEqual(c['bootstrap_means'],32000)
        self.assertGreater(y[0]['effect'],0);self.assertLess(y[1]['effect'],0);self.assertIsNone(y[17]['effect'])
    def test_covariance_corruption_detected(self):
        core=load('bad_cov','analyze_melanoma_IMC_author_graph_states_v3.py');ordinary=core.edge_covariance;core.edge_covariance=lambda *args:ordinary(*args)+.1
        load('audit_bad','melanoma_IMC_numeric_replay_audit_v1.py').attach(core)
        with self.assertRaises(AssertionError):core.edge_covariance(np.arange(4.)[:,None],np.arange(4.)[:,None],np.arange(4),np.arange(4))
    def test_P_corruption_detected(self):
        core=load('bad_p','analyze_melanoma_IMC_author_graph_states_v3.py');ordinary=core.finish_family
        def corrupt(*args,**kwargs):
            rows=ordinary(*args,**kwargs);rows[0]['P']=.123456;return rows
        core.finish_family=corrupt;load('audit_bad_p','melanoma_IMC_numeric_replay_audit_v1.py').attach(core)
        with self.assertRaises(AssertionError):core.finish_family(**self.family())
    def test_CI_corruption_detected(self):
        core=load('bad_ci','analyze_melanoma_IMC_author_graph_states_v3.py');ordinary=core.finish_family
        def corrupt(*args,**kwargs):
            rows=ordinary(*args,**kwargs);rows[0]['case_bootstrap_95_CI']=[-100,100];return rows
        core.finish_family=corrupt;load('audit_bad_ci','melanoma_IMC_numeric_replay_audit_v1.py').attach(core)
        with self.assertRaises(AssertionError):core.finish_family(**self.family())


if __name__=='__main__':unittest.main()
