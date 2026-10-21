import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def load(name):
    s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load('analyze_moldoveanu_nuclear_signal_sensitivity_v1')
h=load('analyze_melanoma_IMC_continuous_state_neighbors_v2')


class NuclearSensitivity(unittest.TestCase):
    def test_rank_regression_normal_equations_and_unit_scores(self):
        rng=np.random.default_rng(52);n=120;hist=rng.normal(size=n);dna=rng.normal(size=n)
        frame=pd.DataFrame({'Ki67':hist+dna+rng.normal(size=n),'Histone.H3':hist,
            'X191Ir':dna+rng.normal(size=n)*.1,'X193Ir':dna+rng.normal(size=n)*.1})
        scores,reason=m.nuclear_scores(frame,h.rank_z)
        self.assertIsNone(reason);self.assertEqual(scores.shape,(n,3))
        rank,_=h.rank_z(frame[['Histone.H3']].to_numpy())
        self.assertAlmostEqual(float(scores[:,1]@rank[:,0]),0.,places=8)
        self.assertAlmostEqual(float(scores[:,1]@scores[:,2]),0.,places=8)
        np.testing.assert_allclose(np.mean(scores**2,axis=0),1.,atol=1e-10)
        self.assertFalse(np.allclose(scores[:,0],scores[:,1]))
    def test_constant_rank_deficient_and_no_residual_rejected(self):
        n=40;x=np.arange(n,dtype=float)
        frame=pd.DataFrame({'Ki67':x,'Histone.H3':x,'X191Ir':x,'X193Ir':x})
        value,reason=m.nuclear_scores(frame,h.rank_z)
        self.assertIsNone(value);self.assertIn('rank_deficient',reason)
        rng=np.random.default_rng(3);frame['X191Ir']=rng.normal(size=n);frame['X193Ir']=rng.normal(size=n)
        value,reason=m.nuclear_scores(frame,h.rank_z)
        self.assertIsNone(value);self.assertEqual(reason,'no_residual_Ki67_variation')
        frame['X191Ir']=0.
        value,reason=m.nuclear_scores(frame,h.rank_z)
        self.assertIsNone(value);self.assertIn('constant',reason)
    def test_full_2000_positive_negative_and_identical_eligibility(self):
        x=np.arange(40,dtype=float);x=(x-x.mean())/np.sqrt(np.mean((x-x.mean())**2))
        z=np.where(np.arange(40)%2,1.,-1.);right=np.column_stack([x,x,z])
        entries=[]
        for _ in range(5):
            entries.append(dict(coords=np.column_stack([np.arange(121),np.zeros(121)]),
                tumor=np.arange(121)==120,positions=[np.arange(40),np.arange(40,80),np.arange(80,120)],
                scores=[right.copy(),-right.copy(),right.copy()],reasons=[None,None,None]))
        edges=np.array([(i,80+i) for i in range(40)]+[(40+i,80+i) for i in range(40)])
        result=m.statistics(entries,lambda coords:edges,lambda c,e:np.zeros(len(c),dtype=int))
        self.assertEqual(len(result['family']),12)
        for row in result['family']:
            self.assertEqual(row['common_source_materials'],5)
            self.assertLess(row['BH'],.01)
            self.assertAlmostEqual(row['effect'],1. if row['parent_relation']=='regulatory_CD8' else -1.)
        self.assertEqual(m.core.METRICS,['continuous_rank_covariance','binary_standardized_covariance','binary_double_positive_edge_fraction'])
        json.dumps(result,allow_nan=False)
    def test_missing_design_excludes_all_six_positions(self):
        entries=[dict(coords=np.zeros((2,2)),tumor=np.zeros(2,dtype=bool),positions=[[],[],[]],
            scores=[None,None,None],reasons=['nuclear_design_rank_deficient']*3)]
        result=m.statistics(entries,lambda c:None,lambda c,e:None)
        self.assertTrue(all(row['common_source_materials']==0 and row['P'] is None for row in result['family']))
        self.assertTrue(all(row['missing_reason_counts']['nuclear_design_rank_deficient']==1 for row in result['family']))


if __name__=='__main__':unittest.main()
