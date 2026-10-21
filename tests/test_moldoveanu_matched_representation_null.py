"""Matched denominator, full permutations and preserved tumor affinity."""
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
m=load('analyze_moldoveanu_matched_representation_null_v1')
h=load('analyze_melanoma_IMC_continuous_state_neighbors_v2')


def fixture(affinity=False):
    rows=[]
    for material in range(5):
        for role,offset in [('Treg',0),('Th.ae',40),('Tc.ae',80)]:
            for i in range(40):
                positive=i>=20 if affinity else i%2==1
                if role=='Th.ae':positive=not positive
                rows.append({'sample.id':f'synthetic_{material}','Cluster':role,
                    'coord.x':float(offset+i),'coord.y':0.,'Ki67':float(39-i if role=='Th.ae' else i),
                    'Ki67.prob':float(positive),'Histone.H3':float(i),'X190BCKG':0.,'X191Ir':0.,'X193Ir':0.})
        rows.append({'sample.id':f'synthetic_{material}','Cluster':'melano','coord.x':120.,'coord.y':0.,
                     'Ki67':0.,'Ki67.prob':0.,'Histone.H3':0.,'X190BCKG':0.,'X191Ir':0.,'X193Ir':0.})
    edges=np.array([(i,80+i) for i in range(40)]+[(40+i,80+i) for i in range(40)]+
                   [(offset+i,120) for offset in (0,40,80) for i in range(20)])
    return pd.DataFrame(rows),edges


class MatchedNull(unittest.TestCase):
    def run_stats(self,frame,edges):
        return m.statistics(frame,lambda coords:edges,h.rank_z,
                            lambda coords,extent:np.zeros(len(coords),dtype=int))
    def test_full_2000_positive_negative_common_denominator_and_joint_scores(self):
        frame,edges=fixture();ordinary=m.three_statistics;checks=[]
        def checked(x,y,bx,by,a,b):
            np.testing.assert_allclose(x[:,1],2*bx-1,atol=1e-12)
            np.testing.assert_allclose(y[:,1],2*by-1,atol=1e-12)
            checks.append(True);return ordinary(x,y,bx,by,a,b)
        m.three_statistics=checked
        try:result=self.run_stats(frame,edges)
        finally:m.three_statistics=ordinary
        rows=result['family'];self.assertEqual(len(rows),12)
        for row in rows:
            self.assertEqual(row['common_source_materials'],5)
            self.assertIsNone(row['patient_CI'])
            self.assertLess(row['BH'],.01)
            self.assertEqual(row['effect']>0,row['parent_relation']=='regulatory_CD8')
        self.assertAlmostEqual(rows[0]['effect'],1.)
        self.assertAlmostEqual(rows[1]['effect'],1.)
        self.assertAlmostEqual(rows[2]['effect'],.25)
        self.assertNotIn('synthetic_',json.dumps(result,allow_nan=False))
        self.assertEqual(len(checks),5*(2+2*2+2*2000*2))
    def test_only_shared_tumor_affinity_is_removed_without_changing_materials(self):
        frame,edges=fixture(affinity=True);rows=self.run_stats(frame,edges)['family']
        for parent in range(2):
            grid,conditioned=rows[parent*6:parent*6+3],rows[parent*6+3:parent*6+6]
            self.assertTrue(all(row['common_source_materials']==5 for row in grid+conditioned))
            self.assertGreater(abs(grid[1]['effect']),.9)
            for row in conditioned[1:]:
                self.assertAlmostEqual(row['effect'],0.)
                self.assertEqual(row['P'],1.)
    def test_pure_qualification_names_constant_channels_and_threshold_ties(self):
        frame,_=fixture();frame.loc[0,'Ki67.prob']=.5
        result=m.profile_background(frame,h.rank_z)
        self.assertTrue(result['no_background_effects_or_clinical_outcomes_computed'])
        for role in result['roles']:
            profile={x['field']:x['counts'] for x in role['per_channel_source_material_counts']}
            self.assertEqual(profile['Histone.H3']['variable'],5)
            self.assertEqual(profile['X190BCKG']['constant'],5)
            self.assertEqual(role['joint_design_counts']['constant_or_nonfinite_required_background_no_covariate_dropped'],5)
        self.assertEqual(result['roles'][0]['threshold_qualification_counts']['source_cells_equal_threshold'],1)
        self.assertNotIn('synthetic_',json.dumps(result))
    def test_constant_binary_excludes_same_material_from_all_six_positions(self):
        frame,edges=fixture();take=frame['sample.id'].eq('synthetic_0') & frame['Cluster'].eq('Treg')
        frame.loc[take,'Ki67.prob']=0.
        result=self.run_stats(frame,edges);rows=result['family']
        self.assertTrue(all(row['common_source_materials']==4 and row['P'] is None for row in rows[:6]))
        self.assertTrue(all(row['common_source_materials']==5 for row in rows[6:]))
        for row in rows[:6]:self.assertEqual(row['missing_reason_counts']['no_binary_variation'],1)
    def test_permutation_indices_never_cross_grid_or_tumor_status(self):
        groups=np.array([0,0,1,1,16,16,17,17])
        indices=m.permute_indices(groups,np.random.default_rng(50))
        np.testing.assert_array_equal(groups[indices],groups)
        self.assertEqual(sorted(indices.tolist()),list(range(8)))


if __name__=='__main__':unittest.main()
