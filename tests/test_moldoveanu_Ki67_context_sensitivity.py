"""Check real inference wiring, context conditioning and rejected designs."""
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/f'{name}.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module
m = load('analyze_moldoveanu_Ki67_context_sensitivity_v1')
h = load('analyze_melanoma_IMC_continuous_state_neighbors_v2')


def fixture(affinity_only=False):
    rng = np.random.default_rng(831)
    rows = []
    for material in range(5):
        for role, offset in [('Treg', 0), ('Th.ae', 40), ('Tc.ae', 80)]:
            for i in range(40):
                positive = i < 20 if affinity_only else i % 2 == 1
                if role == 'Th.ae' and not affinity_only: positive = not positive
                row = {'sample.id':f'synthetic_{material}', 'Cluster':role,
                       'coord.x':float(offset+i), 'coord.y':0.,
                       'Ki67.prob':float(positive), 'Ki67':float(positive)+rng.uniform(0, .01)}
                row.update({field:float(rng.normal()) for field in m.BACKGROUND}); rows.append(row)
        rows.append({'sample.id':f'synthetic_{material}', 'Cluster':'melano',
                     'coord.x':120., 'coord.y':0., 'Ki67.prob':0., 'Ki67':0.,
                     **{field:0. for field in m.BACKGROUND}})
    edges = [(i,80+i) for i in range(40)] + [(40+i,80+i) for i in range(40)]
    edges += [(offset+i,120) for offset in (0,40,80) for i in range(20)]
    return pd.DataFrame(rows), np.array(edges)


class ContextSensitivity(unittest.TestCase):
    def test_full_2000_positive_negative_and_source_summary_only(self):
        frame, edges = fixture()
        result = m.statistics(frame, lambda coords:edges, h.rank_z,
                              lambda coords, extent:np.zeros(len(coords), dtype=int))
        self.assertEqual(len(result['family']), 8)
        for slot in (0,1,2):
            self.assertAlmostEqual(result['family'][slot]['effect'], .25)
            self.assertLess(result['family'][slot]['BH'], .01)
        for slot in (4,5,6):
            self.assertAlmostEqual(result['family'][slot]['effect'], -.25)
            self.assertLess(result['family'][slot]['BH'], .01)
        self.assertNotIn('synthetic_', json.dumps(result, allow_nan=False))
        for row in result['family']:
            self.assertEqual(row['source_materials'], 5)
            self.assertIsNone(row['patient_CI'])
        self.assertEqual(result['source_materials_considered'], 5)

    def test_tumor_affinity_preserved_not_false_pair_state_enrichment(self):
        frame, edges = fixture(affinity_only=True)
        result = m.statistics(frame, lambda coords:edges, h.rank_z,
                              lambda coords, extent:np.zeros(len(coords), dtype=int))
        for slot in (0,4):
            self.assertEqual(result['family'][slot]['effect'], 0.)
            self.assertEqual(result['family'][slot]['P'], 1.)
        for slot in (1,2,5,6):
            self.assertIsNone(result['family'][slot]['P'])
            self.assertEqual(result['family'][slot]['missing_reason_counts']['no_state_variation_in_context_parent'], 5)

    def test_background_collinearity_and_constant_not_dropped(self):
        frame, _ = fixture();part = frame.iloc[:40].copy()
        part[m.BACKGROUND[0]] = 1.
        result, reason = m.nuclear_residual(part, h.rank_z)
        self.assertIsNone(result);self.assertIn('constant', reason)
        part[m.BACKGROUND[0]] = part[m.BACKGROUND[1]]
        result, reason = m.nuclear_residual(part, h.rank_z)
        self.assertIsNone(result);self.assertIn('not_identifiable', reason)

    def test_missing_family_and_quadrant_identity(self):
        rows = m.aggregate([[] for _ in range(8)], [m.Counter() for _ in range(8)], [[],[]])['family']
        self.assertTrue(all(x['P'] is None and x['BH'] is None for x in rows))
        x = np.arange(20, dtype=float); a = np.arange(20)
        quadrants = m.quadrant_contributions(x, x, a, a)
        self.assertAlmostEqual(sum(quadrants.values()), m.covariance(x, x, a, a))
        self.assertGreater(quadrants['below_below'], 0.)
        self.assertGreater(quadrants['above_above'], 0.)


if __name__ == '__main__':unittest.main()
