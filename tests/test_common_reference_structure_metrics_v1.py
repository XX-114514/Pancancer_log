import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from common_reference_structure_metrics_v1 import exact_neighbor_graph,paired_native_states


class StructureChecks(unittest.TestCase):
    def test_tied_coordinates_never_keep_self(self):
        g=exact_neighbor_graph(np.zeros((120,3)),30)
        self.assertTrue(np.all(g.diagonal()==0));self.assertTrue(np.all(np.diff(g.indptr)==30))

    def test_equal_coordinates_pair_and_independent_predictions(self):
        rng=np.random.default_rng(7);z=rng.normal(size=(160,5));u=np.tile(np.arange(8),20);d=rng.normal(size=160);y=(z[:,0]+.1*d)[:,None]
        graph=exact_neighbor_graph(z);metrics,pairs,error=paired_native_states({'a':z,'b':z.copy()},u,d,y,['native'],audit=True,graphs={'a':graph,'b':graph})
        self.assertLess(error,1e-9);self.assertEqual(pairs[0]['status'],'estimated');self.assertAlmostEqual(pairs[0]['second_minus_first'],0.)

    def test_single_unit_and_missing_readout_not_negative_evidence(self):
        rng=np.random.default_rng(8);z=rng.normal(size=(120,3));y=np.full((120,1),np.nan)
        r,p,e=paired_native_states({'a':z,'b':z.copy()},np.zeros(120),np.ones(120),y,['missing'],audit=True)
        self.assertTrue(all(x['n_evaluable']==0 for x in r));self.assertNotEqual(p[0]['status'],'estimated')


if __name__=='__main__':unittest.main()
