import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from standard_study_holdout_readback_v1 import prototype_prediction,independent_scores
class HoldoutReadbackTests(unittest.TestCase):
    def test_inverse_distance_and_exact_center(self):
        centers=np.array([[0.],[10.]]);values=np.array([[0.],[1.]])
        pred=prototype_prediction(centers,values,np.array([[5.],[0.]]),2)
        self.assertAlmostEqual(pred[0,0],.5);self.assertLess(pred[1,0],1e-6)
    def test_positive_negative_depth_adjustment_and_unmeasured(self):
        n=100;t=np.arange(n);y=np.sin(t/5);total=np.exp(4+t/200)-1
        plus=independent_scores(y,y+np.log1p(total),total,30);minus=independent_scores(y,-y+np.log1p(total),total,30)
        self.assertAlmostEqual(plus['depth_adjusted_Pearson'],1.,places=10);self.assertAlmostEqual(minus['depth_adjusted_Pearson'],-1.,places=10)
        missing=independent_scores(np.full(n,np.nan),y,total,30);self.assertEqual(missing,dict(status='insufficient_evaluable_query_cells',n_cells=0))
if __name__=='__main__':unittest.main()
