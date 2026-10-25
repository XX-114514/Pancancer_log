from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import unittest
import numpy as np
from common_reference_sampling_v1 import reserve_then_balance,stratum_reserve,select_reference


class Checks(unittest.TestCase):
    def test_small_source_and_local_rare_reserve(self):
        np.testing.assert_array_equal(reserve_then_balance([100,2,100],[2,2,2],10,7).sum(),10)
        q=reserve_then_balance([100,2,100],[2,2,2],10,7);self.assertTrue((q>=2).all());self.assertEqual(q[1],2)
        unit=np.array(['u1']*501+['u2']*500);strata=np.array(['rare']+['common']*1000);eligible=np.ones(1001,bool)
        out=select_reference(unit,strata,eligible,10,19)
        self.assertIn(0,out);self.assertEqual(len(out),10);self.assertEqual(len(set(out)),10)

    def test_no_silent_loss_when_budget_below_reserves(self):
        with self.assertRaises(ValueError):reserve_then_balance([100,100],[2,2],3,2)
        with self.assertRaises(ValueError):select_reference(['a']*5,['a','b','c','d','e'],np.ones(5,bool),3,2)
        self.assertEqual(stratum_reserve(['a','b','b'],np.array([1,1,0],bool)),2)

    def test_deterministic_eligibility_and_empty_scope(self):
        unit=np.array(['u1','u1','u2','u2','u2']);strata=np.array(['x','x','x','y','y']);eligible=np.array([1,0,1,1,1],bool)
        a=select_reference(unit,strata,eligible,3,18,minimum=1);b=select_reference(unit,strata,eligible,3,18,minimum=1)
        np.testing.assert_array_equal(a,b);self.assertNotIn(1,a)
        self.assertEqual(len(select_reference(unit,strata,np.zeros(5,bool),10,18)),0)


if __name__=='__main__':unittest.main()
