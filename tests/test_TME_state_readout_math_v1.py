import argparse,json,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from validate_TME_state_readouts_v1 import column_spearman,split_means

class ReadoutMath(unittest.TestCase):
    def test_ties_positive_and_negative(self):
        a=np.array([[1,1],[1,1],[2,2],[3,3]],float)
        b=np.column_stack([np.array([5,5,6,7]),np.array([7,7,6,5])])
        np.testing.assert_allclose(column_spearman(a,b),[1,-1],atol=1e-12)
    def test_constant_is_unestimable_not_zero(self):
        a=np.array([[1,1],[2,1],[3,1],[4,1]],float)
        actual=column_spearman(a,np.array([4,3,2,1]))
        self.assertAlmostEqual(actual[0],-1);self.assertTrue(np.isnan(actual[1]))
    def test_split_preserves_even_and_odd_material_means(self):
        for n in [4,5,41]:
            y=np.column_stack([np.arange(n),np.arange(n)**2]).astype(float)
            a,b,na,nb=split_means(y,[20261031,n])
            self.assertEqual(na+nb,n);self.assertLessEqual(abs(na-nb),1)
            np.testing.assert_allclose((na*a+nb*b)/n,y.mean(axis=0),atol=1e-12)
            again=split_means(y,[20261031,n]);np.testing.assert_array_equal(a,again[0]);np.testing.assert_array_equal(b,again[1])
    def test_missing_observation_rejected_in_rank_correlation(self):
        with self.assertRaises(AssertionError):column_spearman(np.array([[1.],[np.nan]]),np.array([1.,2.]))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    result=unittest.main(argv=[sys.argv[0]],exit=False).result
    if args.output and result.wasSuccessful():
        r=Path(args.output);(r/'outputs/readout_math_acceptance_v1.json').write_text(json.dumps(dict(synthetic_only=True,n_tests=result.testsRun,passed=True)))
        (r/'provenance/terminal_v1.json').write_text(json.dumps(dict(status='STATE_READOUT_RANK_AND_MATERIAL_SPLIT_COUNTEREXAMPLES_PASS',synthetic_only=True)))
    sys.exit(0 if result.wasSuccessful() else 1)
