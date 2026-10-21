"""Case-unit identity, sign direction and fixed multiplicity safeguards."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from CM02_marked_spatial_case_inference_v2 import case_signflip_p,cohort_summary,bh_fixed_family


class Tests(unittest.TestCase):
    def test_sign_and_case_order_invariance(self):
        d=np.linspace(.2,.5,25);keys=[f'synthetic_{i}' for i in range(25)]
        p=case_signflip_p(d,keys);self.assertLess(p,.01)
        self.assertEqual(p,case_signflip_p(-d,keys))
        self.assertEqual(p,case_signflip_p(d[::-1],keys[::-1]))
    def test_duplicate_case_block_rejected(self):
        with self.assertRaises(ValueError):case_signflip_p([.1,.2],['synthetic_same','synthetic_same'])
    def test_insufficient_cases_not_promoted(self):
        result=cohort_summary(np.ones(7),np.zeros((7,2000)),list(range(7)))
        self.assertNotIn('patient_signflip_P',result);self.assertIsNone(result['CI95'])
    def test_fixed_BH_and_conservative_BY(self):
        x={'a':{'patient_signflip_P':.01},'b':{'patient_signflip_P':.04},'missing':{}}
        bh_fixed_family(x,12);self.assertAlmostEqual(x['a']['BH_fixed_family'],.12)
        self.assertAlmostEqual(x['b']['BH_fixed_family'],.24)
        self.assertGreater(x['a']['BY_arbitrary_dependence_sensitivity'],.12)
        self.assertNotIn('BH_fixed_family',x['missing'])


if __name__=='__main__':unittest.main()
