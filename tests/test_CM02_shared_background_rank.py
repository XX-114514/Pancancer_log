from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import CM02_shared_background_rank_v1 as m


class SharedBackground(unittest.TestCase):
    def fixture(self,shared_T=False):
        rng=np.random.default_rng(79);n=80;c=rng.normal(size=(n,8));cat=np.array([['D' if i%2 else 'D+T','positive' if i%3 else 'negative'] for i in range(n)])
        background=rng.normal(size=n);t=rng.normal(size=n);x=2*background+.3*rng.normal(size=n);y=2*background+.3*rng.normal(size=n)
        if shared_T:x+=3*t;y+=3*t
        return x,y,c,cat,background+.03*rng.normal(size=n),background+.03*rng.normal(size=n)

    def test_shared_background_can_attenuate_and_independent_paired_bootstrap_replays(self):
        args=self.fixture();a=m.analyze(*args,4,draws=150,required=120);b=m.analyze(*args,4,draws=150,required=120,independent=True);m.assert_replay(a,b)
        self.assertGreater(a['same_case_original_rho'],.7);self.assertLess(a['delta_rho'],-.5);self.assertLess(a['delta_bootstrap95'][1],0)
        self.assertEqual(set(a['point_sensitivities']),{'author_concordant_B_only','author_concordant_myeloid_only'})

    def test_shared_T_component_can_remain_after_background(self):
        a=m.analyze(*self.fixture(True),5,draws=100,required=80);self.assertGreater(a['with_joint_background_rho'],.6)

    def test_complete_controls_and_case_count_are_not_dropped(self):
        x,y,c,cat,b,mm=self.fixture();b[:]=1;mm[:]=1
        a=m.analyze(x,y,c,cat,b,mm,0,draws=10,required=8);self.assertEqual(a['reason'],'CONSTANT_CONTINUOUS_CONTROL')
        a=m.analyze(x[:19],y[:19],c[:19],cat[:19],b[:19],mm[:19],0);self.assertEqual(a['reason'],'INSUFFICIENT_ELIGIBLE_CASES')


if __name__=='__main__':unittest.main()
