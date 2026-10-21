import importlib.util
import ast
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import CM02_RNA_state_stability_helpers_v1 as m


class Stability(unittest.TestCase):
    def test_four_programme_native_reader_replays_closed_first_two_and_dense_halves(self):
        from types import SimpleNamespace
        import ucell_native_sparse_scores_v1 as scoring
        resource=SimpleNamespace(RUSAGE_SELF=0,getrusage=lambda unused:SimpleNamespace(ru_maxrss=1024))
        env=dict(np=np,pd=pd,resource=resource,scoring=scoring,baseline=SimpleNamespace(scoring=scoring))
        for file,name in [('analyze_CM02_RNA_multigene_in_RAM_v1.py','read_scores'),('characterize_CM02_G2M_state_stability_v1.py','read_all_scores')]:
            function=next(n for n in ast.parse((ROOT/'scripts'/file).read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
            exec(compile(ast.Module(body=[function],type_ignores=[]),'synthetic_reader','exec'),env)
        class Raw(dict):attrs={'shape':(5002,5000)}
        ptr=np.zeros(5003,dtype=np.int64);ptr[1:5001]=4;ptr[5001:]=7
        raw=Raw(indptr=ptr,indices=np.array([4,4,2,4999,0,2,14]),data=np.array([2,3,2,500,4,6,1]))
        z={'layers/counts':raw};human=np.ones(5000,bool);human[-1]=False;measured=np.ones(5000,bool)
        p=[np.array([0,2,4]),np.array([1,4,14]),np.array([1,14]),np.array([4])]
        a=env['read_scores'](z,np.array([0,5000]),human,measured,p[:2],lambda v:None)
        b=env['read_all_scores'](z,np.array([0,5000]),human,measured,p,lambda v:None)
        np.testing.assert_array_equal(a[0],b[0][:,:2]);np.testing.assert_array_equal(a[1],b[1]);np.testing.assert_array_equal(a[2],b[2]);self.assertEqual(a[3],b[3])
        self.assertEqual(b[0].shape,(2,4))

    def test_halves_follow_external_positions_exclude_prior_and_do_not_overlap(self):
        genes=[f'synthetic_{j}' for j in range(51)];p=m.split_programmes({'s.genes':['S1','S2'],'g2m.genes':genes})
        self.assertEqual(p['g2m.odd_positions'],genes[::2]);self.assertEqual(p['g2m.even_positions'],genes[1::2])
        self.assertEqual((len(p['g2m.odd_positions']),len(p['g2m.even_positions'])),(26,25))
        with self.assertRaises(AssertionError):m.split_programmes({'g2m.genes':['MKI67']})

    def test_all_LOO_positive_negative_and_singular_replayed_no_case_selection(self):
        rng=np.random.default_rng(64);n=40;c=rng.normal(size=(n,8));cat=np.array([['D' if i%2 else 'D+T','positive' if i%3 else 'negative'] for i in range(n)])
        x=rng.normal(size=n);y=x+.1*rng.normal(size=n)
        for sign in [-1,1]:
            r=m.leave_one_case(x,sign*y,c,cat);self.assertEqual(r['total_case_omissions'],n);self.assertEqual(r['estimable_omissions'],n)
            self.assertEqual(r['positive'] if sign>0 else r['negative'],n)
        c[:,2]=c[:,1];r=m.leave_one_case(x,y,c,cat);self.assertEqual(r['nonestimable_reasons'],{'SINGULAR_FULL_DESIGN':n})

    def test_all52_fine_rows_and_rare_phase_expression_suppression(self):
        tabs=[];frames=[]
        for case in range(5):
            lib=f'synthetic_library_{case}';tabs.append(dict(library_id=lib,patient_registry=f'synthetic_case_{case}',timepoint_registry='pre'))
            frames.append(pd.DataFrame(dict(position=np.arange(case*20,(case+1)*20),library_id=lib,timepoint='pre',major='T',node_id='member',label='Prolif Treg')))
        tab=pd.DataFrame(tabs);frame=pd.concat(frames,ignore_index=True);scores=np.repeat(np.arange(5)/10,20)[:,None]*np.ones((1,4))
        rows=m.fine_profiles(tab,frame,frame,{'member'},scores);self.assertEqual(len(rows),52)
        accepted=next(r for r in rows if r['context']==m.inputs.CONTEXTS[0] and r['phase']=='pre' and r['author_fine_label']=='Prolif Treg')
        self.assertEqual(accepted['cases_atleast20_fine_cells'],5);self.assertEqual(accepted['expression_summaries']['G2M']['case_mean_min_q25_median_q75_max'],[0,.1,.2,.3,.4])
        rows=m.fine_profiles(tab.iloc[:4],frame.iloc[:80],frame.iloc[:80],{'member'},scores[:80])
        self.assertTrue(all(not r['expression_summaries'] for r in rows));self.assertFalse(any(r['phase']=='post2' for r in rows))


if __name__=='__main__':unittest.main()
