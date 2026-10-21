from pathlib import Path
import ast,sys,unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import CM02_G2M_additive_diagnostics_v1 as m
import ucell_native_sparse_scores_v1 as original
import RNA_programme_case_inputs_v1 as inputs


class Diagnostics(unittest.TestCase):
    def test_group_profiles_suppress_rare_and_condition_concentration_correctly(self):
        tabs=[];parts=[]
        for i in range(5):
            library=f'synthetic_library_{i}';tabs.append(dict(library_id=library,patient_registry=f'synthetic_case_{i}'))
            parts.append(pd.DataFrame(dict(position=np.arange(i*20,(i+1)*20),library_id=library,timepoint='pre',major='T',node_id='member',label='Prolif Treg')))
        frame=pd.concat(parts,ignore_index=True);d=np.zeros((100,2,4));a=np.tile([0,1,2,5,10],20);d[:,1,0]=d[:,1,1]=d[:,1,2]=a;d[:,1,3]=np.divide(1.,a,out=np.zeros(100),where=a>0)
        rows=m.support_profiles(pd.DataFrame(tabs),frame,{'member'},d,inputs);self.assertEqual(len(rows),64)
        r=next(r for r in rows if r['context']==inputs.CONTEXTS[0] and r['phase']=='pre' and r['reference_kind']=='fine' and r['reference_label']=='Prolif Treg')
        self.assertAlmostEqual(r['metrics']['active']['case_mean_min_q25_median_q75_max'][2],3.6)
        self.assertAlmostEqual(r['metrics']['maximum_fraction']['case_mean_min_q25_median_q75_max'][2],.45)
        self.assertAlmostEqual(r['metrics']['fraction_active_atleast5']['case_mean_min_q25_median_q75_max'][2],.4)
        self.assertTrue(all(not r['metrics'] for r in rows if r['cases_atleast20_observed_cells']<5))

    def test_extended_CSR_reader_keeps_closed_scores_and_replays_diagnostics(self):
        from types import SimpleNamespace
        import resource
        root=Path(__file__).resolve().parents[1];env=dict(np=np,pd=pd,baseline=SimpleNamespace(scoring=original),diagnostic=m,
            resource=SimpleNamespace(RUSAGE_SELF=0,getrusage=lambda x:SimpleNamespace(ru_maxrss=1024)))
        function=next(n for n in ast.parse((root/'scripts/measure_CM02_G2M_decomposition_v1.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='read_extended')
        exec(compile(ast.Module(body=[function],type_ignores=[]),'extended_reader','exec'),env)
        class Raw(dict):attrs={'shape':(2,5000)}
        raw=Raw(indptr=np.array([0,4,7]),indices=np.array([4,4,2,4999,0,2,14]),data=np.array([2,3,2,500,4,6,1]))
        human=np.ones(5000,bool);human[-1]=False;measured=np.ones(5000,bool);programmes=[np.array([0,2,4]),np.array([1,4,14])]
        expected=original.score_selected_CSR_rows(raw['indptr'],raw['indices'],raw['data'],np.array([0,1]),human,measured,programmes)
        got=env['read_extended']({'layers/counts':raw},np.array([0,1]),human,measured,programmes,lambda x:None)
        for i in range(3):np.testing.assert_array_equal(expected[i],got[i])
        self.assertEqual(got[3].shape,(2,2,4));self.assertEqual(got[4]['independent_dense_boundary_rows'],2)

    def test_gene_support_with_ties_zeros_and_dense_rank_math(self):
        rng=np.random.default_rng(37)
        for npositive in [0,30,1600,2200]:
            values=np.zeros(5000);values[:npositive]=rng.integers(1,12,npositive)
            ids=[np.array([0,2,16,44,102,1800,4000]),np.array([1,5,25,1999,2500,4999])];query=[values[p] for p in ids]
            score,diag=m.signature_diagnostics(values[values>0],query,5000)
            np.testing.assert_allclose(score,original.score_positive_sparse_cell(values[values>0],query,5000),atol=1e-12)
            rank=pd.Series(-values).rank(method='average').to_numpy()
            for i,p in enumerate(ids):
                deficit=1500-np.minimum(rank[p],1500);mass=deficit.sum()
                expected=[np.count_nonzero(values[p]),np.count_nonzero(deficit),mass**2/(deficit@deficit) if mass else 0,deficit.max()/mass if mass else 0]
                np.testing.assert_allclose(diag[i],expected,atol=1e-12)

    def fixture(self):
        rng=np.random.default_rng(22);n=40;lc=rng.integers(20,80,(n,2));rc=rng.integers(20,80,(n,3))
        ls=lc*rng.uniform(.01,.3,lc.shape);rs=rc*rng.uniform(.01,.3,rc.shape)
        z=np.column_stack([np.ones(n),rng.normal(size=(n,4))]);return lc,rc,ls,rs,z

    def test_all_additive_covariance_terms_match_independent_total_residual(self):
        lc,rc,ls,rs,z=self.fixture();result,fine=m.additive_decomposition(lc,rc,ls,rs,z,['L1','L2'],['R1','R2','R3'])
        x=ls.sum(axis=1)/lc.sum(axis=1);y=rs.sum(axis=1)/rc.sum(axis=1)
        rx=x-z@np.linalg.lstsq(z,x,rcond=None)[0];ry=y-z@np.linalg.lstsq(z,y,rcond=None)[0]
        self.assertAlmostEqual(result['residual_raw_score_Pearson'],np.corrcoef(rx,ry)[0,1],places=11)
        self.assertEqual(len(fine),6)
        for r in result['reference_decompositions']:self.assertAlmostEqual(sum(t['normalized_covariance_contribution'] for t in r['terms']),result['residual_raw_score_Pearson'],places=11)

    def test_composition_only_and_within_only_are_not_confused(self):
        lc,rc,ls,rs,z=self.fixture()
        for mode in ['composition','within']:
            if mode=='composition':a,b=lc,rc;left=lc*np.array([.05,.25]);right=rc*np.array([.04,.1,.28])
            else:a=np.full_like(lc,40);b=np.full_like(rc,40);left=a*(ls/lc);right=b*(rs/rc)
            r,f=m.additive_decomposition(a,b,left,right,z,['L1','L2'],['R1','R2','R3'])
            index=0 if mode=='composition' else 3
            for ref in r['reference_decompositions']:
                self.assertAlmostEqual(ref['terms'][index]['normalized_covariance_contribution'],r['residual_raw_score_Pearson'],places=10)
                self.assertLess(sum(abs(t['normalized_covariance_contribution']) for j,t in enumerate(ref['terms']) if j!=index),1e-10)

    def test_captured_zero_and_singular_controls_are_explicit(self):
        lc,rc,ls,rs,z=self.fixture();lc[:,0]=0;ls[:,0]=0
        r,f=m.additive_decomposition(lc,rc,ls,rs,z,['uncaptured','L2'],['R1','R2','R3'])
        self.assertTrue(all(t['residual_covariance']==0 for t in f if t['left_author_fine']=='uncaptured'))
        z[:,2]=z[:,1]
        with self.assertRaises(AssertionError):m.additive_decomposition(lc,rc,ls,rs,z,['L1','L2'],['R1','R2','R3'])


if __name__=='__main__':unittest.main()
