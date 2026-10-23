"""Counterexamples and a complete synthetic source-held-out pipeline replay."""
from pathlib import Path
from types import SimpleNamespace
import os,sys,json,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    import numpy as np
    import scipy.sparse as sp
    from annotation_representation_pilot_v1 import normalized_log_counts,select_features
    from run_annotation_source_holdout_v1 import main,inverse_neighbor_prediction,conditional_scores
    AVAILABLE=True
except ImportError:
    AVAILABLE=False


@unittest.skipUnless(AVAILABLE and 'ATLAS_SYMPHONY_WHEEL' in os.environ,'Run in qualified reference runtime')
class HoldoutTests(unittest.TestCase):
    def test_neighbor_prediction_exact_match_and_axis_failure(self):
        z=np.array([[0.],[1.],[2.]]);target=np.array([[1.],[3.],[5.]])
        y,d=inverse_neighbor_prediction(z,target,z,1)
        np.testing.assert_allclose(y,target);np.testing.assert_allclose(d,0)
        with self.assertRaises(ValueError):inverse_neighbor_prediction(z,target[:-1],z,1)

    def test_scoring_missing_flat_and_variable_signal(self):
        x=np.arange(100,dtype=float)
        v=conditional_scores(x,x,np.ones(100),20)
        self.assertAlmostEqual(v['depth_adjusted_Pearson'],1.)
        self.assertEqual(v['MAE_native_logCP10K'],0)
        flat=conditional_scores(np.ones(100),np.ones(100),np.ones(100),20)
        self.assertEqual(flat['depth_adjusted_structure_status'],'no_variable_observed_or_predicted_signal')
        missing=conditional_scores(x,np.full(100,np.nan),np.ones(100),20)
        self.assertEqual(missing['status'],'insufficient_evaluable_query_cells')

    def test_three_source_end_to_end_fits_exclude_query_rows(self):
        root=Path(tempfile.mkdtemp(prefix='annotation-source-holdout-counterexample-'))
        for d in ['config','provenance','outputs','checkpoints']:(root/d).mkdir()
        manifest=root/'dummy_manifest.json';manifest.write_text('{}\n')
        rng=np.random.default_rng(6)
        count=sp.csr_matrix(rng.poisson(2,(180,60)),dtype=np.int32)
        total=np.asarray(count.sum(axis=1)).ravel()
        x,_=normalized_log_counts(count,total)
        source=np.repeat(['s0','s1','s2'],60);fold=np.repeat(['f0','f1','f2'],60)
        routes=np.tile(np.repeat(['T_NK','Myeloid','Epithelial'],20),3)
        state1=x[:,0].toarray().ravel();state2=x[:,1].toarray().ravel()
        data={'x':x,'genes':['G'+str(i) for i in range(60)],'source':source,'fold':fold,'route':routes,
              'fine':np.tile(['a','b','c','d'],45),'native_total':total,
              'state':{'p1':state1,'p2':state2},'state_genes':{'p1':['G0'],'p2':['G1']}}
        cfg={'wheel_path':os.environ['ATLAS_SYMPHONY_WHEEL'],'wheel_SHA256':os.environ['ATLAS_SYMPHONY_SHA256'],
             'bindings':{},'read_block_rows':128,'expected_fold_groups':['f0','f1','f2'],
             'n_features':20,'seed':1,'components':5,'pca_iterations':2,'theta':1.,'lambda':1.,
             'max_iter_harmony':2,'threads':4,'neighbors':6,'resolution':1.,'minimum_training_group':20,
             'reference_routes':['T_NK','Myeloid','Epithelial'],'anchor_cap_per_route':30,
             'state_neighbors':5,'prototype_neighbors':3,'minimum_query_evaluation':5,'output_cap_bytes':16*1024**2}
        (root/'config/source_holdout_spec_v1.json').write_text(json.dumps(cfg))
        calls=[]
        def observed_selection(matrix,labels,n):
            expected=np.flatnonzero(fold!=cfg['expected_fold_groups'][len(calls)])
            np.testing.assert_array_equal(labels,source[expected])
            np.testing.assert_array_equal(matrix.toarray(),x[expected].toarray())
            calls.append(len(expected))
            return select_features(matrix,labels,n)
        with patch('run_annotation_source_holdout_v1.load_pilot',return_value=data),patch('run_annotation_source_holdout_v1.select_features',side_effect=observed_selection):
            main(SimpleNamespace(output=str(root),pilot='unused_synthetic',manifest=str(manifest)))
        self.assertEqual(calls,[120,120,120])
        terminal=json.loads((root/'provenance/terminal_v1.json').read_text())
        self.assertEqual(terminal['n_folds'],3)
        self.assertFalse(terminal['method_winner_or_identity_accuracy_accepted'])
        scopes=json.loads((root/'outputs/source_heldout_qualification_v1.json').read_text())
        for row in scopes:
            self.assertFalse(set(row['training_sources'])&set(row['query_sources']))
            self.assertTrue(all(q['C_reference_unchanged'] for q in row['routes']))
        with np.load(root/'models/fold_0/A_B_shared_minimal_reference.npz',allow_pickle=False) as model:
            self.assertEqual(model['loadings'].shape,(20,5))
            self.assertEqual(len(model['genes']),20)


if __name__=='__main__':unittest.main(verbosity=2)
