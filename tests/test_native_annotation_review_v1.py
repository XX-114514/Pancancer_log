"""Native evidence, abstention, minority gamma-delta and reference API checks."""
from pathlib import Path
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
try:
    import numpy as np
    import scipy.sparse as sp
    from native_annotation_review_v1 import predict_reference, model_fingerprint, native_group_evidence, marker_calls, classify_group, descriptive_markers
    AVAILABLE = True
except ImportError:
    AVAILABLE = False


@unittest.skipUnless(AVAILABLE, 'Use actual scientific environment')
class NativeReviewTests(unittest.TestCase):
    def test_native_total_zero_RNA_and_group_denominators(self):
        raw = sp.csr_matrix([[2, 0], [0, 4], [0, 0], [1, 2]])
        out = native_group_evidence(raw, ['a', 'b'], np.array([10., 20., 0., 10.]), np.array([0, 0, 1, 1]), 2)
        np.testing.assert_array_equal(out['n'], [2, 1])
        np.testing.assert_array_equal(out['cp_sum'], [[2000, 2000], [1000, 2000]])
        np.testing.assert_array_equal(out['detected'], [[1, 1], [1, 1]])
        for totals, groups in [(np.array([1., 20., 0., 10.]), np.array([0,0,1,1])),
                               (np.array([10., np.nan, 0., 10.]), np.array([0,0,1,1])),
                               (np.array([10., 20., 0., 10.]), np.array([0,0,-1,1]))]:
            with self.assertRaises(ValueError):native_group_evidence(raw, ['a','b'], totals, groups, 2)

    def test_unmeasured_marker_and_single_gene_do_not_support_lineage(self):
        genes = ['CD3D', 'CD3E', 'TRDC', 'TRGC1']
        values = np.array([[1, 0, 1, 1], [1, 1, 0, 0]], bool)
        result = marker_calls(values, genes, np.array([1, 1, 1, 0], bool))
        np.testing.assert_array_equal(result['T'], [False, True])
        self.assertFalse(result['gamma_delta_T_related'].any())
        b = marker_calls(np.ones((2,2), bool), ['CD37','CD74'], np.ones(2,bool))
        self.assertFalse(b['B'].any())

    def test_conflicting_lineages_and_minority_gamma_candidate_not_whole_group_identity(self):
        conflict = classify_group({'T': .6, 'Myeloid': .5}, {}, {}, 100)
        self.assertEqual(conflict['local_major'], 'unresolved')
        minority = classify_group({'T': .7, 'gamma_delta_T_related': .06, 'CD8_T_related': .6}, {}, {'Tem/Trm cytotoxic T cells': 90, 'gamma-delta T cells': 10}, 100)
        self.assertEqual(minority['local_celltype'], 'CD8_T')
        self.assertEqual(minority['gamma_delta_RNA_candidate_fraction_for_review'], .06)
        self.assertEqual(classify_group({'Epithelial': .7}, {}, {}, 100)['malignant_status'], 'not_inferred_from_lineage')

    def test_marker_fallback_is_own_group_expression_not_DEG(self):
        out = descriptive_markers(['RPL3', 'MARKER', 'OTHER'], np.array([20., 5., 1.]), np.array([1., .7, .4]), np.array([1., .5, 2.]), np.array([.5, .2, .8]))
        self.assertEqual([v['gene'] for v in out], ['MARKER'])
        self.assertIn('not_inferential_DEG', out[0]['marker_evidence_type'])

    @unittest.skipUnless('ATLAS_CELLTYPIST_MODEL_DIR' in os.environ, 'Need frozen cached reference models')
    def test_cached_model_prediction_matches_official_API_and_never_mutates_model(self):
        import anndata as ad
        import pandas as pd
        import celltypist
        from celltypist.models import Model
        root = Path(os.environ['ATLAS_CELLTYPIST_MODEL_DIR'])
        rng = np.random.default_rng(412)
        for filename in ['Immune_All_Low.pkl', 'Immune_All_High.pkl', 'Human_Lung_Atlas.pkl', 'Cells_Adult_Breast.pkl']:
            model = Model.load(str(root / filename));fingerprint = model_fingerprint(model)
            genes = list(model.features)[::1][:-50] + ['NATIVE_OTHER_NOT_MODELLED']
            counts = rng.poisson(2, (9, len(genes)))
            values = np.log1p(counts * (10000 / counts.sum(axis=1))[:, None])
            a = ad.AnnData(sp.csr_matrix(values), var=pd.DataFrame(index=genes))
            official = celltypist.annotate(a, model=model, majority_voting=False, mode='prob match', p_thres=.5)
            result = predict_reference(a.X, genes, model)
            np.testing.assert_allclose(result['decision_scores'], official.decision_matrix.to_numpy(), atol=2e-7, rtol=2e-7)
            np.testing.assert_allclose(result['probability_scores'], official.probability_matrix.to_numpy(), atol=2e-8, rtol=2e-8)
            self.assertEqual(fingerprint, model_fingerprint(model))
            missing = predict_reference(a.X[:, :3], genes[:3], model)
            self.assertEqual(missing['status'], 'insufficient_measured_model_features')
            self.assertTrue((missing['unique_code'] == -1).all())
            with self.assertRaises(ValueError):predict_reference(-a.X, genes, model)


if __name__ == '__main__':unittest.main()
