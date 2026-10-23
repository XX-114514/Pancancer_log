import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from cell_identity_confidence_v2 import decide_cell, low_reference_identity


class ConfidenceTests(unittest.TestCase):
    def base(self,**kw):
        a=dict(native_families={'T_NK'},reference_families={'T_NK'},native_major='T',reference_major='T',group_major='T',group_type='CD8_T',own_fine=set(),reference_fine=set())
        a.update(kw);return decide_cell(**a)
    def test_coherent_minority_not_overwritten_by_group(self):
        r=self.base(own_fine={'Treg','CD4_T'},reference_fine={'Treg','CD4_T'})
        self.assertEqual(r['consensus_celltype'],'Treg');self.assertTrue(r['group_type_disagreement'])
    def test_cross_lineage_marker_conflict_not_voted_away(self):
        r=self.base(native_families={'T_NK','Epithelial'})
        self.assertEqual(r['consensus_major'],'unresolved')
    def test_within_family_T_NK_disagreement_retains_broad_identity(self):
        r=self.base(reference_major='NK')
        self.assertEqual(r['consensus_major'],'T_NK');self.assertEqual(r['consensus_celltype'],'T_NK')
    def test_native_dropout_requires_group_corroboration(self):
        r=self.base(native_families=set(),native_major=None,group_major='unresolved')
        self.assertEqual(r['consensus_major'],'unresolved')
        r=self.base(native_families=set(),native_major=None,reference_fine={'CD8_T'})
        self.assertEqual(r['consensus_celltype'],'CD8_T');self.assertIn('inferred',r['coarse_evidence'])
    def test_small_branch_retains_direct_rare_identity_without_fine_fit_claim(self):
        r=self.base(group_major='unresolved',group_type='unresolved',own_fine={'Treg'},reference_fine={'Treg'},fine_fit=False)
        self.assertEqual(r['consensus_celltype'],'Treg');self.assertFalse(r['fine_discovery_fit_available'])
    def test_gamma_reference_only_does_not_prove_TCR_identity(self):
        r=self.base(group_type='T_gamma_delta_RNA_candidate',reference_fine={'T_gamma_delta_RNA_candidate'})
        self.assertEqual(r['consensus_celltype'],'T')
        r=self.base(group_type='T_gamma_delta_RNA_candidate',reference_fine={'T_gamma_delta_RNA_candidate'},own_fine={'T_gamma_delta_RNA_candidate'})
        self.assertIn('requires_orthogonal_review',r['fine_evidence'])
    def test_incompatible_group_preserves_coherent_cell_major(self):
        r=self.base(group_major='Epithelial',group_type='Epithelial')
        self.assertEqual(r['consensus_major'],'T');self.assertTrue(r['group_major_disagreement'])
    def test_zero_RNA_and_unknown_fine_vocabulary_rejected(self):
        self.assertEqual(self.base(positive=False)['consensus_major'],'unresolved')
        with self.assertRaises(ValueError):self.base(own_fine={'invented_state_as_identity'})
        with self.assertRaises(ValueError):self.base(native_families={'invented_family'})
    def test_reference_dropout_against_group_retains_broad_family(self):
        r=self.base(native_families=set(),native_major=None,reference_major='NK',group_major='T')
        self.assertEqual(r['consensus_major'],'T_NK')
        self.assertEqual(r['fine_evidence'],'coarse_identity_retained')
        self.assertEqual(r['coarse_evidence'],'within_family_reference_group_conflict_native_dropout')
    def test_reference_scope_and_precursor_labels(self):
        self.assertEqual(low_reference_identity('Regulatory T cells'),('T',{'Treg','CD4_T'}))
        self.assertEqual(low_reference_identity('Cycling T cells'),('T',set()))
        self.assertEqual(low_reference_identity('pDC precursor'),('Myeloid',set()))
        self.assertEqual(low_reference_identity('Epithelial cells'),(None,set()))


if __name__=='__main__':unittest.main()
