import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from library_source_qualification_v1 import qualify_library,normalize_library_alias,resolve_library_aliases

class LibraryPolicy(unittest.TestCase):
    def base(self):return {'experimental_sorting':'unselected','material':'primary_tumor','RNA_modality':'scRNA','retained_data_scope':'whole_tissue','actual_sample_or_library_bridge':True}
    def test_T_proportion_never_certifies_enrichment(self):
        r=self.base();r['experimental_sorting']='unknown';r['T_cell_fraction']=.99
        self.assertEqual(qualify_library(r)['whole_TME_scRNA_main'],'待核实')
        self.assertFalse(qualify_library(r)['sorting_inferred_from_cell_proportions'])
    def test_CD45_and_CD8_define_different_denominators(self):
        a=self.base();a['experimental_sorting']='CD45_positive';b=self.base();b['experimental_sorting']='CD8_positive'
        self.assertEqual(qualify_library(a)['whole_TME_scRNA_main'],'排除')
        self.assertEqual(qualify_library(a)['T_internal_CM'],'条件候选')
        self.assertIn('限定CD8',qualify_library(b)['T_internal_scope'])
    def test_same_study_can_have_different_library_eligibility(self):
        a=self.base();a['study']='same';b={**a,'experimental_sorting':'other_lineage_selected'}
        self.assertEqual(qualify_library(a)['whole_TME_scRNA_main'],'条件候选')
        self.assertEqual(qualify_library(b)['whole_TME_scRNA_main'],'排除')
    def test_blood_product_nuclei_and_exvivo_remain_separate(self):
        for change in [{'material':'blood'},{'cell_product':True},{'RNA_modality':'snRNA'},{'ex_vivo_stimulation':'yes'}]:
            a={**self.base(),**change};self.assertEqual(qualify_library(a)['whole_TME_scRNA_main'],'排除');self.assertTrue(qualify_library(a)['archive_retained'])
    def test_selected_file_scope_is_not_called_experimental_sorting(self):
        a=self.base();a['experimental_sorting']='unknown';a['retained_data_scope']='selected_nonmalignant_subset'
        q=qualify_library(a);self.assertEqual(q['whole_TME_scRNA_main'],'排除');self.assertIn('具体分选/富集方法',q['critical_missing_fields'])
    def test_mixed_study_missing_library_bridge_stays_unknown(self):
        a=self.base();a['mixed_study_library_branch_unresolved']=True;a['actual_sample_or_library_bridge']=False
        q=qualify_library(a);self.assertEqual(q['whole_TME_scRNA_main'],'待核实');self.assertFalse(q['patient_independence_newly_certified'])
    def test_sorted_fraction_weights_change_the_question_explicitly(self):
        a=self.base();a['experimental_sorting']='mixed_sorted_fractions'
        self.assertEqual(qualify_library(a)['whole_TME_scRNA_main'],'排除')
        a['mixing_weights_reconstructable']=True
        self.assertEqual(qualify_library(a)['whole_TME_scRNA_main'],'条件候选')
    def test_CD45_negative_is_not_whole_TME_or_T_capture(self):
        r={**self.base(),'experimental_sorting':'CD45_negative'}
        self.assertEqual(qualify_library(r)['whole_TME_scRNA_main'],'排除')
        self.assertEqual(qualify_library(r)['T_internal_CM'],'排除')
    def test_all_cell_file_name_does_not_certify_unsorted_tissue(self):
        r={**self.base(),'experimental_sorting':'unknown','retained_data_scope':'whole_tissue','file_title':'All cells'}
        self.assertEqual(qualify_library(r)['whole_TME_scRNA_main'],'待核实')
    def test_identifier_token_boundaries_prevent_false_library_join(self):
        self.assertNotEqual(normalize_library_alias('BC2 1'),normalize_library_alias('BC21'))
        q=resolve_library_aliases(['BC21'],{'a':['BC2 1'],'b':['BC21']})
        self.assertEqual(q,{'reference_keys':['b'],'status':'unique'})
    def test_ambiguous_aliases_and_replicates_stay_unresolved(self):
        self.assertEqual(resolve_library_aliases(['A-1'],{'a':['A_1'],'b':['A 1']})['status'],'ambiguous')
        self.assertEqual(resolve_library_aliases(['P13.post.2'],{'x':['P13.post.3']})['status'],'unmatched')

if __name__=='__main__':unittest.main()
