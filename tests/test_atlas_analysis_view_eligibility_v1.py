from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from atlas_analysis_view_eligibility_v1 import qualify_unit_views


def reviewed():
    names=['archive_membership','source_asset_join','native_measurement','prior_QC_scope','expression_unit_join','assay_scope','annotation_release','source_alias_scope','material_and_anatomy','capture_design','unselected_tissue_sampling','lineage_selection_scope','within_lineage_capture_comparability','product_design','patient_identity','patient_sample_library_join','treatment_phase','clinical_endpoint','endpoint_time_origin','repeated_sampling_plan','clinical_cohort_scope']
    return {**dict.fromkeys(names,'accepted'),'capture_class':'unselected_tissue'}


class ViewChecks(unittest.TestCase):
    def test_full_archive_not_a_claim_of_composition_or_patient_eligibility(self):
        r=qualify_unit_views({'archive_membership':'accepted'})
        self.assertTrue(r['all_archive']['eligible'])
        for view in ['expression_reference','whole_tissue_composition','within_selected_lineage_composition','treatment_clinical']:
            self.assertEqual(r[view]['status'],'unknown');self.assertFalse(r[view]['eligible'])

    def test_enrichment_and_mixed_sorted_fractions_never_native_tissue(self):
        for kind in ['immune_enriched','T_enriched_nuclei','selected_lineage','selected_nonmalignant_subset','mixed_sorted_fractions']:
            r=qualify_unit_views({**reviewed(),'capture_class':kind})
            self.assertEqual(r['whole_tissue_composition']['status'],'rejected')
            self.assertTrue(r['within_selected_lineage_composition']['eligible'])
            self.assertTrue(r['expression_reference']['eligible'])
        r=qualify_unit_views({**reviewed(),'capture_class':'unknown'})
        self.assertEqual(r['whole_tissue_composition']['status'],'unknown')
        r=qualify_unit_views({**reviewed(),'capture_class':'whole_tissue_sampling_candidate'})
        self.assertEqual(r['whole_tissue_composition']['status'],'unknown')
        self.assertTrue(r['within_selected_lineage_composition']['eligible'])

    def test_product_and_clinical_are_separate(self):
        r=qualify_unit_views({**reviewed(),'capture_class':'cell_product','patient_identity':'unknown'})
        self.assertTrue(r['cell_product']['eligible']);self.assertFalse(r['whole_tissue_composition']['eligible'])
        self.assertFalse(r['within_selected_lineage_composition']['eligible']);self.assertFalse(r['treatment_clinical']['eligible'])
        self.assertIn('patient_identity:unknown',r['treatment_clinical']['reasons'])
        r=qualify_unit_views({**reviewed(),'clinical_endpoint':'rejected'})
        self.assertTrue(r['whole_tissue_composition']['eligible']);self.assertEqual(r['treatment_clinical']['status'],'rejected')

    def test_quality_and_annotation_gate_do_not_delete_archive(self):
        r=qualify_unit_views({**reviewed(),'annotation_release':'unknown'})
        self.assertTrue(r['all_archive']['eligible']);self.assertTrue(r['expression_reference']['eligible']);self.assertFalse(r['whole_tissue_composition']['eligible'])
        r=qualify_unit_views({**reviewed(),'native_measurement':'rejected'})
        self.assertTrue(r['all_archive']['eligible']);self.assertFalse(r['expression_reference']['eligible'])
        for bad in [{'archive_membership':'accepted','patient_identity':True},{'archive_membership':'unknown'},{'archive_membership':'accepted','capture_class':'tumor'}]:
            with self.assertRaises(ValueError):qualify_unit_views(bad)


if __name__=='__main__':unittest.main()
