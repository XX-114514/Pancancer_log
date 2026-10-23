import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_source_analysis_roles_v1 import unit_identity_coverage,source_role


class SourceRoleTests(unittest.TestCase):
    def test_identifier_not_automatically_reviewed(self):
        r=unit_identity_coverage(['u','u'],['p','p'],['False','False'])
        self.assertEqual(r['n_units_complete_single_patient_key_unreviewed_or_reviewed'],1)
        self.assertEqual(r['n_units_complete_single_patient_key_and_all_rows_verified'],0)
    def test_multiple_patient_unit_not_promoted(self):
        r=unit_identity_coverage(['u','u','v'],['p','q','p'],['True']*3)
        self.assertEqual(r['n_units_multiple_patient_keys'],1)
        self.assertEqual(r['n_units_complete_single_patient_key_and_all_rows_verified'],1)
    def test_partial_missing_identity_not_complete(self):
        r=unit_identity_coverage(['u','u','v'],['p','unknown','nan'],['True']*3)
        self.assertEqual(r['n_units_partial_patient_key_coverage'],1)
        self.assertEqual(r['n_units_no_usable_patient_key'],1)
        self.assertEqual(r['n_units_complete_single_patient_key_and_all_rows_verified'],0)
    def test_capture_roles_do_not_claim_whole_tissue(self):
        for c in ['immune_enriched','T_enriched_nuclei','cell_product','selected_nonmalignant_subset']:
            self.assertTrue(source_role({'capture_class':c})['whole_tissue_composition'].startswith('not_eligible'))
        self.assertEqual(source_role({'capture_class':'unknown'})['whole_tissue_composition'],'unqualified_capture_unknown')


if __name__=='__main__':unittest.main()
