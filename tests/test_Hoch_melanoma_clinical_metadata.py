"""Patient counts do not become core counts; ambiguous dates stay unresolved."""
import importlib.util
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('Hoch_meta',ROOT/'scripts/qualify_Hoch_melanoma_clinical_metadata_v1.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def row(identity,image,response='PR'):
    return dict(PatientID=identity,ImageNumber=image,Date_surgery='2019-01-01',
        Start_treatment='2019-02-01',Drug_after_surgery='Pembrolizumab',
        Treatment_after_surgery='immunotherapy',Adjuvant='0',Nr_treatments_before_surgery='0',
        Cancer_Stage='IV',Response_at_3m=response)


class ClinicalMetadata(unittest.TestCase):
    def test_only_explicit_ISO_dates_and_known_drug_text(self):
        self.assertIsNotNone(m.iso_date('2019-02-01'))
        for value in ['02/03/2019','2019-02-31','UNKNOWN','']:
            self.assertIsNone(m.iso_date(value))
        self.assertEqual(m.regimen('Nivolumab plus Ipilimumab'),'PD1_and_CTLA4_text_present_not_schedule_equivalence')
        self.assertEqual(m.regimen('unclassified'),'OTHER_OR_UNCLASSIFIED_TEXT')
    def test_multiple_cores_are_one_author_patient(self):
        result=m.summarize([row('synthetic_A','1'),row('synthetic_A','2'),row('synthetic_B','3','')])
        self.assertEqual(result['author_patient_keys'],2)
        self.assertEqual(result['distinct_nonmissing_image_keys'],3)
        self.assertEqual(result['reason_counts']['strict_baseline_nonadjuvant_ICI_metadata_upper_bound'],2)
        self.assertEqual(result['reason_counts']['strict_baseline_nonadjuvant_ICI_with_concordant_nonmissing_3m_response_upper_bound'],1)
        self.assertNotIn('synthetic_A',str(result))
    def test_conflicting_patient_dates_and_ambiguous_image_links_are_flagged(self):
        a=row('synthetic_A','1');b=row('synthetic_A','2');b['Date_surgery']='2020-01-01'
        result=m.summarize([a,b,row('synthetic_B','1')])
        self.assertEqual(result['patients_with_context_conflict'],1)
        self.assertEqual(result['image_keys_linked_to_multiple_patients'],1)
    def test_same_day_or_prior_treatment_does_not_pass_naive_baseline_upper_bound(self):
        a=row('synthetic_A','1');a['Start_treatment']='2019-01-01'
        b=row('synthetic_B','2');b['Nr_treatments_before_surgery']='1'
        result=m.summarize([a,b]);self.assertEqual(result['reason_counts'].get('strict_baseline_nonadjuvant_ICI_metadata_upper_bound',0),0)


if __name__=='__main__':unittest.main()
