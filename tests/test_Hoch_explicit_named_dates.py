import importlib.util
from datetime import date
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('named_dates',ROOT/'scripts/qualify_Hoch_melanoma_clinical_metadata_v2.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class NamedDates(unittest.TestCase):
    def test_named_month_calendar_and_iso(self):
        self.assertEqual(m.explicit_date('29-Feb-2024'),date(2024,2,29))
        self.assertEqual(m.explicit_date('03-jAn-2020'),date(2020,1,3))
        self.assertEqual(m.explicit_date('2020-01-03'),date(2020,1,3))
        for value in ['29-Feb-2023','03/04/2020','03-Apr-20','03-Abc-2020',None]:
            self.assertIsNone(m.explicit_date(value))

    def test_summary_preserves_patient_conflicts_and_no_record_export(self):
        def row(patient,surgery,start):return dict(PatientID=patient,ImageNumber=patient+str(surgery),
            Date_surgery=surgery,Start_treatment=start,Drug_after_surgery='pembrolizumab',
            Treatment_after_surgery='anti-PD1',Adjuvant='0',Nr_treatments_before_surgery='0',
            Cancer_Stage='IV',Response_at_3m='synthetic_response')
        rows=[row('synthetic_A','01-Jan-2020','01-Feb-2020'),row('synthetic_A','01-Jan-2020','01-Feb-2020'),
              row('synthetic_B','01-Mar-2020','01-Feb-2020'),
              row('synthetic_C','01-Jan-2020','01-Feb-2020'),row('synthetic_C','02-Jan-2020','01-Feb-2020')]
        result=m.summarize(rows)
        self.assertEqual(result['author_patient_keys'],3)
        self.assertEqual(result['patients_with_context_conflict'],1)
        self.assertEqual(result['reason_counts']['strict_baseline_nonadjuvant_ICI_metadata_upper_bound'],1)
        self.assertEqual(result['source_date_format_counts']['Date_surgery']['explicit_ISO_or_named_calendar_month_parsed'],2)
        self.assertNotIn('synthetic_',str(result))


if __name__=='__main__':unittest.main()
