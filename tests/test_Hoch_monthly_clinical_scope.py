from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import qualify_Hoch_monthly_clinical_scope_v2 as m
class MonthlyScope(unittest.TestCase):
    def row(self):
        r={f:'' for f in m.FIELDS};r.update(PatientID='synthetic_A',BlockID='synthetic_block',ImageNumber='synthetic_image',Date_surgery='01-Jan-2020',Start_treatment='01-Feb-2020',Date_death='01-Mar-2020',Last_PET='01-Mar-2020',Time_to_death_or_last_PET='35',Time_to_progression_or_last_PET='35',censoring_death='1',censoring_progression='0');return r
    def test_month_compatible_reported_duration_is_not_exact_daily_replay(self):
        r=self.row();self.assertEqual(m.duration_check(r,'Date_death','Time_to_death_or_last_PET'),'COMPATIBLE_WITH_EXPORTED_MONTH_INTERVALS')
        v=m.summarize([r]);self.assertEqual(v['preindex_reported_OS_metadata_only_opportunities']['MONTH_PREINDEX_REPORTED_OS_METADATA_OPPORTUNITY'],1)
        self.assertTrue(v['monthly_compatibility_not_exact_daily_verification'])
    def test_same_month_is_unresolved_and_hidden_days_not_guessed(self):
        r=self.row();r['Date_surgery']='01-Feb-2020';v=m.summarize([r]);self.assertEqual(v['selected_tissue_index_month_relationship']['SAME_MONTH_ORDER_UNRESOLVED'],1)
        self.assertNotIn('MONTH_PREINDEX_REPORTED_OS_METADATA_OPPORTUNITY',v['preindex_reported_OS_metadata_only_opportunities'])
    def test_outside_interval_and_nonfirst_day_rejected_without_repair(self):
        r=self.row();r['Time_to_death_or_last_PET']='1000';self.assertEqual(m.duration_check(r,'Date_death','Time_to_death_or_last_PET'),'OUTSIDE_EXPORTED_MONTH_INTERVALS')
        r['Date_surgery']='15-Jan-2020'
        with self.assertRaises(ValueError):m.summarize([r])
if __name__=='__main__':unittest.main()
