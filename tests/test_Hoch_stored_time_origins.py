from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import audit_Hoch_stored_time_origins_v1 as m
class TimeOrigins(unittest.TestCase):
    def row(self):
        r={f:'' for f in m.scope.FIELDS};r.update(PatientID='synthetic',BlockID='synthetic_block',ImageNumber='synthetic_image',Date_surgery='01-Jan-2020',Start_treatment='02-Jan-2020',Date_death='01-Feb-2020',Last_PET='01-Feb-2020',Time_to_death_or_last_PET='31',Time_to_progression_or_last_PET='31');return r
    def test_explicit_two_anchors_not_a_winning_replacement(self):
        a=m.audit([self.row()]);death=a['rows'][:2];self.assertEqual(death[0]['audit_counts']['EXACT_DAYS_AGREEMENT'],1);self.assertEqual(death[1]['audit_counts']['DAYS_MISMATCH'],1)
        self.assertTrue(a['no_new_unit_guesses_or_stored_value_replacement'])
    def test_missing_death_uses_last_PET_but_unparsed_death_does_not(self):
        r=self.row();r['Date_death']='';a=m.audit([r]);self.assertEqual(a['rows'][0]['audit_counts']['EXACT_DAYS_AGREEMENT'],1)
        r['Date_death']='ambiguous';a=m.audit([r]);self.assertEqual(a['rows'][0]['audit_counts']['UNAVAILABLE_DATE_OR_STORED_VALUE'],1)
if __name__=='__main__':unittest.main()
