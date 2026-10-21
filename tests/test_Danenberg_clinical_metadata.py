from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import qualify_Danenberg_clinical_metadata_v1 as m
class Metadata(unittest.TestCase):
    def row(self):
        r={f:'' for f in m.FIELDS};r.update(metabric_id='synthetic_A',yearsToStatus='3.5',DeathBreast='1',isValidation='TRUE',ERStatus='reported_positive',ERBB2_pos='FALSE',Grade='reported_grade',LymphNodesOrdinal='reported_nodes',sizeOrdinal='reported_size');return r
    def test_duplicate_identity_not_multiple_patients_or_conflicting_first_row_selection(self):
        r=self.row();v=m.summarize([r,r]);self.assertEqual(v['distinct_nonmissing_author_patient_keys'],1);self.assertEqual(v['duplicate_patient_rows'],1)
        a=self.row();a['DeathBreast']='0';v=m.summarize([r,a]);self.assertEqual(v['discordant_patient_groups'],1)
    def test_unknown_event_or_flag_or_nonpositive_time_remains_unavailable(self):
        for field,value,reason in [('DeathBreast','2','UNKNOWN_REPORTED_EVENT_CODE'),('isValidation','unreported','UNKNOWN_AUTHOR_VALIDATION_FLAG'),('yearsToStatus','0','NO_POSITIVE_FINITE_REPORTED_FOLLOWUP')]:
            r=self.row();r[field]=value;v=m.summarize([r]);self.assertEqual(v['patient_qualification_reasons'][reason],1)
    def test_validation_stratum_and_metadata_are_not_claimed_cell_or_clinical_validation(self):
        v=m.summarize([self.row()]);self.assertEqual(v['by_author_validation_flag'][0]['reported_breast_death_events'],1);self.assertTrue(v['cell_state_capture_and_complete_model_eligibility_pending']);self.assertTrue(v['author_flag_not_automatically_independent_cohort'])
if __name__=='__main__':unittest.main()
