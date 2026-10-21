from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import qualify_Hoch_block_landmark_scope_v1 as m


def row(patient='synthetic_A',block='block_A',image='image_A',surgery='01-Jan-2020',death='01-Feb-2020'):
    value={f:'' for f in m.FIELDS}
    value.update(PatientID=patient,BlockID=block,ImageNumber=image,Date_surgery=surgery,
        Start_treatment='02-Jan-2020',Date_death=death,Last_PET='01-Feb-2020',
        Time_to_death_or_last_PET='30',Time_to_progression_or_last_PET='30',
        Nr_treatments_before_surgery='0')
    return value


class BlockLandmark(unittest.TestCase):
    def test_repeated_images_collapse_but_different_blocks_are_not_patient_conflict(self):
        a=row();b=row(image='image_B');c=row(block='block_B',image='image_C',surgery='15-Jan-2020')
        c['Drug_after_surgery']='different reported episode'
        v=m.summarize([a,b,c]);self.assertEqual(v['accepted_block_contexts'],2)
        self.assertEqual(v['tissue_landmark_eligible_metadata_upper_bound'],1)
        self.assertEqual(v['source_treatment_origin_duration_replay']['TREATMENT_ORIGIN_DURATION_MATCH'],4)

    def test_earliest_selection_precedes_followup_not_replaced_with_later_eligible_tissue(self):
        a=row(death='01-Jan-2020');b=row(block='block_B',image='image_B',surgery='15-Jan-2020')
        v=m.summarize([a,b]);self.assertEqual(v['tissue_landmark_eligible_metadata_upper_bound'],0)
        self.assertEqual(v['earliest_tissue_selection_reasons']['ENDPOINT_NOT_STRICTLY_AFTER_SELECTED_TISSUE'],1)

    def test_invalid_death_not_converted_to_censoring_and_missing_prior_unknown(self):
        a=row(death='unknown calendar');v=m.summarize([a])
        self.assertEqual(v['earliest_tissue_selection_reasons']['UNPARSED_DEATH_NOT_CENSORING'],1)
        a=row();a['Nr_treatments_before_surgery']='';v=m.summarize([a])
        self.assertEqual(v['selected_phase_context_groups'][0]['prior_treatment'],'PRIOR_COUNT_UNKNOWN')

    def test_tied_earliest_blocks_are_not_arbitrarily_selected(self):
        v=m.summarize([row(),row(block='block_B',image='image_B')])
        self.assertEqual(v['earliest_tissue_selection_reasons']['TIED_EARLIEST_BLOCKS_NOT_ARBITRARILY_CHOSEN'],1)

    def test_within_block_conflict_and_cross_patient_image_collisions_preserved(self):
        a=row();b=row(image='image_B');b['Date_surgery']='02-Jan-2020';v=m.summarize([a,b])
        self.assertEqual(v['unresolved_block_reasons']['WITHIN_BLOCK_CONTEXT_CONFLICT'],1)
        v=m.summarize([row(),row(patient='synthetic_B',block='block_B')])
        self.assertEqual(v['image_identity_collision_keys'],1)
        self.assertEqual(v['tissue_landmark_eligible_metadata_upper_bound'],0)


if __name__=='__main__':unittest.main()
