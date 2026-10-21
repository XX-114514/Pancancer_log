"""Block-aware clinical opportunity audit, no clinical association or cell read.

Source rows and all identifiers remain RAM only. Multiple blocks within a
patient are not automatically conflicting episodes. The earliest observed
dated image-linked block is selected before follow-up eligibility is inspected.
Post-biopsy follow-up is a different target from pretreatment ICI prediction.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import argparse,csv,hashlib,importlib.util,io,json,os,math,zipfile
import qualify_Hoch_melanoma_clinical_metadata_v2 as dates

DIRECTORY_SHA='dff220e35d38da2b60c42950190d71b2259f12aae3c21ad4bc6c05ad5c5d2893'
CLINICAL_SHA='ac4a02e4d0e23504aa980b21585de74c4d7390d2e6e0ec2399303c5c29ca5b25'
FIELDS=['PatientID','BlockID','ImageNumber','Date_surgery','Start_treatment',
        'Last_PET','Date_death','Date_progression','Nr_treatments_before_surgery',
        'Drug_after_surgery','Treatment_after_surgery','Adjuvant','Response_at_3m',
        'Status_at_3m','Cancer_Stage','Age_range','Gender',
        'Time_to_death_or_last_PET','Time_to_progression_or_last_PET']
DATE_FIELDS=['Date_surgery','Start_treatment','Last_PET','Date_death','Date_progression']


def canonical(value):return '' if dates.core.missing(value) else str(value).strip()


def endpoint(context,anchor):
    """Unknown/unparsed death is not censoring; never choose a later block."""
    if canonical(context['Date_death']):
        end=dates.explicit_date(context['Date_death']);event=True
        if end is None:return None,'UNPARSED_DEATH_NOT_CENSORING'
    else:
        end=dates.explicit_date(context['Last_PET']);event=False
        if end is None:return None,'NO_ACCEPTED_LAST_PET_CENSOR_DATE'
    if end<=anchor:return None,'ENDPOINT_NOT_STRICTLY_AFTER_SELECTED_TISSUE'
    return (end,event),'POSITIVE_POST_TISSUE_ENDPOINT_OPPORTUNITY'


def summarize(rows):
    assert rows and set(FIELDS)<=set(rows[0])
    selected=[{k:canonical(row.get(k)) for k in FIELDS} for row in rows]
    images=defaultdict(set);block_patients=defaultdict(set);grouped=defaultdict(list);reasons=Counter();date_formats={f:Counter() for f in DATE_FIELDS}
    for row in selected:
        for f in DATE_FIELDS:date_formats[f]['MISSING' if not row[f] else 'ACCEPTED_CALENDAR' if dates.explicit_date(row[f]) else 'UNACCEPTED_CALENDAR']+=1
        if not row['PatientID'] or not row['BlockID'] or not row['ImageNumber']:
            reasons['ROW_WITHOUT_PATIENT_BLOCK_IMAGE_KEY']+=1;continue
        key=(row['PatientID'],row['BlockID']);images[row['ImageNumber']].add(key);block_patients[row['BlockID']].add(row['PatientID']);grouped[key].append(row)
    bad_images={k for k,v in images.items() if len(v)>1};bad_blocks={k for k,v in block_patients.items() if len(v)>1}
    contexts={};conflicts=Counter();unresolved=Counter();patients=defaultdict(list);duration=Counter()
    context_fields=[f for f in FIELDS if f not in ['PatientID','BlockID','ImageNumber']]
    for key,records in grouped.items():
        patient,block=key;patients[patient].append(key)
        if block in bad_blocks or any(r['ImageNumber'] in bad_images for r in records):
            unresolved['IDENTITY_COLLISION_BLOCK']+=1;continue
        values={f:{r[f] for r in records if r[f]} for f in context_fields}
        divergent=[f for f,v in values.items() if len(v)>1]
        if divergent:
            conflicts.update(divergent);unresolved['WITHIN_BLOCK_CONTEXT_CONFLICT']+=1;continue
        context={f:next(iter(v)) if v else '' for f,v in values.items()};contexts[key]=context
        # The source code subtracts Start_treatment, despite its surgery comment.
        start=dates.explicit_date(context['Start_treatment'])
        for endfield,timefield in [('Date_death','Time_to_death_or_last_PET'),('Date_progression','Time_to_progression_or_last_PET')]:
            raw=context[endfield] or context['Last_PET'];end=dates.explicit_date(raw);stored=context[timefield]
            if start is None or end is None or not stored:duration['DURATION_REPLAY_UNAVAILABLE']+=1;continue
            try:number=float(stored)
            except ValueError:duration['STORED_DURATION_UNPARSEABLE']+=1;continue
            if not math.isfinite(number):duration['STORED_DURATION_UNPARSEABLE']+=1;continue
            duration['TREATMENT_ORIGIN_DURATION_MATCH' if abs(number-(end-start).days)<1e-8 else 'TREATMENT_ORIGIN_DURATION_MISMATCH']+=1
    decision=Counter();phase_groups=Counter();events=Counter();coverage=Counter();block_multiplicity=Counter(len(v) for v in patients.values())
    for patient,keys in patients.items():
        if any(k not in contexts for k in keys):decision['UNRESOLVED_OBSERVED_BLOCK_CONTEXT']+=1;continue
        anchors={k:dates.explicit_date(contexts[k]['Date_surgery']) for k in keys}
        if any(v is None for v in anchors.values()):decision['EARLIEST_OBSERVED_TISSUE_DATE_UNRESOLVED']+=1;continue
        first=min(anchors.values());earliest=[k for k,v in anchors.items() if v==first]
        if len(earliest)!=1:decision['TIED_EARLIEST_BLOCKS_NOT_ARBITRARILY_CHOSEN']+=1;continue
        context=contexts[earliest[0]]
        value,reason=endpoint(context,first);decision[reason]+=1
        if value is None:continue
        end,event=value;events['DEATH_EVENT' if event else 'LAST_PET_CENSORED']+=1
        start=dates.explicit_date(context['Start_treatment'])
        phase='INDEX_TREATMENT_DATE_UNRESOLVED' if start is None else 'TISSUE_BEFORE_REPORTED_INDEX_START' if first<start else 'SAME_CALENDAR_DAY_ORDER_UNRESOLVED' if first==start else 'TISSUE_AFTER_REPORTED_INDEX_START'
        prior=dates.core.integer(context['Nr_treatments_before_surgery']);prior_state='PRIOR_COUNT_UNKNOWN' if prior is None else 'NO_PRIOR_TREATMENT_REPORTED' if prior==0 else 'PRIOR_TREATMENT_REPORTED'
        phase_groups[phase,prior_state]+=1
        for field in ['Cancer_Stage','Age_range','Gender','Response_at_3m','Status_at_3m']:
            coverage[field+'_NONMISSING' if context[field] else field+'_MISSING']+=1
    assert sum(decision.values())==len(patients)
    eligible=decision['POSITIVE_POST_TISSUE_ENDPOINT_OPPORTUNITY'];assert sum(events.values())==eligible
    return dict(status='COMPLETE_BLOCK_AND_TISSUE_LANDMARK_OPPORTUNITY_NOT_CLINICAL_ASSOCIATION',source_metadata_rows=len(rows),
        author_patients_with_image_linked_blocks=len(patients),image_linked_patient_block_groups=len(grouped),distinct_images=len(images),
        image_identity_collision_keys=len(bad_images),block_linked_to_multiple_patients=len(bad_blocks),
        patient_observed_block_multiplicity=[dict(observed_blocks=k,author_patients=v) for k,v in sorted(block_multiplicity.items())],
        unresolved_block_reasons=dict(unresolved),within_block_conflicting_field_counts=dict(conflicts),
        accepted_block_contexts=len(contexts),source_row_exclusion_reasons=dict(reasons),date_format_counts={f:dict(v) for f,v in date_formats.items()},
        source_treatment_origin_duration_replay=dict(duration),earliest_tissue_selection_reasons=dict(decision),
        selected_tissue_followup_events=dict(events),selected_metadata_nonmissing_coverage=dict(coverage),
        selected_phase_context_groups=[dict(phase=k[0],prior_treatment=k[1],author_patients=v) for k,v in sorted(phase_groups.items())],
        tissue_landmark_eligible_metadata_upper_bound=eligible,selection_before_endpoint_or_expression_eligibility=True,
        repeated_blocks_not_automatically_patient_conflict=True,no_time_values_or_individual_keys_exported=True,
        no_cell_objects_or_new_clinical_associations_read_or_computed=True,
        not_pretreatment_ICI_prediction_or_CM02_projection=True,cell_capture_and_complete_model_eligibility_pending=True,
        last_PET_censoring_is_author_convention_not_independently_verified_last_alive=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);parser.add_argument('--scripts-root',type=Path,required=True);args=parser.parse_args();assert not args.output.exists()
    path=args.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py';assert hashlib.sha256(path.read_bytes()).hexdigest()=='546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec=importlib.util.spec_from_file_location('bounded_small_clinical_reader',path);transport=importlib.util.module_from_spec(spec);spec.loader.exec_module(transport);transport.PROXY_PORT=int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    rf=transport.RangeFile(None,'https://zenodo.org/api/records/5994136/files/data_for_analysis.zip/content',7123744479);rf.max_chunk=64*1024
    def directory(z):return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC) for x in z.infolist()],separators=(',',':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as z:
            assert directory(z)==DIRECTORY_SHA
            raw=z.read('data_for_analysis_zenodo/protein/clinical_data_protein.csv');assert len(raw)<256*1024 and hashlib.sha256(raw).hexdigest()==CLINICAL_SHA
            reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')));assert set(FIELDS)<=set(reader.fieldnames)
            rows=[{k:row[k] for k in FIELDS} for row in reader];result=summarize(rows)
        with zipfile.ZipFile(rf) as replay:assert directory(replay)==DIRECTORY_SHA
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),clinical_member_SHA256=CLINICAL_SHA,source_directory_SHA256=DIRECTORY_SHA,
            complete_selected_member_CRC_and_fresh_directory_replay=True,wire_bytes=rf.transferred,producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    finally:rf.close()


if __name__=='__main__':main()
