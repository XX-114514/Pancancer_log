"""Honor the author's first-day monthly export; never infer hidden day values.

The source prepares treatment-origin day durations before formatting public
dates as 01-%b-20%y. Month-interval compatibility is weaker than an independent
daily replay. Same exported date is same month, not verified same-day ordering.
No century correction, exact tissue-origin duration, model or cell read.
"""
from collections import Counter,defaultdict
from datetime import date,datetime,timezone
from pathlib import Path
import calendar,argparse,csv,hashlib,importlib.util,io,json,os,zipfile
import qualify_Hoch_block_landmark_scope_v1 as previous

FIELDS=previous.FIELDS+['censoring_death','censoring_progression']


def bounds(value):
    d=previous.dates.explicit_date(value)
    if d is None:return None
    if d.day!=1:raise ValueError('NONFIRST_DAY_NOT_COMPATIBLE_WITH_FIXED_AUTHOR_EXPORT')
    return d,date(d.year,d.month,calendar.monthrange(d.year,d.month)[1])


def duration_check(context,endfield,timefield):
    start=bounds(context['Start_treatment']);end=bounds(context[endfield] or context['Last_PET'])
    if start is None or end is None or not context[timefield]:return 'UNAVAILABLE_REPORTED_TIME_OR_MONTHS'
    try:value=float(context[timefield])
    except ValueError:return 'UNPARSEABLE_REPORTED_TIME'
    if not value.is_integer():return 'NONINTEGER_REPORTED_DAY_DURATION'
    low=(end[0]-start[1]).days;high=(end[1]-start[0]).days
    return 'COMPATIBLE_WITH_EXPORTED_MONTH_INTERVALS' if low<=value<=high else 'OUTSIDE_EXPORTED_MONTH_INTERVALS'


def summarize(rows):
    assert rows and set(FIELDS)<=set(rows[0])
    old=previous.summarize(rows);groups=defaultdict(list);patients=defaultdict(list);formats=Counter()
    for row in rows:
        r={k:previous.canonical(row.get(k)) for k in FIELDS}
        for field in previous.DATE_FIELDS:
            d=previous.dates.explicit_date(r[field])
            formats['MISSING' if not r[field] else 'UNPARSED' if d is None else 'FIRST_DAY_MONTH_EXPORT' if d.day==1 else 'OTHER_DAY']+=1
        if all(r[k] for k in ['PatientID','BlockID','ImageNumber']):groups[r['PatientID'],r['BlockID']].append(r)
    assert old['image_identity_collision_keys']==old['block_linked_to_multiple_patients']==0
    contexts={};conflicts=0;time_checks=Counter();events=Counter()
    for key,records in groups.items():
        patients[key[0]].append(key)
        values={f:{r[f] for r in records if r[f]} for f in FIELDS if f not in ['PatientID','BlockID','ImageNumber']}
        if any(len(v)>1 for v in values.values()):conflicts+=1;continue
        c={f:next(iter(v)) if v else '' for f,v in values.items()};contexts[key]=c
        for endfield,timefield,eventfield in [('Date_death','Time_to_death_or_last_PET','censoring_death'),('Date_progression','Time_to_progression_or_last_PET','censoring_progression')]:
            time_checks[endfield+'__'+duration_check(c,endfield,timefield)]+=1
            flag=previous.dates.core.boolean(c[eventfield]);expected=bool(c[endfield])
            events[eventfield+'__'+('MISSING_OR_UNPARSED_FLAG' if flag is None else 'PRESENCE_AGREEMENT' if flag==expected else 'PRESENCE_MISMATCH')]+=1
    selected=Counter();phase=Counter();monthly_future=Counter();preindex=Counter()
    for patient,keys in patients.items():
        if any(k not in contexts for k in keys):selected['WITHIN_BLOCK_CONTEXT_UNRESOLVED']+=1;continue
        anchors={k:bounds(contexts[k]['Date_surgery']) for k in keys}
        if any(v is None for v in anchors.values()):selected['EARLIEST_OBSERVED_MONTH_UNRESOLVED']+=1;continue
        first=min(v[0] for v in anchors.values());ties=[k for k,v in anchors.items() if v[0]==first]
        if len(ties)!=1:selected['MULTIPLE_BLOCKS_IN_EARLIEST_MONTH_UNRESOLVED']+=1;continue
        selected['UNIQUE_EARLIEST_OBSERVED_BLOCK_MONTH']+=1;c=contexts[ties[0]];surgery=anchors[ties[0]]
        start=bounds(c['Start_treatment']);end=bounds(c['Date_death'] or c['Last_PET'])
        context_phase='INDEX_MONTH_UNKNOWN' if start is None else 'TISSUE_MONTH_STRICTLY_BEFORE_INDEX_MONTH' if surgery[1]<start[0] else 'TISSUE_MONTH_STRICTLY_AFTER_INDEX_MONTH' if surgery[0]>start[1] else 'SAME_MONTH_ORDER_UNRESOLVED'
        phase[context_phase]+=1
        endpoint_phase='END_MONTH_UNKNOWN' if end is None else 'END_MONTH_STRICTLY_AFTER_TISSUE_MONTH' if end[0]>surgery[1] else 'END_MONTH_STRICTLY_BEFORE_TISSUE_MONTH' if end[1]<surgery[0] else 'END_AND_TISSUE_SAME_MONTH_ORDER_UNRESOLVED'
        monthly_future[endpoint_phase]+=1
        if context_phase!='TISSUE_MONTH_STRICTLY_BEFORE_INDEX_MONTH':continue
        reason=duration_check(c,'Date_death','Time_to_death_or_last_PET')
        if reason!='COMPATIBLE_WITH_EXPORTED_MONTH_INTERVALS':preindex[reason]+=1;continue
        if float(c['Time_to_death_or_last_PET'])<=0:preindex['NONPOSITIVE_REPORTED_SURVIVAL']+=1;continue
        flag=previous.dates.core.boolean(c['censoring_death'])
        if flag is None or flag!=bool(c['Date_death']):preindex['EVENT_FLAG_NOT_ACCEPTED']+=1;continue
        preindex['MONTH_PREINDEX_REPORTED_OS_METADATA_OPPORTUNITY']+=1
        preindex['DEATH_EVENTS' if flag else 'AUTHOR_LAST_PET_CENSORED']+=1
    assert sum(selected.values())==len(patients)
    return dict(status='COMPLETE_AUTHOR_MONTH_PRECISION_AND_REPORTED_DURATION_COMPATIBILITY',
        patients_with_image_linked_blocks=len(patients),observed_patient_block_groups=len(groups),accepted_block_contexts=len(contexts),within_block_conflicts=conflicts,
        date_field_representation_counts=dict(formats),source_reported_duration_month_interval_checks=dict(time_checks),source_event_flag_checks=dict(events),
        earliest_month_selection_reasons=dict(selected),selected_tissue_index_month_relationship=dict(phase),selected_tissue_endpoint_month_relationship=dict(monthly_future),
        preindex_reported_OS_metadata_only_opportunities=dict(preindex),
        source_date_export_format='01-%b-20%y',daily_duration_precedes_monthly_date_export=True,
        previous_exact_day_replay_and_same_day_phase_not_accepted=True,monthly_compatibility_not_exact_daily_verification=True,
        no_hidden_day_or_original_century_reconstruction=True,no_exact_tissue_origin_survival_reconstructed=True,
        metadata_not_cell_capture_or_full_clinical_model_eligibility=True,no_new_P_or_clinical_association_or_CM02_projection=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--scripts-root',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    path=a.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py';assert hashlib.sha256(path.read_bytes()).hexdigest()=='546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec=importlib.util.spec_from_file_location('bounded_month_scope_reader',path);t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);t.PROXY_PORT=int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    rf=t.RangeFile(None,'https://zenodo.org/api/records/5994136/files/data_for_analysis.zip/content',7123744479);rf.max_chunk=64*1024
    def directory(z):return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC) for x in z.infolist()],separators=(',',':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as z:
            assert directory(z)==previous.DIRECTORY_SHA
            raw=z.read('data_for_analysis_zenodo/protein/clinical_data_protein.csv');assert len(raw)<256*1024 and hashlib.sha256(raw).hexdigest()==previous.CLINICAL_SHA
            reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')));assert set(FIELDS)<=set(reader.fieldnames)
            rows=[{k:r[k] for k in FIELDS} for r in reader];result=summarize(rows)
        with zipfile.ZipFile(rf) as replay:assert directory(replay)==previous.DIRECTORY_SHA
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),clinical_member_SHA256=previous.CLINICAL_SHA,full_small_member_CRC_and_fresh_directory_replay=True,wire_bytes=rf.transferred,producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    finally:rf.close()


if __name__=='__main__':main()
