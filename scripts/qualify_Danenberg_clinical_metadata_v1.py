"""Small author-reported follow-up availability only; no phenotype or model.

Patient keys and records are RAM only. Unknown event/validation codes remain
unknown; no positive-time offset, subgroup selection or endpoint substitution.
Metadata opportunities are upper bounds before cell capture qualification.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path,PurePosixPath
import argparse,csv,hashlib,importlib.util,io,json,math,os,zipfile,time

DIRECTORY_SHA='17911e718220292268a8cddfac32fb4bfa6c021a615ead0e011464c4a3658814'
FIELDS=['metabric_id','ERStatus','LymphNodesOrdinal','sizeOrdinal','Grade',
        'ERBB2_pos','yearsToStatus','DeathBreast','isValidation','PAM50','IntClust']
MISSING={'','na','nan','n/a','none','null'}


def clean(v):return '' if v is None or str(v).strip().lower() in MISSING else str(v).strip()


def binary(v):
    text=clean(v).lower()
    if text in {'true','t','1','1.0'}:return True
    if text in {'false','f','0','0.0'}:return False
    return None


def positive_time(v):
    try:x=float(clean(v))
    except ValueError:return None
    return x if math.isfinite(x) and x>0 else None


def summarize(rows):
    assert rows and set(FIELDS)<=set(rows[0]);groups=defaultdict(list);excluded=Counter()
    for row in rows:
        r={f:clean(row[f]) for f in FIELDS}
        if not r['metabric_id']:excluded['MISSING_PATIENT_KEY_ROWS']+=1;continue
        groups[r['metabric_id']].append(r)
    decisions=Counter();strata=defaultdict(Counter);coverage=Counter();duplicate=0;disagreement=0
    for records in groups.values():
        duplicate+=len(records)-1
        if len({tuple(r[f] for f in FIELDS) for r in records})>1:
            disagreement+=1;decisions['DUPLICATE_PATIENT_CONTEXT_UNRESOLVED']+=1;continue
        r=records[0]
        for f in FIELDS[1:]:coverage[f+'_NONMISSING' if r[f] else f+'_MISSING']+=1
        event=binary(r['DeathBreast']);validation=binary(r['isValidation']);followup=positive_time(r['yearsToStatus'])
        if followup is None:decisions['NO_POSITIVE_FINITE_REPORTED_FOLLOWUP']+=1;continue
        if event is None:decisions['UNKNOWN_REPORTED_EVENT_CODE']+=1;continue
        if validation is None:decisions['UNKNOWN_AUTHOR_VALIDATION_FLAG']+=1;continue
        decisions['REPORTED_FOLLOWUP_EVENT_AND_FLAG_METADATA_OPPORTUNITY']+=1
        key='AUTHOR_VALIDATION_TRUE' if validation else 'AUTHOR_VALIDATION_FALSE';bucket=strata[key]
        bucket['metadata_patients']+=1;bucket['reported_breast_death_events' if event else 'reported_non_event_records']+=1
        cov=['ERStatus','LymphNodesOrdinal','sizeOrdinal','Grade','ERBB2_pos']
        if all(r[f] for f in cov):bucket['five_reported_baseline_covariates_nonmissing']+=1
        if all(r[f] for f in cov) and binary(r['ERBB2_pos']) is not None:bucket['five_covariates_nonmissing_and_HER2_code_recognized']+=1
    assert sum(decisions.values())==len(groups)
    return dict(status='COMPLETE_REPORTED_FOLLOWUP_METADATA_OPPORTUNITIES_NOT_CLINICAL_EFFECT',source_metadata_rows=len(rows),
        distinct_nonmissing_author_patient_keys=len(groups),duplicate_patient_rows=duplicate,discordant_patient_groups=disagreement,
        excluded_row_reasons=dict(excluded),patient_qualification_reasons=dict(decisions),field_nonmissing_coverage=dict(coverage),
        by_author_validation_flag=[dict(author_flag=k,**v) for k,v in sorted(strata.items())],
        no_patient_rows_or_times_or_phenotypes_exported=True,no_new_survival_or_state_association=True,
        author_flag_not_automatically_independent_cohort=True,time_origin_and_treatment_history_pending_primary_review=True,
        no_CM02_projection_or_ICI_response_equivalence=True,cell_state_capture_and_complete_model_eligibility_pending=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--scripts-root',type=Path,required=True);a=p.parse_args();assert not a.output.exists()
    path=a.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py';assert hashlib.sha256(path.read_bytes()).hexdigest()=='546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec=importlib.util.spec_from_file_location('bounded_small_breast_metadata',path);t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t);t.PROXY_PORT=int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    class Limited(t.RangeFile):
        def _read_chunk(self,n):
            assert self.requested+n<=4*1024**2 and time.monotonic()-self.started<=180
            return super()._read_chunk(n)
    rf=Limited(None,'https://zenodo.org/api/records/7324285/files/MBTMEStrPublicProcessedDataCode.zip/content',910666172);rf.max_chunk=64*1024
    def directory(z):return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC) for x in z.infolist()],separators=(',',':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as z:
            assert directory(z)==DIRECTORY_SHA
            selected=[x for x in z.infolist() if PurePosixPath(x.filename).name=='IMCClinical.csv'];assert len(selected)==1 and selected[0].file_size==51316 and selected[0].CRC==1165356369
            raw=z.read(selected[0]);reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')));assert reader.fieldnames==FIELDS
            rows=[{k:r[k] for k in FIELDS} for r in reader];result=summarize(rows)
            # Separate pandas grouping/conversion replay, no individual outputs.
            import pandas as pd
            frame=pd.DataFrame(rows)
            for column in FIELDS:
                frame[column]=frame[column].str.strip()
                frame[column]=frame[column].mask(frame[column].str.lower().isin(MISSING),'')
            frame=frame.loc[frame.metabric_id.ne('')]
            assert frame.metabric_id.nunique()==result['distinct_nonmissing_author_patient_keys']
            frame=frame.drop_duplicates(FIELDS)
            disputed=set(frame.loc[frame.metabric_id.duplicated(keep=False),'metabric_id'])
            assert len(disputed)==result['discordant_patient_groups']
            frame=frame.loc[~frame.metabric_id.isin(disputed)]
            mapping={'true':True,'t':True,'1':True,'1.0':True,'false':False,'f':False,'0':False,'0.0':False}
            events=frame.DeathBreast.str.lower().map(mapping);flags=frame.isValidation.str.lower().map(mapping);followup=pd.to_numeric(frame.yearsToStatus,errors='coerce')
            good=followup.gt(0)&followup.map(math.isfinite)&events.notna()&flags.notna()
            assert int(good.sum())==result['patient_qualification_reasons'].get('REPORTED_FOLLOWUP_EVENT_AND_FLAG_METADATA_OPPORTUNITY',0)
            for r in result['by_author_validation_flag']:
                mask=good&flags.eq(r['author_flag']=='AUTHOR_VALIDATION_TRUE')
                assert int(mask.sum())==r['metadata_patients'] and int(events.loc[mask].eq(True).sum())==r.get('reported_breast_death_events',0)
        with zipfile.ZipFile(rf) as replay:assert directory(replay)==DIRECTORY_SHA
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),clinical_member_SHA256=hashlib.sha256(raw).hexdigest(),full_small_member_CRC_and_fresh_directory_replay=True,
            independent_group_count_and_event_conversion_replay=True,wire_bytes=rf.transferred,producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    finally:rf.close()


if __name__=='__main__':main()
