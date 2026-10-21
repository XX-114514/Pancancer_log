"""Outcome-blind patient/phase availability from fixed small public metadata.

Metadata-only upper bounds before cell capture/QC. No clinical effects, patient
rows, raw matrices, source model objects or large processed SCE are exported.
Dates accept explicit ISO only; ambiguous dates remain unknown.
"""
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import zipfile

DIRECTORY_SHA='dff220e35d38da2b60c42950190d71b2259f12aae3c21ad4bc6c05ad5c5d2893'
CLINICAL_SHA='ac4a02e4d0e23504aa980b21585de74c4d7390d2e6e0ec2399303c5c29ca5b25'
PANEL_SHA='3133cfc5c21314379a577fc6f566f61e0c9d934987bf1f970a4435f7d79a3e3a'
MISSING={'','na','nan','n/a','none','null'}


def missing(value):return value is None or str(value).strip().lower() in MISSING


def iso_date(value):
    if missing(value):return None
    text=str(value).strip()
    if not re.fullmatch(r'\d{4}-\d{2}-\d{2}(?:[ T]\d{2}:\d{2}:\d{2})?',text):return None
    try:return date.fromisoformat(text[:10])
    except ValueError:return None


def integer(value):
    if missing(value):return None
    try:
        number=float(value)
        return int(number) if math.isfinite(number) and number>=0 and number.is_integer() else None
    except ValueError:return None


def boolean(value):
    value=str(value).strip().lower()
    if value in {'0','false','no','n','nein'}:return False
    if value in {'1','true','yes','y','ja'}:return True
    return None


def regimen(value):
    text=str(value).lower();pd1=bool(re.search(r'nivolumab|pembrolizumab|\banti[- _]?pd[- _]?1\b',text))
    ctla4=bool(re.search(r'ipilimumab|\banti[- _]?ctla[- _]?4\b',text))
    if pd1 and ctla4:return 'PD1_and_CTLA4_text_present_not_schedule_equivalence'
    if pd1:return 'PD1_text_present'
    if ctla4:return 'CTLA4_text_present'
    return 'OTHER_OR_UNCLASSIFIED_TEXT'


def summarize(rows):
    needed={'PatientID','ImageNumber','Date_surgery','Start_treatment','Drug_after_surgery',
        'Treatment_after_surgery','Adjuvant','Nr_treatments_before_surgery','Cancer_Stage','Response_at_3m'}
    assert rows and needed<=set(rows[0])
    grouped={};image_patients={};no_patient=0
    for row in rows:
        if missing(row['PatientID']):no_patient+=1;continue
        identity=row['PatientID'].strip();grouped.setdefault(identity,[]).append(row)
        if not missing(row['ImageNumber']):image_patients.setdefault(row['ImageNumber'].strip(),set()).add(identity)
    fields=sorted(needed-{'PatientID','ImageNumber','Response_at_3m'})
    reasons=Counter();groups=Counter();dates={field:Counter() for field in ['Date_surgery','Start_treatment']}
    response_counts=Counter();complete_context=0;inconsistent=0
    for records in grouped.values():
        values={field:{r[field].strip() for r in records if not missing(r[field])} for field in fields}
        if any(len(v)>1 for v in values.values()):
            inconsistent+=1;reasons['conflicting_context_within_author_patient']+=1;continue
        context={field:next(iter(values[field])) if values[field] else None for field in fields}
        complete_context+=1
        for field in dates:
            value=context[field]
            dates[field]['missing' if missing(value) else 'explicit_ISO_parsed' if iso_date(value) else 'unparsed_nonISO_or_invalid']+=1
        surgery=iso_date(context['Date_surgery']);start=iso_date(context['Start_treatment'])
        if surgery is None or start is None:phase='DATE_ORDER_UNRESOLVED'
        elif surgery<start:phase='SURGERY_STRICTLY_BEFORE_REPORTED_AFTER_SURGERY_TREATMENT'
        elif surgery==start:phase='SAME_DAY_ORDER_UNRESOLVED'
        else:phase='SURGERY_AFTER_REPORTED_TREATMENT_START'
        prior=integer(context['Nr_treatments_before_surgery'])
        prior_state='NO_PRIOR_TREATMENTS_REPORTED' if prior==0 else 'PRIOR_TREATMENTS_REPORTED' if prior is not None else 'PRIOR_TREATMENTS_UNKNOWN'
        adjuvant=boolean(context['Adjuvant']);adj_state='ADJUVANT' if adjuvant is True else 'NOT_ADJUVANT_CODE' if adjuvant is False else 'ADJUVANT_UNKNOWN'
        drug=regimen(context['Drug_after_surgery'])
        groups[(phase,prior_state,adj_state,drug)]+=1
        outcome_values={r['Response_at_3m'].strip() for r in records if not missing(r['Response_at_3m'])}
        response_counts['conflicting_nonmissing_strings' if len(outcome_values)>1 else 'available_nonmissing_string' if outcome_values else 'missing']+=1
        if phase=='SURGERY_STRICTLY_BEFORE_REPORTED_AFTER_SURGERY_TREATMENT' and prior==0 and adjuvant is False and drug!='OTHER_OR_UNCLASSIFIED_TEXT':
            reasons['strict_baseline_nonadjuvant_ICI_metadata_upper_bound']+=1
            if len(outcome_values)==1:reasons['strict_baseline_nonadjuvant_ICI_with_concordant_nonmissing_3m_response_upper_bound']+=1
    return dict(status='COMPLETE_AUTHOR_PATIENT_PHASE_AVAILABILITY_NOT_CLINICAL_OR_CELL_QUALIFICATION',
        metadata_records=len(rows),author_patient_keys=len(grouped),rows_without_patient_key=no_patient,
        distinct_nonmissing_image_keys=len(image_patients),image_keys_linked_to_multiple_patients=sum(len(v)>1 for v in image_patients.values()),
        patients_without_observed_context_conflict=complete_context,patients_with_context_conflict=inconsistent,
        source_date_format_counts={k:dict(v) for k,v in dates.items()},
        three_month_response_availability_only=dict(response_counts),reason_counts=dict(reasons),
        phase_context_groups=[dict(phase=k[0],prior_treatment=k[1],adjuvant=k[2],drug_text_group=k[3],author_patients=n) for k,n in sorted(groups.items())],
        metadata_only_upper_bounds_before_cell_capture_QC_or_complete_clinical_design=True,
        response_values_or_patient_rows_not_exported=True,
        coded_adjuvant_and_drug_text_interpretations_require_semantic_review_before_model=True,
        no_baseline_age_or_response_criteria_equivalence_assumed=True,
        no_clinical_association_or_patient_expression_object_read=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--scripts-root',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists()
    path=a.scripts_root/'melanoma_IMC_existing_clash_transport_v1.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='546b491e20e6797e5e255bfa226fa568a83dca08eb01334a27128fb0587755ac'
    spec=importlib.util.spec_from_file_location('small_public_meta_proxy',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.PROXY_PORT=int(os.environ['PANCANCER_EXISTING_CLASH_PORT'])
    rf=m.RangeFile(None,'https://zenodo.org/api/records/5994136/files/data_for_analysis.zip/content',7123744479);rf.max_chunk=64*1024
    def directory(archive):return hashlib.sha256(json.dumps([(x.filename,x.file_size,x.compress_size,x.CRC) for x in archive.infolist()],separators=(',',':')).encode()).hexdigest()
    try:
        with zipfile.ZipFile(rf) as archive:
            assert directory(archive)==DIRECTORY_SHA
            clinical=archive.read('data_for_analysis_zenodo/protein/clinical_data_protein.csv');assert hashlib.sha256(clinical).hexdigest()==CLINICAL_SHA
            panel=archive.read('data_for_analysis_zenodo/protein/melanoma_1.06_protein.csv');assert hashlib.sha256(panel).hexdigest()==PANEL_SHA
            rows=list(csv.DictReader(io.StringIO(clinical.decode('utf-8-sig'))));result=summarize(rows)
            panel_rows=list(csv.DictReader(io.StringIO(panel.decode('utf-8-sig'))))
            markers=[]
            for row in panel_rows:
                target=row.get('Target','');normalized=re.sub('[^a-z0-9]','',target.lower())
                if normalized in {'cd3','cd4','cd8a','cd8','foxp3','ki67','tcf1tcf7'} or normalized.startswith('ki67'):
                    markers.append(dict(public_target=target,metal_tag=row.get('Metal Tag'),clean_target=row.get('clean_target')))
            result['declared_panel_targets']=markers
        with zipfile.ZipFile(rf) as replay:assert directory(replay)==DIRECTORY_SHA
        result.update(completed_UTC=datetime.now(timezone.utc).isoformat(),clinical_member_SHA256=CLINICAL_SHA,
            panel_member_SHA256=PANEL_SHA,source_directory_SHA256=DIRECTORY_SHA,selected_full_member_CRCs_checked=True,
            fresh_wire_directory_replay=True,wire_bytes=rf.transferred,
            producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    finally:rf.close()


if __name__=='__main__':main()
