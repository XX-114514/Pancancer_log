"""Read an existing official workbook; return source aggregates, never rows.

Expression/frequency sheets contribute only their patient key column. No
derived cellular frequencies, clinical effects, thresholds or models are read.
"""
from collections import Counter
from pathlib import Path
import hashlib
import math
import re
import openpyxl

EXPECTED_SHA='f33f58f219adcecfc5355b1ce8b526dfd36d7eca4555a320ba1deb185d122982'


def key(value):
    if isinstance(value,bool): return None
    if isinstance(value,int) and value>0: return str(value)
    if isinstance(value,float) and math.isfinite(value) and value>0 and value.is_integer(): return str(int(value))
    if isinstance(value,str) and re.fullmatch(r'[1-9][0-9]*',value.strip()): return value.strip()
    return None


def numeric(value):
    if isinstance(value,bool) or value is None: return None
    try: x=float(value)
    except (ValueError,TypeError): return None
    return x if math.isfinite(x) else None


def category(value,allowed):
    x=str(value).strip() if value is not None else ''
    return x if x in allowed else 'MISSING_OR_UNRECOGNIZED'


def summarize(rows):
    response=Counter(category(r['RECIST'],['CR','PR','SD','PD']) for r in rows)
    complete={}
    fields=['Age','Sex','Smoking (pack years)','Was chemo given before PD1-blockade?',
        'PD1-blockade given as line of therapy #','PD1-blockade agent','Clinical PDL1 IHC (%E1L3N+ of tumor cells)']
    for name in fields:
        complete[name]=sum(r.get(name) is not None and str(r[name]).strip() not in ['', 'NA','N/A','Unknown','unknown'] for r in rows)
    survival={}
    for name,duration,event in [('PFS','PFS (days)','PFS censored (0=yes)'),('OS','OS (days)','OS Censored (0=not dead)')]:
        valid=[r for r in rows if numeric(r.get(duration)) is not None and numeric(r[duration])>0 and numeric(r.get(event)) in [0,1]]
        survival[name]=dict(positive_reported_duration_and_binary_flag=len(valid),events_code1=sum(numeric(r[event])==1 for r in valid),
            missing_nonpositive_or_unrecognized=len(rows)-len(valid),time_origin_not_qualified=True)
    return dict(patients=len(rows),RECIST_counts=dict(response),RECIST_evaluable=sum(response[v] for v in ['CR','PR','SD','PD']),
        CR_PR=sum(response[v] for v in ['CR','PR']),SD_PD=sum(response[v] for v in ['SD','PD']),
        field_nonblank_nonliteralNA_counts=complete,reported_survival_metadata=survival,
        no_sample_level_or_clinical_effects_returned=True)


def qualify(path):
    path=Path(path);digest=hashlib.sha256(path.read_bytes()).hexdigest();assert digest==EXPECTED_SHA
    book=openpyxl.load_workbook(path,read_only=True,data_only=True)
    stream=iter(book['B_Patient_Metadata'].values);header=next(stream)
    assert len(header)==len(set(header)) and header[0]=='Patient #'
    cases={};unparsed=0
    for raw in stream:
        k=key(raw[0])
        if k is None:
            unparsed+=int(any(v is not None for v in raw));continue
        assert k not in cases,'DUPLICATE_CLINICAL_PATIENT_KEY'
        cases[k]=dict(zip(header,raw))
    modality={};keys={}
    for name,label in [('C_ISH_Panel_Freq_of_DAPI','co_registered_RNA'),('D_Ab_Panel_Freq_of_DAPI','co_registered_protein')]:
        # Only the ID column; never read the published derived frequencies.
        stream=iter(book[name].iter_rows(min_col=1,max_col=1,values_only=True));assert next(stream)[0]=='Patient #'
        rows=[r[0] for r in stream if r[0] is not None]
        kk=[key(r) for r in rows]
        good=[k for k in kk if k is not None]
        assert len(good)==len(set(good)),'DUPLICATE_MODALITY_PATIENT_KEY'
        ids=set(good);keys[label]=ids
        modality[label]=dict(published_patient_rows=len(rows),strictly_parsed_unique_patient_keys=len(ids),unparsed_key_rows=sum(k is None for k in kk),
            missing_author_clinical_key=len(ids-set(cases)),clinical=summarize([cases[k] for k in sorted(ids&set(cases))]))
    both=keys['co_registered_RNA']&keys['co_registered_protein']
    # Generic staining protocol contains no human observations or clinical rows.
    protocol=[list(r) for r in book['P_Stain_Protocol'].values if any(v is not None for v in r)]
    assert len(protocol)<=100
    result=dict(status='COMPLETE_EXISTING_WORKBOOK_SOURCE_CLINICAL_MODALITY_METADATA_ONLY',source='Chen_NatImmunol2024_Zenodo11198494',
        original_workbook_SHA256=digest,clinical_sheet_nonpatient_note_rows=unparsed,
        all_author_clinical=summarize(list(cases.values())),modalities=modality,
        RNA_protein_shared_published_patient_keys=len(both),RNA_only_patient_keys=len(keys['co_registered_RNA']-keys['co_registered_protein']),
        protein_only_patient_keys=len(keys['co_registered_protein']-keys['co_registered_RNA']),
        both_modalities_clinical=summarize([cases[k] for k in sorted(both&set(cases))]),
        staining_protocol=protocol,individual_specimen_phase_and_survival_time_origin_not_yet_verified=True,
        published_patient_tables_not_full_raw_file_patient_key_verification=True,
        no_cell_frequency_expression_coordinate_or_clinical_effects_read=True,
        no_workbook_or_individual_rows_written=True)
    book.close();return result
