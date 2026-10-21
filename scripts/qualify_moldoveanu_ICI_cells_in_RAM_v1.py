"""Outcome-blind exact source-key/coordinate/type/marker qualification in RAM.

No sample is promoted to an independent patient. Clinical table headers only;
no response values, new tests, state effects or learned models are read.
"""
from pathlib import Path, PurePosixPath
import csv
import gzip
import hashlib
import io
import json
import numpy as np
import pandas as pd

REQUIRED = ['Cell_ID','Sample_ID','Location_Center_X','Location_Center_Y','Ki67']
ST4_REQUIRED = ['obj.id','sample.id','coord.x','coord.y','Cluster','Ki67.prob']


def table_bytes(archive, basename):
    found=[x for x in archive.getmembers() if x.isfile() and PurePosixPath(x.name).name==basename]
    assert len(found)==1
    raw=archive.extractfile(found[0]).read()
    assert len(raw)==found[0].size
    return raw


def decoded(raw):
    return gzip.decompress(raw) if raw.startswith(b'\x1f\x8b') else raw


def header(raw):
    line=raw.split(b'\n',1)[0].decode('utf-8-sig').rstrip('\r')
    delimiter='\t' if '\t' in line else ',' if ',' in line else ' '
    fields=next(csv.reader([line],delimiter=delimiter,skipinitialspace=True))
    return fields,delimiter


def bind(st4, raw):
    assert set(ST4_REQUIRED)<=set(st4) and set(REQUIRED)<=set(raw)
    left=st4.copy();right=raw.copy()
    assert not left[[x for x in ST4_REQUIRED if x!='Ki67.prob']].isna().any().any()
    assert not right[REQUIRED].isna().any().any()
    assert not left.duplicated(['sample.id','obj.id']).any()
    assert not right.duplicated(['Sample_ID','Cell_ID']).any()
    samples=set(right['Sample_ID'])
    left=left.loc[left['sample.id'].isin(samples)]
    joined=left.merge(right,left_on=['sample.id','obj.id'],right_on=['Sample_ID','Cell_ID'],
        how='left',validate='one_to_one',indicator=True)
    assert len(joined)>0 and joined['_merge'].eq('both').all()
    assert np.isfinite(joined[['coord.x','coord.y','Location_Center_X','Location_Center_Y','Ki67']].to_numpy(dtype=float)).all()
    a=joined[['coord.x','coord.y']].to_numpy(dtype=float)
    b=joined[['Location_Center_X','Location_Center_Y']].to_numpy(dtype=float)
    maximum=float(np.max(np.abs(a-b)))
    assert maximum<=1e-6
    assert joined['Ki67.prob'].dropna().between(0,1).all()
    assert all(isinstance(v,str) and v for v in joined['Cluster'])
    return joined,dict(all_selected_author_rows_exactly_joined=True,
        exact_composite_keys_unique_in_both_tables=True,
        coordinate_maximum_absolute_difference=maximum,
        fixed_coordinate_serialization_tolerance=1e-6,
        source_samples=len(samples),selected_source_samples=joined['sample.id'].nunique(),
        all_raw_ICI_cells=len(raw),author_selected_ICI_cells=len(joined),
        author_Ki67_probability_missing_not_filled=int(joined['Ki67.prob'].isna().sum()),
        raw_unselected_cells_remain_excluded=len(raw)-len(joined))


def analyze(bundle, parameters):
    assert parameters=={'coordinate_tolerance':1e-6,'minimum_parent_cells':20}
    archive=bundle['archive'];proof=bundle['qualification']
    assert proof['archive_MD5_verified'] and proof['full_gzip_CRC_checked']
    raw_bytes=table_bytes(archive,'ICI_meanIntensity.tsv.gz')
    content=decoded(raw_bytes);fields,delimiter=header(content)
    assert delimiter=='\t' and set(REQUIRED)<=set(fields)
    selected=REQUIRED+[x for x in ['FoxP3','CD3','CD4','CD8a','Histone.H3','Area','Cell_Area','OX40','LAG3','PD1','CD127'] if x in fields]
    raw=pd.read_csv(io.BytesIO(content),sep='\t',usecols=selected,dtype={'Cell_ID':str,'Sample_ID':str})
    # Full scalar CSV replay of selected numerical fields and exact keys.
    ordinary_keys=set();count=0;total=0.;compensation=0.
    for row in csv.DictReader(io.StringIO(content.decode('utf-8-sig')),delimiter='\t'):
        key=(row['Sample_ID'],row['Cell_ID']);assert key not in ordinary_keys;ordinary_keys.add(key)
        value=float(row['Ki67']);assert np.isfinite(value)
        adjusted=value-compensation;updated=total+adjusted;compensation=(updated-total)-adjusted;total=updated;count+=1
    assert count==len(raw) and ordinary_keys==set(zip(raw['Sample_ID'],raw['Cell_ID']))
    assert abs(float(raw['Ki67'].sum())-total)<=1e-9*max(1.,abs(total))
    st4_bytes=table_bytes(archive,'ST4_cell_data.txt');st4_content=decoded(st4_bytes)
    st4_fields,separator=header(st4_content)
    st4=pd.read_csv(io.BytesIO(st4_content),sep=separator,dtype={'obj.id':str,'sample.id':str},skipinitialspace=True)
    joined,qualification=bind(st4,raw)
    clinical_bytes=table_bytes(archive,'ST1_sample_info.txt');clinical_fields,_=header(decoded(clinical_bytes))
    roles={'author_regulatory':['Treg'],'author_conventional_CD4':['Th.ae','Th.naive'],
        'author_CD8':['Tc.ae','Tc.naive']}
    role_coverage=[]
    for role,labels in roles.items():
        subset=joined.loc[joined['Cluster'].isin(labels)]
        counts=subset.groupby('sample.id').size()
        role_coverage.append(dict(role=role,source_labels=labels,source_cells=len(subset),
            source_samples_present=len(counts),source_samples_parent_ge20=int((counts>=20).sum())))
    # Accepted objects remain only in the existing bounded process RAM.
    archive._own_ICI_qualified_cells=joined
    archive._own_ICI_qualification=dict(qualification)
    return dict(status='COMPLETE_OUTCOME_BLIND_ICI_SOURCE_CELL_QUALIFICATION_NOT_PATIENT_OR_EFFECT_ACCEPTANCE',
        qualification=qualification,ICI_field_names=fields,author_cell_table_field_names=st4_fields,
        clinical_header_field_names=clinical_fields,role_coverage=role_coverage,
        actual_author_type_names=sorted(joined['Cluster'].unique()),
        independent_full_CSV_key_count_and_Ki67_scalar_replay=True,
        source_member_SHA256=dict(ICI=hashlib.sha256(raw_bytes).hexdigest(),ST4=hashlib.sha256(st4_bytes).hexdigest(),ST1=hashlib.sha256(clinical_bytes).hexdigest()),
        no_clinical_values_read_or_patient_independence_assumed=True,
        no_new_state_tests_or_individual_outputs=True,
        producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
