"""Nonsequence TCR state measurement and exact RNA join; all records RAM only.

Raw CDR3/clonotype/receptor gene fields are never selected by the CSV parser.
Source-reported sizes are not independent sequence-based clone reconstruction.
"""
from collections import Counter
from pathlib import Path
import gzip
import hashlib
import io
import json
import re
import resource
import numpy as np
import pandas as pd

SELECTED=['sampleID','cellID','sub_cell_type','expansion','clonotype_number','T_new_name']
MISSING={'','NA','NaN','nan','None','none','N/A'}
FLAG={'expanded':True,'expansion':True,'yes':True,'true':True,'nonexpanded':False,
      'notexpanded':False,'nonexpansion':False,'noexpansion':False,'no':False,'false':False,'singleton':False,'single':False}


def parse_selected(compressed,plan):
    assert len(compressed)==plan['source']['bytes'] and compressed[:3]==b'\x1f\x8b\x08'
    stated=int.from_bytes(compressed[-4:],'little');assert 0<stated<=plan['max_uncompressed_bytes']
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as f:raw=f.read(plan['max_uncompressed_bytes']+1)
    assert len(raw)==stated and len(raw)<=plan['max_uncompressed_bytes']
    header=raw.split(b'\n',1)[0].decode('utf-8-sig').strip('\r')
    import csv
    assert next(csv.reader([header]))==plan['expected_TCR_header']
    frame=pd.read_csv(io.BytesIO(raw),usecols=SELECTED,dtype=str,keep_default_na=False)
    del raw
    assert list(frame.columns)==SELECTED and 0<len(frame)<=plan['max_TCR_rows']
    return frame


def expansion_qualification(frame):
    raw=frame.clonotype_number
    number=pd.to_numeric(raw,errors='coerce').to_numpy(float)
    valid=np.isfinite(number)&(number>=1)&(number==np.floor(number))&(number<=len(frame))
    explicit_missing=(raw.isin(MISSING)|raw.isin(['0','0.0'])).to_numpy()
    e_numeric=pd.to_numeric(frame.expansion,errors='coerce').to_numpy(float)
    normalized=frame.expansion.str.lower().str.replace(r'[^a-z]','',regex=True)
    flags=normalized.map(FLAG).to_numpy(object).copy();known=np.array([isinstance(v,(bool,np.bool_)) for v in flags])
    binary=np.isfinite(e_numeric)&np.isin(e_numeric,[0.,1.])
    flags[binary]=e_numeric[binary].astype(bool);known|=binary
    reported_size_equals_expansion_numeric=valid&np.isfinite(e_numeric)&(number==e_numeric)
    flags[reported_size_equals_expansion_numeric]=number[reported_size_equals_expansion_numeric]>=2
    known|=reported_size_equals_expansion_numeric
    flag_values=np.array([bool(v) if good else False for v,good in zip(flags,known)])
    missing_expansion=frame.expansion.isin(MISSING).to_numpy()
    compare=valid&known;flag_conflicts=int((compare&(flag_values!=(number>=2))).sum())
    counts=frame.loc[valid].assign(reported_size=number[valid].astype(np.int64)).groupby(['sampleID','reported_size']).size()
    divisibility=sum(int(int(count)%int(size)!=0) for (_,size),count in counts.items())
    by_case=frame.groupby('sampleID').size();too_big=0
    for case,part in frame.loc[valid].assign(reported_size=number[valid]).groupby('sampleID'):
        too_big+=int((part.reported_size>by_case.loc[case]).sum())
    valid_state=valid&known
    accepted=bool((~valid&~explicit_missing).sum()==0 and (~known&~missing_expansion).sum()==0
                  and flag_conflicts==0 and divisibility==0 and too_big==0
                  and valid_state.any() and len(set((number[valid_state]>=2).tolist()))==2)
    profile=dict(source_rows=len(frame),reported_size_positive_integer_rows=int(valid.sum()),explicit_missing_size_rows=int(explicit_missing.sum()),
        unsupported_size_rows=int((~valid&~explicit_missing).sum()),recognized_expansion_flag_rows=int(known.sum()),
        explicit_missing_expansion_rows=int(missing_expansion.sum()),
        expansion_category_count=int(frame.expansion.nunique()),unrecognized_expansion_category_count=int(frame.loc[~known,'expansion'].nunique()),
        size_flag_compared_rows=int(compare.sum()),size_flag_conflict_rows=flag_conflicts,
        material_size_histogram_cells_not_divisible_by_reported_size=divisibility,reported_size_exceeds_material_TCR_rows=too_big,
        operational_reported_expansion_qualification=accepted,
        no_sequence_based_clone_reconstruction_or_identity_inference=True,
        source_size_scope_necessary_checks_not_sequence_proof=True)
    return number,flag_values,valid_state,profile


def match_keys(tcr,md):
    keys=['sampleID','cellID'];assert not md.duplicated(keys).any()
    index=pd.MultiIndex.from_frame(md[keys]);take=index.get_indexer(pd.MultiIndex.from_frame(tcr[keys]))
    return take,take>=0


def qualify(compressed,root,plan):
    root=Path(root);frame=parse_selected(compressed,plan)
    original_rows=len(frame);frame=frame.drop_duplicates(SELECTED).copy();identical_duplicates=original_rows-len(frame)
    ambiguous=frame.duplicated(['sampleID','cellID'],keep=False)
    ambiguous_rows=int(ambiguous.sum());ambiguous_cells=int(frame.loc[ambiguous,['sampleID','cellID']].drop_duplicates().shape[0])
    frame=frame.loc[~ambiguous].copy().reset_index(drop=True)
    number,expanded,known,semantic=expansion_qualification(frame)
    base=read_RNA_source(root,plan['RNA_reference'])
    proof=base['proof']['all243_source_QC_parent_cell_totals']
    assert proof['native_QC500_250_MT20']==proof['all_author_retained'],'ORIGINAL_T_LINEAGE_QC_SET_IDENTITY_REQUIRED'
    mdpath=root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(mdpath)==plan['RNA_reference']['native_metadata_SHA256']
    md=pd.read_csv(mdpath,usecols=['sampleID','cellID','sub_cell_type'],dtype=str,keep_default_na=False)
    assert len(md)==plan['expected_RNA_rows'] and md.sampleID.nunique()==243
    take,matched=match_keys(frame,md);native_type=md.sub_cell_type.to_numpy()[np.maximum(take,0)]
    subtype_mismatch=int((matched&(frame.sub_cell_type.to_numpy()!=native_type)).sum())
    newname_matches=int((matched&(frame.T_new_name.to_numpy()==native_type)).sum())
    cases=sorted(md.sampleID.unique());case_code=pd.Index(cases).get_indexer(frame.sampleID)
    eligible_record=matched&known&(case_code>=0)
    native=root/'NSCLC_full_native_gene_units_20261016_v1';terminal=json.loads((native/'private_terminal.json').read_text())
    assert sha(native/'private_terminal.json')==plan['RNA_reference']['native_terminal_SHA256']
    for name in ['private_full_RNA.npy','private_detected_genes.npy']:
        assert sha(native/name)==terminal['outputs'][name]['sha256']
    rna=np.load(native/'private_full_RNA.npy',allow_pickle=False);genes=np.load(native/'private_detected_genes.npy',allow_pickle=False)
    assert len(rna)==len(genes)==len(md)
    defs=plan['RNA_reference']['author_lineage_names'];groups={name:[name] for name in defs['Treg']+defs['CD8T']}
    groups.update(Treg_POOL=defs['Treg'],CD8_POOL=defs['CD8T'],other_CD4_POOL=defs['CD4T_nonTreg'])
    ncase=len(cases);group_arrays={};source_coverage=[]
    for group,labels in groups.items():
        use=eligible_record&np.isin(native_type,labels);ii=case_code[use]
        nrna=md.sub_cell_type.isin(labels).to_numpy();native_case=pd.Index(cases).get_indexer(md.sampleID)
        parent=np.bincount(native_case[nrna],minlength=ncase)
        called=np.bincount(ii,minlength=ncase);positive=np.bincount(ii,weights=expanded[use],minlength=ncase).astype(int)
        oracle=frame.loc[use,['sampleID']].assign(expanded=expanded[use].astype(int)).groupby('sampleID').expanded.sum().reindex(cases,fill_value=0).to_numpy(int)
        np.testing.assert_array_equal(positive,oracle)
        assert (positive<=called).all() and (called<=parent).all()
        qualities=[]
        for vector in [rna,genes]:
            sums=np.bincount(ii,weights=np.log1p(vector[take[use]]),minlength=ncase)
            qualities.append(sums/np.maximum(called,1))
        group_arrays[group]=dict(RNA=parent,TCR=called,expanded=positive,called_mean_logRNA=qualities[0],called_mean_loggenes=qualities[1])
        source_coverage.append(dict(author_group=group,RNA_rows=int(parent.sum()),TCR_records_matched_with_known_expansion=int(called.sum()),
            expanded_records=int(positive.sum()) if semantic['operational_reported_expansion_qualification'] else None,
            materials_RNA_ge20=int((parent>=20).sum()),materials_TCR_ge20=int((called>=20).sum()),
            materials_any_TCR=int((called>0).sum()),unknown_is_missing_not_nonexpanded=True))
    rows=[r for r in base['rows'] if r['QC']=='native_QC500_250_MT20'];assert len(rows)==243
    for row in rows:
        i=row['case_index'];row['TCR_groups']={group:{key:(float(v[i]) if 'mean_' in key else int(v[i])) for key,v in arrays.items()} for group,arrays in group_arrays.items()}
    cover=[]
    for r in ['Treg_POOL']+defs['Treg']:
        for e in defs['CD8T']:
            chosen=[v for v in rows if v['clinical_eligible'] and min(v['TCR_groups'][k]['TCR'] for k in [r,e,'other_CD4_POOL'])>=20]
            drug=lambda v:re.sub(r'\(C\d+(?:-C\d+)?\)','',str(v)).strip()
            dc=Counter(drug(v['drug']) for v in chosen);chosen=[v for v in chosen if dc[drug(v['drug'])]>=10]
            cover.append(dict(R_group=r,E_group=e,clinical_TCR_and_drug_qualified_materials=len(chosen),
                response1=sum(v['response']==1 for v in chosen),response0=sum(v['response']==0 for v in chosen),
                LUAD=sum(v['cancer_type']=='LUAD' for v in chosen),LUSC=sum(v['cancer_type']=='LUSC' for v in chosen),
                reviewed_independent_patients=0,measurement_semantics_qualified=semantic['operational_reported_expansion_qualification']))
    summary=dict(status='COMPLETE_NONSEQUENCE_TCR_STRUCTURE_AND_MEASUREMENT_PREFLIGHT_NO_EFFECTS',source='GSE243013',
        original_TCR_rows=original_rows,exact_selected_column_duplicate_rows_collapsed=identical_duplicates,
        ambiguous_key_rows_quarantined=ambiguous_rows,ambiguous_cells_quarantined=ambiguous_cells,
        TCR_rows_with_exact_RNA_keys=int(matched.sum()),unmatched_TCR_rows=int((~matched).sum()),
        TCR_subtype_differs_from_RNA_authoritative_rows=subtype_mismatch,TCR_new_name_matches_RNA_authoritative_rows=newname_matches,
        semantics=semantic,all_author_group_capture=source_coverage,all40_pair_measurement_coverage=cover,
        complete_actual_expanded_counts_independent_groupby=True,exact_RNA_and_clinical_SHA_and_original_QC_parent_equivalence=True,
        raw_CDR3_clonotype_and_gene_usage_fields_not_selected=True,
        no_TRG_TRD_annotation_fields_therefore_not_gammadelta_validation=True,
        no_patient_identity_or_same_phase_spatial_or_original_CM_validation=True,
        no_hypothesis_or_clinical_effect_tests=True,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    assert summary['peak_RSS_KiB']<plan['RSS_cap_KiB']
    payload=dict(rows=rows,measurement_semantics_qualified=semantic['operational_reported_expansion_qualification'],RNA_proof=base['proof'])
    return payload,summary


def synthetic():
    header=['sampleID','cellID','sub_cell_type','TRA_v_gene','TRA_j_gene','TRA_c_gene','TRA_cdr3','TRB_v_gene','TRB_j_gene','TRB_c_gene','TRB_cdr3','clonotype','expansion','clonotype_number','T_new_name']
    records=[]
    for i,size in enumerate([1,2,2,3,3,3]):
        record={name:'DO_NOT_SELECT_SEQUENCE' for name in header}
        record.update(sampleID='SYNTHETIC',cellID='cell'+str(i),sub_cell_type='T',expansion='Yes' if size>1 else 'No',clonotype_number=str(size),T_new_name='T')
        records.append(record)
    raw=pd.DataFrame(records,columns=header).to_csv(index=False).encode();compressed=gzip.compress(raw)
    plan=dict(source=dict(bytes=len(compressed)),max_uncompressed_bytes=100000,max_TCR_rows=10,expected_TCR_header=header)
    frame=parse_selected(compressed,plan);assert 'DO_NOT_SELECT_SEQUENCE' not in json.dumps(frame.to_dict())
    sizes,flag,known,p=expansion_qualification(frame);assert p['operational_reported_expansion_qualification'] and known.all() and flag.sum()==5
    broken=frame.iloc[:-1].copy();assert not expansion_qualification(broken)[3]['operational_reported_expansion_qualification']
    missing=frame.copy();missing.loc[0,'clonotype_number']='NA';assert not expansion_qualification(missing)[3]['operational_reported_expansion_qualification']
    extra=frame.iloc[[0]].copy();extra['cellID']='extra';extra['clonotype_number']='NA';extra['expansion']='NA'
    assert expansion_qualification(pd.concat([frame,extra],ignore_index=True))[3]['operational_reported_expansion_qualification']
    md=frame[['sampleID','cellID','sub_cell_type']].copy();idx,valid=match_keys(frame,md);assert valid.all() and idx.tolist()==list(range(6))
    moved=frame.copy();moved.loc[0,'sampleID']='OTHER';assert not match_keys(moved,md)[1][0]


if __name__=='__main__':synthetic();print('PASS_NONSEQUENCE_FIELDS_SIZE_FLAG_UNIT_CHECKS_AND_EXACT_KEYS')
