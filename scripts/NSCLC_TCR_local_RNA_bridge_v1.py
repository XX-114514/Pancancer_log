"""Local-only nonsequence reported expansion and original RNA/clinical bridge.

All identifiers/records remain local RAM. No receptor sequences or clonotype
IDs selected; no statistical effects until a separate frozen contract.
"""
from collections import Counter
import io
import json
import re
import numpy as np
import pandas as pd
import qualify_NSCLC_reported_TCR_expansion_RAM_v1 as original

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
FIELDS=['sampleID','age','gender','smoking_history','cancer_type','pre_treatment_staging','anti-PD1_therapy','chemotherapy','targeted_therapy','cycles','pathological_response','pathological_response_rate','radiological_response']
RNA_COLUMNS=list(dict.fromkeys(['sampleID','cellID','sub_cell_type','total_counts','n_genes_by_counts','total_counts_mt']+FIELDS))


def semantics(frame):
    number,flag,known,previous=original.expansion_qualification(frame)
    valid=known&np.isfinite(number)&(number>=1)&(number==np.floor(number))
    has_true=bool((valid&flag).any());has_false=bool((valid&~flag).any())
    lo=int(number[valid&flag].min()) if has_true else None
    hi=int(number[valid&~flag].max()) if has_false else None
    monotone=has_true and has_false and hi<lo
    accepted=bool(previous['unsupported_size_rows']==0 and previous['unrecognized_expansion_category_count']==0
        and previous['material_size_histogram_cells_not_divisible_by_reported_size']==0
        and previous['reported_size_exceeds_material_TCR_rows']==0 and monotone)
    bins=[('size1',number==1),('size2',number==2),('size3to5',(number>=3)&(number<=5)),
          ('size6to20',(number>=6)&(number<=20)),('size21to100',(number>=21)&(number<=100)),('sizeAbove100',number>100)]
    table=[dict(reported_size_bin=name,source_records=int((valid&mask).sum()),author_flag_true=int((valid&mask&flag).sum()),
                author_flag_false=int((valid&mask&~flag).sum())) for name,mask in bins]
    profile=dict(status='QUALIFIED_SOURCE_REPORTED_EXPANSION_NOT_SEQUENCE_RECONSTRUCTION' if accepted else 'REPORTED_EXPANSION_SEMANTICS_PENDING',
        previous_size_ge2_diagnostic=previous,source_flag_unchanged=True,observed_flag_true_min_reported_size=lo,
        observed_flag_false_max_reported_size=hi,ordered_boundary_interval=[hi+1,lo] if monotone else None,
        unique_integer_boundary=lo if monotone and lo==hi+1 else None,
        source_size_flag_bin_table=table,
        diagnostic_disagreements_ge2=int((valid&(flag!=(number>=2))).sum()),
        diagnostic_disagreements_ge3=int((valid&(flag!=(number>=3))).sum()),
        operational_reported_expansion_qualified=accepted,no_outcome_or_effect_used_to_qualify_source_field=True,
        necessary_within_material_count_consistency_not_independent_clone_assignment_proof=True)
    return number,flag,valid,profile


def clinical(md):
    raw=md[FIELDS].drop_duplicates().set_index('sampleID',drop=False)
    assert raw.index.is_unique and len(raw)==243
    raw=raw.loc[sorted(raw.index)].copy();c=pd.DataFrame(index=raw.index)
    c['response']=raw.pathological_response.map({'pCR':1.,'MPR':1.,'non-MPR':0.})
    c['age_numeric']=pd.to_numeric(raw.age,errors='coerce');c['male']=raw.gender.map({'M':1.,'F':0.})
    c['smoker']=raw.smoking_history.map({'Y':1.,'N':0.});c['LUSC']=raw.cancer_type.map({'LUSC':1.,'LUAD':0.})
    c['stage_ordinal']=raw.pre_treatment_staging.map(lambda s:next((v for prefix,v in [('IV',4.),('III',3.),('II',2.),('I',1.)] if s.startswith(prefix)),np.nan))
    c['cycles_numeric']=pd.to_numeric(raw.cycles,errors='coerce')
    c['chemo_binary']=raw.chemotherapy.map(lambda s:0. if s=='No' else 1. if s=='Yes' or any(v in s for v in ['Carboplatin','Cisplatin','Nedaplatin','Lobaplatin','Abraxane','Pemetrexed','Paclitaxel','Gemcitabine','Docetaxel']) else np.nan)
    c['targeted_binary']=raw.targeted_therapy.map(lambda s:0. if s=='No' else 1. if any(v in s for v in ['Apatinib','Bevacizumab','Avastin','Endostar','Anlotinib']) else np.nan)
    drugs='(?:Tislelizumab|Pembrolizumab|Camrelizumab|Sintilimab|Nivolumab|Toripalimab)';token=drugs+r'(?:\(C\d+(?:-C\d+)?\))?'
    explicit=re.compile('^'+token+'(?: \\+ '+token+')*$')
    roles=raw['anti-PD1_therapy'].map(lambda s:'EXPLICIT_ICI_AUTHOR_CASE' if s=='Yes' or explicit.fullmatch(s) else 'NO_ICI' if s=='No' else 'AMBIGUOUS_OR_UNKNOWN_ICI')
    c['eligible']=roles.eq('EXPLICIT_ICI_AUTHOR_CASE')&np.isfinite(c[COMMON+['response']].to_numpy(float)).all(axis=1)
    return raw,c,roles


def measure(tcr,md,plan):
    before=len(tcr);tcr=tcr.drop_duplicates(original.SELECTED).copy()
    duplicate=tcr.duplicated(['sampleID','cellID'],keep=False);ambiguous=int(duplicate.sum())
    tcr=tcr.loc[~duplicate].reset_index(drop=True)
    number,flag,valid_state,sem=semantics(tcr)
    assert len(md)==plan['expected_RNA_rows'] and not md.duplicated(['sampleID','cellID']).any()
    take,matched=original.match_keys(tcr,md);native=md.sub_cell_type.to_numpy()[np.maximum(take,0)]
    raw,c,roles=clinical(md);cases=raw.index.tolist()
    assert roles.value_counts().to_dict()==plan['expected_clinical_role_counts']
    all_case=pd.Index(cases).get_indexer(md.sampleID);tc_case=pd.Index(cases).get_indexer(tcr.sampleID)
    total=pd.to_numeric(md.total_counts,errors='coerce').to_numpy(float)
    genes=pd.to_numeric(md.n_genes_by_counts,errors='coerce').to_numpy(float)
    mito=pd.to_numeric(md.total_counts_mt,errors='coerce').to_numpy(float)
    finite=np.isfinite(total)&np.isfinite(genes)&np.isfinite(mito)
    qc=finite&(total>=500)&(genes>=250)&(mito<=.2*total)
    defs=plan['author_lineage_names'];typed=md.sub_cell_type.isin(sum(defs.values(),[])).to_numpy()
    assert (qc[typed]).all(),'REPORTED_QC_SET_NOT_EQUAL_TO_PREVIOUS_NATIVE_AUTHOR_LINEAGE_SET'
    groups={name:[name] for name in defs['Treg']+defs['CD8T']}
    groups.update(Treg_POOL=defs['Treg'],CD8_POOL=defs['CD8T'],other_CD4_POOL=defs['CD4T_nonTreg'])
    eligible_record=matched&valid_state&(tc_case>=0)&qc[np.maximum(take,0)]
    arrays={};coverage=[]
    for name,labels in groups.items():
        use=eligible_record&np.isin(native,labels);ii=tc_case[use]
        parent=qc&md.sub_cell_type.isin(labels).to_numpy()
        nRNA=np.bincount(all_case[parent],minlength=243);called=np.bincount(ii,minlength=243)
        positive=np.bincount(ii,weights=flag[use],minlength=243).astype(int)
        ge2=np.bincount(ii,weights=number[use]>=2,minlength=243).astype(int)
        oracle=tcr.loc[use,['sampleID']].assign(flag=flag[use].astype(int)).groupby('sampleID').flag.sum().reindex(cases,fill_value=0).to_numpy(int)
        np.testing.assert_array_equal(positive,oracle);assert (positive<=called).all() and (called<=nRNA).all()
        qr=[]
        for value in [total,genes]:
            qr.append(np.bincount(ii,weights=np.log1p(value[take[use]]),minlength=243)/np.maximum(1,called))
        arrays[name]=dict(RNA=nRNA,TCR=called,expanded=positive,expanded_ge2=ge2,called_mean_logRNA=qr[0],called_mean_loggenes=qr[1])
        coverage.append(dict(author_group=name,RNA_cells=int(nRNA.sum()),known_TCR_expansion_calls=int(called.sum()),
            source_reported_expanded_cells=int(positive.sum()) if sem['operational_reported_expansion_qualified'] else None,
            materials_TCR_ge20=int((called>=20).sum()),materials_any_TCR=int((called>0).sum()),
            source_TCR_label_differs_from_RNA_rows=int((use&(tcr.sub_cell_type.to_numpy()!=native)).sum())))
    rows=[]
    for i in range(243):
        rows.append(dict(case_index=i,clinical=[float(v) if np.isfinite(v) else None for v in c[COMMON].iloc[i]],
            clinical_eligible=bool(c.eligible.iloc[i]),response=float(c.response.iloc[i]) if np.isfinite(c.response.iloc[i]) else None,
            drug=raw['anti-PD1_therapy'].iloc[i],cancer_type=raw.cancer_type.iloc[i],
            TCR_groups={name:{key:float(v[i]) if 'mean_' in key else int(v[i]) for key,v in a.items()} for name,a in arrays.items()}))
    pairs=[];agent=lambda s:re.sub(r'\(C\d+(?:-C\d+)?\)','',str(s)).strip()
    for r in ['Treg_POOL']+defs['Treg']:
        for e in defs['CD8T']:
            selected=[v for v in rows if v['clinical_eligible'] and min(v['TCR_groups'][name]['TCR'] for name in [r,e,'other_CD4_POOL'])>=20]
            dc=Counter(agent(v['drug']) for v in selected);selected=[v for v in selected if dc[agent(v['drug'])]>=10]
            pairs.append(dict(R_group=r,E_group=e,qualified_author_materials=len(selected),response1=sum(v['response']==1 for v in selected),response0=sum(v['response']==0 for v in selected),
                LUAD=sum(v['cancer_type']=='LUAD' for v in selected),LUSC=sum(v['cancer_type']=='LUSC' for v in selected),reviewed_independent_patients=0))
    summary=dict(status='COMPLETE_LOCAL_NONSEQUENCE_TCR_EXACT_RNA_AND_FROZEN_CLINICAL_CODING_MEASUREMENT_NO_EFFECTS',source='GSE243013',
        original_TCR_rows=before,exact_selected_duplicates_collapsed=before-ambiguous-len(tcr),ambiguous_rows_quarantined=ambiguous,
        exact_matched_RNA_cells=int(matched.sum()),unmatched_TCR_cells=int((~matched).sum()),
        TCR_original_label_differs_from_RNA_rows=int((matched&(tcr.sub_cell_type.to_numpy()!=native)).sum()),
        TCR_new_name_equals_RNA_rows=int((matched&(tcr.T_new_name.to_numpy()==native)).sum()),
        semantic_profile=sem,all_author_group_capture=coverage,all40_clinical_capture_positions=pairs,
        all_original_RNA_author_T_lineage_reported_QC_sets_equal_retained=True,
        strict_original_clinical_roles_reproduced=True,clinical_coding_source='Same original published metadata SHA and frozen original explicit coding, no private-case file transfer',
        all_expanded_source_counts_independent_groupby_checked=True,
        authoritative_RNA_cell_type_via_exact_keys_no_string_alias_or_favorable_source_label_selection=True,
        no_sequence_or_raw_clonotype_fields_selected=True,no_remote_process_or_individual_network_export=True,
        clinical_or_coordination_effects_not_computed=True,not_independent_patient_baseline_or_original_CM_validation=True)
    return dict(rows=rows,semantics_qualified=sem['operational_reported_expansion_qualified']),summary


def synthetic():
    frame=pd.DataFrame(dict(sampleID=['CASE']*6,cellID=list('abcdef'),sub_cell_type=['T']*6,
        expansion=['No','No','No','Yes','Yes','Yes'],clonotype_number=['1','2','2','3','3','3'],T_new_name=['T']*6))
    n,f,k,p=semantics(frame)
    assert p['operational_reported_expansion_qualified'] and p['unique_integer_boundary']==3
    assert p['diagnostic_disagreements_ge2']==2 and p['diagnostic_disagreements_ge3']==0
    assert f.sum()==3 and (n>=2).sum()==5
    broken=frame.iloc[:-1].copy();assert not semantics(broken)[3]['operational_reported_expansion_qualified']
    crossing=frame.copy();crossing.loc[0,'expansion']='Yes';assert not semantics(crossing)[3]['operational_reported_expansion_qualified']
    original.synthetic()


if __name__=='__main__':synthetic();print('PASS_SOURCE_FLAG_NOT_GE2_AND_EXACT_FIELDS')
