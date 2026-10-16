"""Reviewed clinical keys and complete frozen functional-change rows, private I/O."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd

QC=['n_cells','fraction_of_all_source_library_cells','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo']
NAMES=['intercept','delta_logN','delta_fraction','delta_logRNA','delta_loggenes','delta_pct_mito','delta_pct_ribo','baseline_function','treatment_DplusT','HPV_positive','delta_function']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(root):
    r=Path(root);native=r/'HNSCC_native_multicell_functions_20261016_v1';paired=r/'HNSCC_paired_multicell_function_changes_20261016_v2';clinical=r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv';capture=r/'complete_external_lineage_opportunities_20261016_v2'
    paths=dict(question=r/'HNSCC_longitudinal_function_clinical_question_20261016_v1/PRIVATE_NEXT_QUESTION_AND_SCOPE.json',paired_audit=paired/'private_independent_paired_change_audit.json',coordination_audit=r/'HNSCC_complete_multicell_function_review_20261016_v3/private_independent_complete_review_audit.json',clinical=clinical,author_join=r/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',axis=native/'private_all228_pool_library_axis_and_quality.tsv',scores=native/'private_all_multicell_function_scores_and_QA.npz',masks=native/'private_all16_mask_axis.tsv',references=native/'private_all50_reference_axis.tsv',differences=paired/'private_all_slot_exact_patient_pair_differences.tsv.gz',node_counts=capture/'private_all_original_node_library_counts.tsv.gz',pool_definitions=capture/'private_all_pool_node_definitions.tsv')
    question=json.loads(paths['question'].read_text());assert sha(paths['question'])=='66d35c47ecc401d6eb1ef5481f41a5b94e9ee4e8fca1201a394f19d08a611320'
    for key in ['paired_audit','coordination_audit']:
        assert json.loads(paths[key].read_text())['status'].startswith('PASS_ALL') and sha(paths[key])==question['input_sha256'][key]
    assert sha(clinical)==question['input_sha256']['clinical_definition']
    for directory in [native,paired,capture]:
        terminal=json.loads((directory/'private_terminal.json').read_text())
        for name,p in paths.items():
            if p.parent==directory and p.name in terminal['files']:assert sha(p)==terminal['files'][p.name]
    c=pd.read_csv(clinical,sep='\t',keep_default_na=False,dtype={'patient_id':str});assert len(c)==29 and c.patient_id.is_unique;c=c.set_index('patient_id')
    c['regression']=pd.to_numeric(c.regression,errors='coerce');c['HPV_numeric']=c.HPV.str.strip().str.lower().map({'negative':0.,'positive':1.,'neg':0.,'pos':1.,'-':0.,'+':1.});c['treatment_numeric']=c.treatment.map({'D':0.,'D+T':1.})
    h=pd.read_csv(paths['author_join'],sep='\t',keep_default_na=False,dtype={'library_id':str,'patient_registry':str});assert len(h)==57 and h.library_id.is_unique;assert set(h.patient_registry)==set(c.index);h=h.set_index('library_id')
    axis=pd.read_csv(paths['axis'],sep='\t',dtype={'library_id':str,'patient_block':str});assert len(axis)==228 and axis.pool_index.is_unique;axis=axis.set_index('pool_index')
    for row in axis.itertuples():
        author=h.loc[row.library_id];assert row.patient_block=='GSE286827::author_patient::'+str(author.patient_registry) and row.timepoint==author.timepoint_registry and row.treatment_arm==c.loc[author.patient_registry,'treatment']
    axis['patient_id']=axis.library_id.map(h.patient_registry)
    differences=pd.read_csv(paths['differences'],sep='\t',dtype={'patient_block':str});assert not differences.duplicated(['mask_id','reference','pool','patient_block']).any()
    masks=pd.read_csv(paths['masks'],sep='\t');references=pd.read_csv(paths['references'],sep='\t');assert len(masks)==16 and len(references)==50
    archive=np.load(paths['scores'],allow_pickle=False);scores={key:archive[key] for key in archive.files};groups={key:g.set_index('patient_block') for key,g in differences.groupby(['mask_id','reference','pool'],sort=False)}
    assert len(groups)==3200;mi={v:i for i,v in enumerate(masks.mask_id)};ri={v:i for i,v in enumerate(references.reference)}
    return dict(paths=paths,axis=axis,clinical=c,author_join=h,differences=differences,masks=masks,references=references,scores=scores,groups=groups,mask_index=mi,reference_index=ri,input_sha256={key:sha(p) for key,p in paths.items()})


def clinical_complete(data,patients):
    a=data['axis'];block_to_patient=a[['patient_block','patient_id']].drop_duplicates().set_index('patient_block').patient_id;assert block_to_patient.index.is_unique
    c=data['clinical'].reindex(block_to_patient.reindex(patients).to_numpy());c.index=patients;ok=np.isfinite(c[['regression','HPV_numeric','treatment_numeric']].to_numpy(float)).all(axis=1)
    return [p for p,keep in zip(patients,ok) if keep]


def raw_design(data,mask,reference,pool,patients):
    frame=data['groups'][(mask,reference,pool)].loc[patients];pre=data['axis'].loc[frame.pre_pool_index];post=data['axis'].loc[frame.post_pool_index];assert pre.patient_block.tolist()==post.patient_block.tolist()==patients
    c=data['clinical'].loc[pre.patient_id.to_numpy()];quality=post[QC].to_numpy(float)-pre[QC].to_numpy(float);quality[:,0]=np.log1p(post.n_cells.to_numpy(float))-np.log1p(pre.n_cells.to_numpy(float));i=data['mask_index'][mask];j=data['reference_index'][reference]
    baseline=data['scores']['ULM_score'][i,frame.pre_pool_index.to_numpy(int),j];delta=frame.rounded_difference.to_numpy(float);x=np.column_stack([np.ones(len(frame)),quality,baseline,c.treatment_numeric.to_numpy(float),c.HPV_numeric.to_numpy(float),delta]);y=c.regression.to_numpy(float)
    assert x.shape==(len(patients),11) and len(y)==len(patients)
    strata=(2*c.treatment_numeric+c.HPV_numeric).to_numpy(int)
    return x,y,strata
