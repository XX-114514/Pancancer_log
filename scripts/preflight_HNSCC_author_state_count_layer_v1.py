"""All published author states as a separate cohort count layer, no cell relabeling."""
import hashlib,itertools,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

STATES={'author_T':['Tn/Tmem','4-1BB+ Treg','HLA-DR+ Treg','IFN-stim CD4T','LAG3+ Treg','Prolif Treg','GZMK+T','Tmem','Teff','Tstr','Tprolif'],'author_Myeloid':['LA_TAMs','Angio_TAMs','Inflm_TAMs','DCs','Neutrophils','unknown_TAMs','cDC1','HSP+TAM'],'author_B_Plasma':['APB','ASC_B','plasmaB','matureB','MZB','GCB_DLZ','GCB_DZ']}
NESTED={'author_Treg':['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg'],'author_TAM':['LA_TAMs','Angio_TAMs','Inflm_TAMs','unknown_TAMs','HSP+TAM']}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);p=r/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv';clinical=r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv';atlas=r/'complete_external_lineage_opportunities_20261016_v2/private_all_library_latest_metadata.tsv';table=pd.read_csv(p,sep='\t',keep_default_na=False,dtype={'library_id':str,'patient_registry':str});assert len(table)==57 and table.library_id.is_unique
    labels=['malignant']+sum(STATES.values(),[]);assert len(labels)==27 and len(set(labels))==27;counts=table[labels].apply(pd.to_numeric,errors='raise');total=pd.to_numeric(table['Cells analyzed'],errors='raise');assert np.isfinite(counts).all().all() and (counts>=0).all().all() and np.equal(counts,np.floor(counts)).all().all() and np.array_equal(counts.sum(axis=1),total)
    c=pd.read_csv(clinical,sep='\t',keep_default_na=False,dtype={'patient_id':str});assert c.patient_id.is_unique and set(c.patient_id)==set(table.patient_registry);c=c.set_index('patient_id');assert table.neo_type.eq(table.patient_registry.map(c.treatment)).all() and np.array_equal(pd.to_numeric(table['Tumor.regression']),pd.to_numeric(table.patient_registry.map(c.regression)))
    a=pd.read_csv(atlas,sep='\t',keep_default_na=False);a=a[a.arm.eq('archived_immune_atlas')&a.dataset.eq('GSE286827')].set_index('library_id');assert set(a.index)==set(table.library_id)
    out=r/'HNSCC_author_state_count_layer_preflight_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256={'author_join':sha(p),'clinical_definition':sha(clinical),'latest_library_metadata':sha(atlas)},question='Do source-author immune states and their within-parent proportions carry coordinated cohort/clinical information obscured by one broad Treg node? Counts only: not gene expression or exact per-cell barcode annotations.',selection='All27 published cell/state count columns, all57 libraries, all29 author patient keys. All26 immune states in3 author broad parents, plus all4 within-Treg and5 within-TAM fractions; malignant retained as context, not a novel immune mechanism.',frozen_scope='All221 cross-broad-parent immune state pairs; all20 nested Treg/TAM pairs as separate supplementary opportunities. All26 broad immune single-state and9 nested single-state opportunities. No state chosen by clinical effects; actual inference and calibration remain to be specified.',capture='Both author parents>=20 cells at pre/post; zero child counts retained as observed zero, never require each state>=20 for abundance proportions. Gene-expression measurements would require their own cell-support rules.',limits='Author annotation/count/QC version differs from atlas270 and old172. No assigning author per-cell states without exact barcode mapping; no changing authority/CM. Nested denominators are different hypotheses, not independent replication or a rescue family. The Treg atlas TV=0 is structurally induced by one node, not proof of unchanged biological states.',new_author_state_clinical_effects_computed=False,original_identity_or_endpoint_acceptance_changed=False)
    (out/'plan_before_any_author_state_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n');features=[];long=[];parents=[]
    for parent,names in list(STATES.items())+list(NESTED.items()):
        denominator=counts[names].sum(axis=1)
        for name in names:
            feature=parent+'::'+name;features.append(dict(feature_id=feature,parent=parent,author_state=name,denominator='sum_of_explicit_author_states',role='primary_broad' if parent in STATES else 'supplementary_nested',per_cell_mapping_accepted=False))
        for index,row in table.iterrows():
            n=int(denominator.iloc[index]);parents.append(dict(library_id=row.library_id,patient_id=row.patient_registry,timepoint=row.timepoint_registry,parent=parent,n_cells=n,total_author_cells=int(total.iloc[index]),atlas_processed_cells=int(a.loc[row.library_id,'annotated_scope_cells']),denominators_interchangeable=False))
            for name in names:long.append(dict(library_id=row.library_id,patient_id=row.patient_registry,timepoint=row.timepoint_registry,feature_id=parent+'::'+name,parent=parent,author_state=name,state_cells=int(counts.loc[index,name]),parent_cells=n,fraction=float(counts.loc[index,name]/n) if n else np.nan,zero_child_is_retained=bool(counts.loc[index,name]==0 and n>=20),parent_ge20=n>=20,cohort_count_not_per_cell_label=True))
    features=pd.DataFrame(features);long=pd.DataFrame(long);parent_table=pd.DataFrame(parents);assert len(features)==35 and len(long)==57*35 and len(parent_table)==57*5
    paired={};single=[]
    for feature,g in long.groupby('feature_id'):
        eligible=set()
        for patient,gg in g.groupby('patient_id'):
            pre=gg[gg.timepoint.eq('pre')];post=gg[gg.timepoint.eq('post')]
            if len(pre)==len(post)==1 and min(pre.parent_cells.iloc[0],post.parent_cells.iloc[0])>=20:eligible.add(patient)
        paired[feature]=eligible;f=features.set_index('feature_id').loc[feature];single.append(dict(feature_id=feature,parent=f.parent,author_state=f.author_state,role=f.role,reviewed_paired_patients=len(eligible),clinical_test_accepted=False))
    pair_rows=[]
    for role,pp in [('primary_cross_major',list(STATES)),('supplementary_nested_Treg_TAM',list(NESTED))]:
        for A,B in itertools.combinations(pp,2):
            aa=features.loc[features.parent.eq(A),'feature_id'];bb=features.loc[features.parent.eq(B),'feature_id']
            for x,y in itertools.product(aa,bb):pair_rows.append(dict(role=role,parent_A=A,parent_B=B,feature_A=x,feature_B=y,reviewed_paired_patients=len(paired[x]&paired[y]),clinical_test_accepted=False))
    pairs=pd.DataFrame(pair_rows);assert len(pairs)==241 and pairs.role.value_counts().to_dict()=={'primary_cross_major':221,'supplementary_nested_Treg_TAM':20}
    for name,t in [('private_all35_author_immune_state_fraction_definitions.tsv',features),('private_all1995_library_author_state_fractions.tsv',long),('private_all285_library_parent_denominators.tsv',parent_table),('private_all35_state_paired_opportunities.tsv',pd.DataFrame(single)),('private_all241_state_pair_opportunities.tsv',pairs)]:t.to_csv(out/name,sep='\t',index=False,na_rep='NA')
    original=table[['library_id','patient_registry','timepoint_registry','Cells analyzed']+labels].copy();original.to_csv(out/'private_all57_original_author_counts.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL_AUTHOR_STATE_COUNT_SUMS_DENOMINATORS_ZERO_CHILDREN_AND_PAIRED_OPPORTUNITIES_NO_CLINICAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},per_cell_identity_or_new_clinical_statistics=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
