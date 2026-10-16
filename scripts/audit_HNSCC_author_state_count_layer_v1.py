"""Independently rebuild every author-count denominator and opportunity, no effects."""
import hashlib,itertools,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def compare(actual,expected,keys):
    assert list(actual.columns)==list(expected.columns)
    pd.testing.assert_frame_equal(actual.sort_values(keys).reset_index(drop=True),expected.sort_values(keys).reset_index(drop=True),check_dtype=False,rtol=1e-12,atol=1e-12)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'HNSCC_author_state_count_layer_preflight_20261016_v1';terminal=json.loads((out/'private_terminal.json').read_text())
    for name,h in terminal['files'].items():assert sha(out/name)==h
    p=r/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv';raw=pd.read_csv(p,sep='\t',keep_default_na=False,dtype={'library_id':str,'patient_registry':str});assert len(raw)==57 and raw.library_id.is_unique
    broad={'author_T':['Tn/Tmem','4-1BB+ Treg','HLA-DR+ Treg','IFN-stim CD4T','LAG3+ Treg','Prolif Treg','GZMK+T','Tmem','Teff','Tstr','Tprolif'],'author_Myeloid':['LA_TAMs','Angio_TAMs','Inflm_TAMs','DCs','Neutrophils','unknown_TAMs','cDC1','HSP+TAM'],'author_B_Plasma':['APB','ASC_B','plasmaB','matureB','MZB','GCB_DLZ','GCB_DZ']}
    nested={'author_Treg':['4-1BB+ Treg','HLA-DR+ Treg','LAG3+ Treg','Prolif Treg'],'author_TAM':['LA_TAMs','Angio_TAMs','Inflm_TAMs','unknown_TAMs','HSP+TAM']};all_labels=['malignant']+sum(broad.values(),[]);assert len(all_labels)==len(set(all_labels))==27
    counts=raw[all_labels].apply(pd.to_numeric,errors='raise');total=pd.to_numeric(raw['Cells analyzed'],errors='raise');assert np.isfinite(counts).all().all() and (counts>=0).all().all() and np.equal(counts,counts.astype('int64')).all().all();assert np.array_equal(counts.sum(axis=1),total)
    clinical=pd.read_csv(r/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv',sep='\t',keep_default_na=False,dtype={'patient_id':str}).set_index('patient_id');assert clinical.index.is_unique and set(clinical.index)==set(raw.patient_registry)
    assert np.array_equal(pd.to_numeric(raw['Tumor.regression']),pd.to_numeric(raw.patient_registry.map(clinical.regression)));assert raw.neo_type.eq(raw.patient_registry.map(clinical.treatment)).all()
    atlas=pd.read_csv(r/'complete_external_lineage_opportunities_20261016_v2/private_all_library_latest_metadata.tsv',sep='\t',keep_default_na=False,dtype={'library_id':str});atlas=atlas[atlas.arm.eq('archived_immune_atlas')&atlas.dataset.eq('GSE286827')].set_index('library_id');assert atlas.index.is_unique and set(atlas.index)==set(raw.library_id)
    definitions=[];long=[];parents=[];sets={}
    for parent,labels in list(broad.items())+list(nested.items()):
        for label in labels:definitions.append(dict(feature_id=parent+'::'+label,parent=parent,author_state=label,denominator='sum_of_explicit_author_states',role='primary_broad' if parent in broad else 'supplementary_nested',per_cell_mapping_accepted=False))
        for row in raw.itertuples(index=False,name=None):
            x=dict(zip(raw.columns,row));n=sum(int(x[label]) for label in labels);base=dict(library_id=x['library_id'],patient_id=x['patient_registry'],timepoint=x['timepoint_registry'],parent=parent)
            parents.append(dict(base,n_cells=n,total_author_cells=int(x['Cells analyzed']),atlas_processed_cells=int(atlas.loc[x['library_id'],'annotated_scope_cells']),denominators_interchangeable=False))
            for label in labels:long.append(dict(base,feature_id=parent+'::'+label,author_state=label,state_cells=int(x[label]),parent_cells=n,fraction=int(x[label])/n if n else np.nan,zero_child_is_retained=int(x[label])==0 and n>=20,parent_ge20=n>=20,cohort_count_not_per_cell_label=True))
    defs=pd.DataFrame(definitions);ll=pd.DataFrame(long);pp=pd.DataFrame(parents)
    expected_original=raw[['library_id','patient_registry','timepoint_registry','Cells analyzed']+all_labels]
    compare(pd.read_csv(out/'private_all57_original_author_counts.tsv',sep='\t',dtype={'library_id':str,'patient_registry':str}),expected_original,['library_id'])
    for file,expected,keys in [('private_all35_author_immune_state_fraction_definitions.tsv',defs,['feature_id']),('private_all1995_library_author_state_fractions.tsv',ll,['library_id','feature_id']),('private_all285_library_parent_denominators.tsv',pp,['library_id','parent'])]:
        actual=pd.read_csv(out/file,sep='\t',dtype={'library_id':str,'patient_id':str});compare(actual,expected[actual.columns],keys)
    single=[]
    for row in defs.itertuples():
        g=ll[ll.feature_id.eq(row.feature_id)];eligible=set()
        for patient,rows in g.groupby('patient_id'):
            pre=rows[rows.timepoint.eq('pre')];post=rows[rows.timepoint.eq('post')]
            if len(pre)==len(post)==1 and pre.parent_cells.iloc[0]>=20 and post.parent_cells.iloc[0]>=20:eligible.add(patient)
        sets[row.feature_id]=eligible;single.append(dict(feature_id=row.feature_id,parent=row.parent,author_state=row.author_state,role=row.role,reviewed_paired_patients=len(eligible),clinical_test_accepted=False))
    compare(pd.read_csv(out/'private_all35_state_paired_opportunities.tsv',sep='\t'),pd.DataFrame(single),['feature_id'])
    pairs=[]
    for role,groups in [('primary_cross_major',broad),('supplementary_nested_Treg_TAM',nested)]:
        for a,b in itertools.combinations(groups,2):
            for x,y in itertools.product(groups[a],groups[b]):
                fa=a+'::'+x;fb=b+'::'+y;pairs.append(dict(role=role,parent_A=a,parent_B=b,feature_A=fa,feature_B=fb,reviewed_paired_patients=len(sets[fa].intersection(sets[fb])),clinical_test_accepted=False))
    compare(pd.read_csv(out/'private_all241_state_pair_opportunities.tsv',sep='\t'),pd.DataFrame(pairs),['feature_A','feature_B'])
    proof=dict(status='PASS_ALL57_AUTHOR_COUNTS_35_FEATURES_1995_FRACTIONS_285_DENOMINATORS_35_SINGLE_241_PAIR_OPPORTUNITIES',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),input_author_join_sha256=sha(p),counts_integer_nonnegative_sums_exact=True,zero_children_preserved=True,all_patient_sets_rebuilt=True,per_cell_state_mapping_or_clinical_effects_computed=False)
    dest=out/'private_independent_author_state_count_audit.json';assert not dest.exists();dest.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
