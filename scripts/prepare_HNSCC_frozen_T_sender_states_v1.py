"""Freeze all existing T state candidates against the fixed C1QC receiver.
No programme learning, LR selection, clinical endpoint test or type relabeling.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs,score_functions


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    dictionary_path=root/'state_programme_dictionary_20261014_v2/frozen_programme_dictionary.json'
    inventory=json.loads((dictionary_path.parent/'inventory.json').read_text())
    assert hashlib.sha256(dictionary_path.read_bytes()).hexdigest()==inventory['dictionary_sha256']
    dictionary=json.loads(dictionary_path.read_text())
    programmes=[r for r in dictionary['programmes'] if 'T cell' in r['allowed_majors']]
    assert len(programmes)==51
    out=root/'HNSCC_frozen_T_sender_states_20261014_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      dictionary_sha256=inventory['dictionary_sha256'],receiver='previous frozen partial CRC macrophage programme within all eligible C1QC candidate baseline patients',
      target_selected_after_previous_clinical_analysis=True,interpretation='exploratory sender-receiver association, not independent clinical or causal validation',
      candidates='all51 frozen T programmes; primary stable and not foreignRNA/basic-expression flagged; other candidates retained separately',
      no_training_or_signature_changes=True,no_manual_LR_selection=True,no_new_clinical_P_values=True,
      score='fixed primary gene sets, CPU average-tie UCell maxrank1500, measured>=80%,atleast10genes',
      score_sensitivity='separate removal of frozen receiver programme genes from sender signature AND rank axis',
      sample_selection='exact previous baseline CD8-candidate deterministic up-to50cell draw; receiver exact prior subtype draw',
      primary_nuisance='both mean logRNA/logdetected; both candidate within-major fractions; both native log cellcounts,8covariates',
      demographic_sensitivity='same8 plus treatment arm,HPV,sex,age; no auto-drop if design fails',
      min_patients=20,min_residual_df=10,permutations=20000,bootstrap=2000,minimum_valid_bootstrap=1800,
      duplicates='identical observed physical gene sets counted once within each primary family; all original programme IDs retained',
      main_family='all admissible unique physical signatures including unmeasurable/unestimable internally P1',
      sensitivity_families='all51 unique signatures and each sensitivity model separated from primary; not additional confirmations',
      fixed_seed=20261014,original172_and_CM_labels_changed=False,formal_state_labels=False)
    (out/'plan_before_sender_scores_or_associations.json').write_text(json.dumps(plan,indent=2)+'\n')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    chosen=cells[cells.library_id.isin(eligible.library_id)&cells.major.eq('T')&cells.subtype.str.startswith('CD8')&cells.evidence.eq('multigene_subtype_candidate')]
    rng=np.random.default_rng(20261014);positions=[]
    for _,group in chosen.groupby('library_id',sort=True):
        ix=cells.index.get_indexer(group.index);draw=rng.choice(ix,min(50,len(ix)),replace=False)
        if group.timepoint.eq('pre').all():positions.extend(draw)
    positions=np.sort(np.asarray(positions,int));md=cells.iloc[positions]
    raw=ad.io.sparse_dataset(g['layers/counts'])[source_rows[positions],:].tocsr()
    assert np.isfinite(raw.data).all() and (raw.data>=0).all()
    lookup={}
    for i,(name,row) in enumerate(var.iterrows()):
        if not measured[i] or str(name).startswith('AMBIGUOUS::'):continue
        known={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':known.update(str(row.get('source_gene_aliases','')).split('|'))
        for alias in known:
            if alias not in ['','nan','None']:lookup.setdefault(alias,set()).add(i)
    excluded=set(np.flatnonzero(~measured|var.index.astype(str).str.startswith('AMBIGUOUS::')))
    references={r['programme_id']:r['primary_genes'] for r in programmes}
    primary,names,coverage=score_functions(raw,lookup,excluded,references)
    receiver_mapping=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv',sep='\t')
    receiver_columns=set(receiver_mapping.loc[receiver_mapping.native_column.ge(0),'native_column'].astype(int))
    sensitivity,sens_names,sens_coverage=score_functions(raw,lookup,excluded|receiver_columns,references)
    receiver=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    receiver=receiver[receiver.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id').sort_index()
    assert receiver.index.is_unique
    primary_frame=pd.DataFrame(primary,columns=names);primary_frame['patient_id']=md.patient_id.to_numpy()
    primary_means=primary_frame.groupby('patient_id').mean().loc[receiver.index]
    sensitivity_frame=pd.DataFrame(sensitivity,columns=sens_names);sensitivity_frame['patient_id']=md.patient_id.to_numpy()
    sensitivity_means=sensitivity_frame.groupby('patient_id').mean().loc[receiver.index]
    old_means=pd.read_csv(root/'frozen_MCP_two_end_function_review_20261014_v1/private_patient_programme_function_means.tsv',sep='\t')
    sender=old_means[old_means.source.eq('GSE286827')&old_means.endpoint.eq('Tnode9')].set_index('block').loc[receiver.index]
    assert np.array_equal(sender.native_cell_count.to_numpy(),chosen[chosen.timepoint.eq('pre')].groupby('patient_id').size().loc[receiver.index].to_numpy())
    source_major=cells[cells.timepoint.eq('pre')].groupby(['patient_id','major']).size()
    metadata=pd.DataFrame(index=receiver.index)
    metadata['receiver_programme']=receiver.score
    metadata['T_logRNA']=sender.logRNA;metadata['T_loggenes']=sender.loggenes
    metadata['Mac_logRNA']=receiver.mean_logRNA;metadata['Mac_loggenes']=receiver.mean_log_genes
    metadata['T_fraction']=[sender.loc[p,'native_cell_count']/source_major.loc[p,'T'] for p in receiver.index]
    metadata['Mac_fraction']=[receiver.loc[p,'native_subtype_cells']/source_major.loc[p,'Myeloid'] for p in receiver.index]
    metadata['T_logcells']=np.log1p(sender.native_cell_count);metadata['Mac_logcells']=np.log1p(receiver.native_subtype_cells)
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    pre=author[author.time.eq('Pre')].set_index('patient_registry').loc[receiver.index]
    metadata['D_plus_T']=pre.neo_type.eq('D+T').astype(int);metadata['HPV_positive']=pre['HPV/p16'].eq('positive').astype(int)
    metadata['male']=pre.Gender.str.strip().str.lower().map({'m':1,'male':1,'f':0,'female':0});metadata['age']=pd.to_numeric(pre.Age)
    assert metadata.notna().all().all()
    candidate_rows=[]
    for r in programmes:
        columns=sorted({next(iter(lookup[gene])) for gene in r['primary_genes'] if len(lookup.get(gene,set()))==1})
        key=hashlib.sha256(json.dumps(columns).encode()).hexdigest()
        row={k:r[k] for k in ['programme_id','method','stability_main_candidate','foreign_RNA_sensitive','basic_expression_review_flag']}
        row['primary_admissible']=bool(r['stability_main_candidate'] and not r['foreign_RNA_sensitive'] and not r['basic_expression_review_flag'])
        row['observed_signature_group']=key;row['measurable']=r['programme_id'] in names
        row['orthogonal_measurable']=r['programme_id'] in sens_names
        candidate_rows.append(row)
    candidates=pd.DataFrame(candidate_rows)
    candidates.to_csv(out/'private_frozen_candidate_contract.tsv',sep='\t',index=False)
    metadata.to_csv(out/'private_patient_receiver_covariates.tsv',sep='\t')
    primary_means.to_csv(out/'private_patient_T_state_means.tsv',sep='\t')
    sensitivity_means.to_csv(out/'private_patient_T_orthogonal_state_means.tsv',sep='\t')
    pd.DataFrame(coverage).to_csv(out/'private_primary_signature_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(sens_coverage).to_csv(out/'private_orthogonal_signature_coverage.tsv',sep='\t',index=False)
    terminal=dict(status='ALL51_FROZEN_T_SENDER_STATES_AND_COVARIATES_PREPARED_NO_ASSOCIATION_TESTS',
      finished_UTC=datetime.now(timezone.utc).isoformat(),patient_blocks=len(metadata),candidate_ids=len(candidates),
      primary_admissible_ids=int(candidates.primary_admissible.sum()),primary_unique_groups=candidates[candidates.primary_admissible].observed_signature_group.nunique(),
      all_unique_groups=candidates.observed_signature_group.nunique(),measurable_candidates=len(names),orthogonal_measurable_candidates=len(sens_names),
      sender_selected_cells=len(md),scalar_UCell_parity=True,new_clinical_P_values=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
