"""Create a versioned positive-evidence freeze; never force uncertain parents into labels."""
import json,re,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';a=r/'stage_a_full_axis';b=r/'stage_b3_refined_naming_proposals';o=r/'stage_c_marker_identity_freeze'
assert json.loads((a/'independent_profile_audit.json').read_text())['status']=='PASS'
assert json.loads((r/'stage_b2_refined_components/independent_refinement_audit.json').read_text())['status']=='PASS'
policy=json.loads((r/'identity_policy_v2.json').read_text());rules=policy['rules'];q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');cur=pd.read_pickle(PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');assert q.index.equals(cur.index)
decision=pd.read_csv(r/'review_decisions/group_name_decisions.tsv',sep='\t',keep_default_na=False);assert not decision.duplicated(['audit_group','rule_id']).any();assert decision.manual_decision.isin(['accept_component_route','accept_prior_reviewed_context','accept_component_and_prior_context','accept_manually_reviewed_coherent_context','do_not_release']).all()
props=pd.read_csv(b/'all_group_name_proposals.tsv',sep='\t',keep_default_na=False).set_index(['audit_group','rule_id']);raw=np.load(a/'marker_counts_PRIVATE.npy',mmap_mode='r');genes=np.load(a/'gene_symbols.npy').tolist();gx={g:i for i,g in enumerate(genes)}
t=cur[['frozen_parent','frozen_label','frozen_confidence','frozen_strict_major_eligible','frozen_strict_subtype_eligible']+[c for c in cur.columns if c.startswith('frozen_') and ('malignan' in c or 'evidence_version' in c)]].copy();assert t.columns.is_unique
t['QC_keep']=cur.current_primary_QC_keep.to_numpy();t['QC_action']=cur.current_QC_action.astype(object);t['previous_working_parent']=cur.current_parent.astype(object);t['audit_group']=q.audit_group.astype(object)
t['annotation_frozen']=False;t['final_parent']=pd.Series(pd.NA,index=t.index,dtype=object);t['final_annotation']=pd.Series(pd.NA,index=t.index,dtype=object);t['representative_markers']=pd.Series(pd.NA,index=t.index,dtype=object);t['freeze_basis']='not_accepted_in_this_identity_release';t['freeze_status']=np.where(t.QC_keep,'identity_evidence_incomplete','outside_primary_QC');t['accepted_rule_id']=pd.Series(pd.NA,index=t.index,dtype='Int64');t['fine_subtype']=pd.Series(pd.NA,index=t.index,dtype=object);t['cluster_state_context']=cur.cluster_RNA_context.astype(object);t['malignancy_context_status']=cur.malignancy_axis_review_status.astype(object);t['new_fine_subtype_claim']=False
eligible_by_rule=np.zeros((len(q),len(rules)),bool);prior_by_rule=np.zeros_like(eligible_by_rule);context_by_rule=np.zeros_like(eligible_by_rule)
group_cat=pd.Categorical(q.audit_group);group_codes=group_cat.codes;group_map={str(g):i for i,g in enumerate(group_cat.categories)};parent_values=q.current_parent.astype(str).to_numpy();aligned_by_rule=[np.isin(parent_values,p['aliases']) for p in rules]
for d in decision.to_dict('records'):
 if d['manual_decision']=='do_not_release':continue
 j=int(d['rule_id']);p=rules[j];key=(d['audit_group'],j);assert key in props.index;pr=props.loc[key];assert d['final_label']==pr.suggested_label and d['final_label'] and d['manual_reason'];assert not any(word in d['final_label'].lower() for word in policy['final_label_prohibited_words'])
 group=(group_codes==group_map[d['audit_group']])&q.current_QC_keep.to_numpy();direct=q['positive_'+str(j)].to_numpy();oldmatch=aligned_by_rule[j];prior=q.already_reviewed_broad.to_numpy()&oldmatch
 if d['manual_decision'] in ['accept_component_route','accept_component_and_prior_context']:eligible_by_rule[:,j]|=group&direct
 if d['manual_decision'] in ['accept_prior_reviewed_context','accept_component_and_prior_context']:eligible_by_rule[:,j]|=group&prior;prior_by_rule[:,j]|=group&prior
 if d['manual_decision']=='accept_manually_reviewed_coherent_context':
  assert float(pr.parent_matching_fraction)>=.95 and float(pr.fraction_direct_positive)>=.20 and int(pr.n_sources)>=3 and d['context_review_rationale']
  eligible_by_rule[:,j]|=group&oldmatch;context_by_rule[:,j]|=group&oldmatch
# Cross-family positive programmes stay outside this freeze. No fresh doublet/QC deletion is made here.
positive=q[['positive_'+str(j) for j in range(len(rules))]].to_numpy(bool);blocked=np.zeros_like(eligible_by_rule)
for j,p in enumerate(rules):
 competitors=[]
 for k,z in enumerate(rules):
  if p['family']==z['family'] or (j==5 and z['parent']=='Plasma'):continue
  if {p['parent'],z['parent']} in [{'HSPC','Myeloid'},{'HSPC','Granulocyte'}]:continue
  competitors.append(k)
 blocked[:,j]=eligible_by_rule[:,j]&positive[:,competitors].any(axis=1)
 eligible_by_rule[:,j]&=~blocked[:,j]
# Do not call an NK record when a full positive T apparatus is observed, or force a B/ASC boundary.
t_ids=[j for j,p in enumerate(rules) if p['parent']=='Tcell'];nk_ids=[j for j,p in enumerate(rules) if p['parent']=='NK'];b_ids=[j for j,p in enumerate(rules) if p['parent']=='Bcell'];asc_ids=[j for j,p in enumerate(rules) if p['parent']=='Plasma']
for j in nk_ids:eligible_by_rule[positive[:,t_ids].any(axis=1),j]=False
basc_shared=positive[:,b_ids].any(axis=1)&positive[:,asc_ids].any(axis=1)
for j in b_ids+asc_ids:eligible_by_rule[basc_shared,j]=False
for j in asc_ids:eligible_by_rule[positive[:,5],j]=False
# Broad myeloid ancestry is compatible with a more specific positive DC or granulocyte apparatus.
specific_myeloid=[j for j,p in enumerate(rules) if p['parent'] in ['Dendritic','Granulocyte']]
for j,p in enumerate(rules):
 if p['parent']=='Myeloid':eligible_by_rule[eligible_by_rule[:,specific_myeloid].any(axis=1),j]=False
# T/NK, B/ASC and DC branches can share programmes. Do not select the strongest score to force a label.
parents=np.asarray([x['parent'] for x in rules]);nparent=np.zeros(len(q),np.int16)
for par in np.unique(parents):nparent+=eligible_by_rule[:,parents==par].any(axis=1)
ambiguous=nparent>1;eligible_by_rule[ambiguous]=False
for j,p in enumerate(rules):
 # Multiple DC branch rules remain broad DC but their marker label needs one reviewed route.
 same=parents==p['parent'];multiple_same=eligible_by_rule[:,same].sum(axis=1)>1;eligible_by_rule[multiple_same,:]=False
for d in decision.to_dict('records'):
 if d['manual_decision']=='do_not_release':continue
 j=int(d['rule_id']);p=rules[j];take=(group_codes==group_map[d['audit_group']])&eligible_by_rule[:,j];assert not t.annotation_frozen.iloc[np.flatnonzero(take)].any()
 t.loc[take,'final_parent']=p['parent'];t.loc[take,'final_annotation']=d['final_label'];t.loc[take,'representative_markers']=props.loc[(d['audit_group'],j),'representative_markers'];t.loc[take,'accepted_rule_id']=j;t.loc[take,'annotation_frozen']=True;t.loc[take,'freeze_status']='broad_identity_marker_frozen';t.loc[take,'freeze_basis']=np.select([context_by_rule[take,j],prior_by_rule[take,j]],['manually_reviewed_coherent_multi_marker_context','previous_independent_cluster_review_plus_current_marker_review'],default='positive_independent_identity_components_plus_manual_group_review')
t.loc[t.QC_keep&~t.annotation_frozen&blocked.any(axis=1),'freeze_status']='positive_competing_lineage_programmes';t.loc[t.QC_keep&~t.annotation_frozen&ambiguous,'freeze_status']='related_parent_identity_not_separated'
t['raw_IFN_programme_positive']=((raw[:,[gx[g] for g in ['ISG15','IFIT1','IFIT3','MX1']]]>=2).sum(axis=1)>=2);t['raw_cycling_programme_positive']=((raw[:,[gx[g] for g in ['MKI67','TOP2A','UBE2C','CENPF']]]>=2).sum(axis=1)>=2);t['raw_heat_shock_programme_positive']=((raw[:,[gx[g] for g in ['DNAJB1','HSPA1A','HSPA1B']]]>=2).sum(axis=1)>=2)
t['state_semantics']='positive RNA programme evidence; not a functional state confirmation';t['release_version']='V7_marker_identity_20261106_v2'
changed=np.zeros(len(t),bool);selected=t.annotation_frozen.to_numpy();changed[selected]=t.loc[selected,'final_parent'].astype(str).to_numpy()!=t.loc[selected,'previous_working_parent'].astype(str).to_numpy();t['new_parent_vs_previous_working']=changed
t.loc[t.new_parent_vs_previous_working,'malignancy_context_status']='identity_parent_or_grain_changed_recheck_frozen_CNV_context_only'
t['annotation_resolution']='major_identity_plus_markers';t['identity_freeze_accepted']=t.annotation_frozen
assert not t.final_annotation.dropna().str.contains('|'.join(policy['final_label_prohibited_words']),case=False,regex=True).any();assert t.loc[t.annotation_frozen,'QC_keep'].all() and t.final_annotation.notna().equals(t.annotation_frozen)
o.mkdir(exist_ok=False)
for col in t.select_dtypes(include=object):t[col]=t[col].astype('category')
p=o/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz';t.to_pickle(p);back=pd.read_pickle(p);pd.testing.assert_frame_equal(t,back)
frozen=t.loc[t.annotation_frozen].copy();p2=o/'marker_identity_frozen_cells_PRIVATE.pkl.gz';frozen.to_pickle(p2);pd.testing.assert_frame_equal(frozen,pd.read_pickle(p2));frozen.to_csv(o/'marker_identity_frozen_cells_PRIVATE.tsv.gz',sep='\t',compression={'method':'gzip','compresslevel':1},na_rep='NA')
cols=['final_parent','final_annotation','representative_markers','freeze_basis'];frozen.groupby(cols,observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'frozen_identity_counts.tsv',sep='\t',index=False)
t.groupby(['freeze_status','QC_action'],observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'all_cell_release_status_counts.tsv',sep='\t',index=False)
t.groupby(['frozen_parent','frozen_label','final_parent','final_annotation','freeze_status'],observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'historical_to_new_identity_counts.tsv',sep='\t',index=False,na_rep='')
frozen.groupby(['audit_group','accepted_rule_id','final_parent','final_annotation'],observed=True,dropna=False).size().rename('n_frozen_cells').reset_index().to_csv(o/'frozen_group_counts.tsv',sep='\t',index=False)
receipt={'status':'versioned_marker_identity_subset_written_not_complete_atlas_publication_acceptance','UTC':datetime.now(timezone.utc).isoformat(),'n_original_axis':len(t),'n_QC_keep':int(t.QC_keep.sum()),'n_marker_identity_frozen':len(frozen),'n_QC_keep_not_yet_in_freeze':int((t.QC_keep&~t.annotation_frozen).sum()),'n_final_label_values':frozen.final_annotation.nunique(),'n_final_parent_values':frozen.final_parent.nunique(),'n_formally_new_fine_subtype_claims':0,'forbidden_progress_tokens_in_final_labels':False,'frozen_old_V7_changed':False,'current_QC_changed':False,'new_malignancy_calls':False,'all_atlas_publication_accepted':False,'files':[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir() if p.is_file()]};(o/'complete.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))
