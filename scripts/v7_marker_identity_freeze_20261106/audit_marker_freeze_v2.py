import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';o=r/'stage_c_marker_identity_freeze';a=r/'stage_a_full_axis';policy=json.loads((r/'identity_policy_v2.json').read_text());rules=policy['rules']
v=pd.read_pickle(PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');t=pd.read_pickle(o/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');f=pd.read_pickle(o/'marker_identity_frozen_cells_PRIVATE.pkl.gz');assert t.index.equals(v.index) and t.index.equals(q.index) and t.index.is_unique
assert np.array_equal(t.annotation_frozen.to_numpy(),t.final_annotation.notna().to_numpy());assert t.loc[t.annotation_frozen].index.equals(f.index);assert t.loc[t.annotation_frozen,'QC_keep'].all();assert np.array_equal(t.QC_keep.to_numpy(),v.current_primary_QC_keep.to_numpy())
for col in [c for c in t if c.startswith('frozen_')]:
 x=t[col].astype(object).to_numpy();y=v[col].astype(object).to_numpy();missing=pd.isna(x);assert np.array_equal(missing,pd.isna(y)) and np.array_equal(x[~missing],y[~missing]),col
assert not f.final_annotation.isna().any() and not f.final_annotation.astype(str).str.contains('|'.join(policy['final_label_prohibited_words']),case=False,regex=True).any();assert not f.new_fine_subtype_claim.any()
dec=pd.read_csv(r/'review_decisions/group_name_decisions.tsv',sep='\t',keep_default_na=False).set_index(['audit_group','rule_id']);props=pd.read_csv(r/'stage_b3_refined_naming_proposals/all_group_name_proposals.tsv',sep='\t',keep_default_na=False).set_index(['audit_group','rule_id'])
for (group,j),g in f.groupby(['audit_group','accepted_rule_id'],observed=True):
 j=int(j);d=dec.loc[(group,j)];pr=props.loc[(group,j)];assert d.manual_decision!='do_not_release' and g.final_annotation.astype(str).eq(d.final_label).all() and d.final_label==pr.suggested_label
 p=rules[j];direct=q.loc[g.index,'positive_'+str(j)].to_numpy();aligned=q.loc[g.index,'current_parent'].astype(str).isin(p['aliases']).to_numpy();prior=q.loc[g.index,'already_reviewed_broad'].to_numpy()&aligned
 if d.manual_decision=='accept_manually_reviewed_coherent_context':
  assert aligned.all() and float(pr.parent_matching_fraction)>=.95 and float(pr.fraction_direct_positive)>=.20 and int(pr.n_sources)>=3 and d.context_review_rationale
 else:assert (direct|prior).all()
 if d.manual_decision=='accept_component_route':assert direct.all()
 if d.manual_decision=='accept_prior_reviewed_context':assert prior.all()
 for k,z in enumerate(rules):
  if z['family']==p['family'] or (j==5 and z['parent']=='Plasma') or {z['parent'],p['parent']} in [{'HSPC','Myeloid'},{'HSPC','Granulocyte'}]:continue
  assert not q.loc[g.index,'positive_'+str(k)].any(),(group,j,k)
for j,p in enumerate(rules):
 idx=f.index[f.accepted_rule_id.eq(j)]
 if p['parent']=='NK':assert not q.loc[idx,'positive_0'].any()
 if p['parent'] in ['Bcell','Plasma']:assert not (q.loc[idx,'positive_2']&q.loc[idx,'positive_3']).any()
 if p['parent']=='Plasma':assert not q.loc[idx,'positive_5'].any()
count=pd.read_csv(o/'frozen_identity_counts.tsv',sep='\t');assert count.n_cells.sum()==len(f);status=pd.read_csv(o/'all_cell_release_status_counts.tsv',sep='\t');assert status.n_cells.sum()==len(v)
cross=pd.read_csv(o/'historical_to_new_identity_counts.tsv',sep='\t',keep_default_na=False);assert cross.n_cells.sum()==len(v)
ci=pd.read_csv(o/'frozen_group_counts.tsv',sep='\t');assert ci.n_frozen_cells.sum()==len(f)
pos=0
for chunk in pd.read_csv(o/'marker_identity_frozen_cells_PRIVATE.tsv.gz',sep='\t',usecols=['cell_id','final_annotation','final_parent'],chunksize=100000):
 n=len(chunk);expect=f.iloc[pos:pos+n];assert np.array_equal(chunk.cell_id.to_numpy(),expect.index.to_numpy());assert np.array_equal(chunk.final_annotation.to_numpy(),expect.final_annotation.astype(str).to_numpy());assert np.array_equal(chunk.final_parent.to_numpy(),expect.final_parent.astype(str).to_numpy());pos+=n
assert pos==len(f)
manifest=json.loads((o/'complete.json').read_text())
for z in manifest['files']:
 p=o/z['name'];assert p.stat().st_size==z['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256']
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_all_axis':len(v),'n_frozen':len(f),'n_QC_keep_not_frozen':int((t.QC_keep&~t.annotation_frozen).sum()),'n_labels':f.final_annotation.nunique(),'all_frozen_labels_have_no_progress_placeholders':True,'old_fields_CNV_and_QC_preserved':True,'all_frozen_records_have_reviewed_positive_identity_or_previously_audited_cluster_context':True,'all_observed_cross_family_competitors_remain_outside_freeze':True,'complete_private_TSV_axis_and_values_verified':True,'all_artifact_hashes_verified':True,'all_atlas_publication_acceptance':False};(o/'independent_freeze_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
