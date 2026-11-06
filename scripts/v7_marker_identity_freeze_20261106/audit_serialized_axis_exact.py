"""Exact, vectorized audit of a completed serialized full-axis freeze sidecar."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';p=r/'stage_c_marker_identity_freeze/annotation_all_cells_with_release_status_PRIVATE.pkl.gz';o=r/'stage_c2_exact_audit_resume';o.mkdir(exist_ok=False);start=datetime.now(timezone.utc).isoformat();before=p.stat().st_size
t=pd.read_pickle(p);assert before==p.stat().st_size
q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');v=pd.read_pickle(PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');assert t.index.equals(q.index) and t.index.equals(v.index) and t.index.is_unique
for col in [x for x in t if x.startswith('frozen_')]:
 x=t[col].astype(object).to_numpy();y=v[col].astype(object).to_numpy();na=pd.isna(x);assert np.array_equal(na,pd.isna(y)) and np.array_equal(x[~na],y[~na]),col
assert np.array_equal(t.QC_keep.to_numpy(),q.current_QC_keep.to_numpy()) and np.array_equal(t.QC_action.astype(str).to_numpy(),v.current_QC_action.astype(str).to_numpy());assert np.array_equal(t.audit_group.astype(str).to_numpy(),q.audit_group.astype(str).to_numpy())
policy=json.loads((r/'identity_policy_v2.json').read_text());rules=policy['rules'];dec=pd.read_csv(r/'review_decisions/group_name_decisions.tsv',sep='\t',keep_default_na=False);N=len(q);positive=q[['positive_'+str(j) for j in range(22)]].to_numpy(bool);allowed=np.zeros((N,22),bool);gc=pd.Categorical(q.audit_group);gcodes=gc.codes;gmap={str(g):i for i,g in enumerate(gc.categories)};pc=pd.Categorical(q.current_parent);pcodes=pc.codes;keep=q.current_QC_keep.to_numpy();prior=q.already_reviewed_broad.to_numpy();names=np.full((len(gmap),22),'',object)
for d in dec.to_dict('records'):
 if d['manual_decision']=='do_not_release':continue
 j=int(d['rule_id']);p0=rules[j];group=(gcodes==gmap[d['audit_group']])&keep;align=pd.Index(pc.categories).isin(p0['aliases'])[pcodes];basis=d['manual_decision'];choice=np.zeros(N,bool)
 if basis in ['accept_component_route','accept_component_and_prior_context']:choice|=positive[:,j]
 if basis in ['accept_prior_reviewed_context','accept_component_and_prior_context']:choice|=prior&align
 if basis=='accept_manually_reviewed_coherent_context':choice|=align
 allowed[:,j]|=group&choice;names[gmap[d['audit_group']],j]=d['final_label']
for j,p0 in enumerate(rules):
 competitor=[k for k,z in enumerate(rules) if p0['family']!=z['family'] and not (j==5 and z['parent']=='Plasma') and {p0['parent'],z['parent']} not in [{'HSPC','Myeloid'},{'HSPC','Granulocyte'}]];allowed[:,j]&=~positive[:,competitor].any(axis=1)
allowed[positive[:,0],1]=False;both=positive[:,2]&positive[:,3];allowed[both,2]=False;allowed[both,3]=False;allowed[positive[:,5],3]=False;allowed[allowed[:,[5,6,7,10]].any(axis=1),4]=False
parents=np.asarray([x['parent'] for x in rules]);nparent=np.zeros(N,np.int8)
for par in np.unique(parents):nparent+=allowed[:,parents==par].any(axis=1)
allowed[nparent>1]=False
for par in np.unique(parents):allowed[allowed[:,parents==par].sum(axis=1)>1]=False
expected=allowed.any(axis=1);assert np.array_equal(expected,t.annotation_frozen.to_numpy());assert np.array_equal(expected,t.final_annotation.notna().to_numpy());rule=allowed.argmax(axis=1);assert np.array_equal(t.loc[expected,'accepted_rule_id'].astype(int).to_numpy(),rule[expected]);assert np.array_equal(t.loc[expected,'final_parent'].astype(str).to_numpy(),parents[rule[expected]]);assert np.array_equal(t.loc[expected,'final_annotation'].astype(str).to_numpy(),names[gcodes[expected],rule[expected]])
assert not t.loc[expected,'final_annotation'].astype(str).str.contains('|'.join(policy['final_label_prohibited_words']),case=False,regex=True).any();assert not t.new_fine_subtype_claim.any()
h=hashlib.sha256(p.read_bytes()).hexdigest();assert before==p.stat().st_size
out={'status':'PASS_full_axis_exact_value_and_complete_membership_audit','UTC_start':start,'UTC_end':datetime.now(timezone.utc).isoformat(),'n_all_axis':N,'n_QC_keep':int(keep.sum()),'n_frozen':int(expected.sum()),'n_QC_keep_outside_identity_freeze':int((keep&~expected).sum()),'n_final_labels':t.loc[expected,'final_annotation'].nunique(),'n_final_parents':t.loc[expected,'final_parent'].nunique(),'all_original_fields_and_CNV_exact':True,'all_names_rule_ids_parents_and_inclusion_mask_exactly_recomputed_from_all_manual_decisions':True,'source_full_axis_file_bytes':before,'source_full_axis_sha256':h,'source_file_unchanged_during_audit':True,'no_missing_identity_assigned_a_biological_name':True};(o/'independent_serialized_axis_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
