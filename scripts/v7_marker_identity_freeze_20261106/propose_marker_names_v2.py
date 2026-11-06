"""Prepare evidence-backed broad marker names and explicit release eligibility."""
import json,re,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';a=r/'stage_a_full_axis';o=r/'stage_b3_refined_naming_proposals';o.mkdir(exist_ok=False)
assert json.loads((r/'launch/full_axis_profile_terminal.json').read_text())['exit_code']==0
policy=json.loads((r/'identity_policy_v2.json').read_text());rules=policy['rules'];q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');summary=pd.read_csv(r/'stage_b2_refined_components/all_group_component_summary.tsv',sep='\t');m=pd.read_csv(a/'all_group_marker_expression.tsv',sep='\t');cur=pd.read_pickle(PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');assert cur.index.equals(q.index)
rows=[];decisions=[]
for row in summary.to_dict('records'):
 group=row['audit_group'];t=q.loc[q.audit_group.eq(group)&q.current_QC_keep];e=m[m.audit_group.eq(group)].set_index('gene');parent=str(row['current_parent_mode']);oldparent=str(row['old_major_mode']);parent_counts=t.current_parent.astype(str).value_counts();n=len(t)
 for j,p in enumerate(rules):
  pos=t['positive_'+str(j)].to_numpy();npositive=int(pos.sum());support=npositive/n
  chosen=[]
  # A representative name must draw from two independent identity components, not generic state/HLA/Ig genes.
  skip={'XBP1','IRF4','SPI1','TYROBP','HLA-DRA','CD74','GZMB','NKG7','PRF1','PTPRC','ACTA2','GATA2','MEIS1'}
  for comp in p['groups']:
   options=[g for g in comp if g not in skip and g in e.index and e.loc[g,'coverage_fraction']>=.8 and e.loc[g,'pct_expressing_measured']>=.05]
   if options:
    best=max(options,key=lambda g:float(e.loc[g,'pct_expressing_measured']));chosen.append(best)
   if len(chosen)==2:break
  # Fall back only to a measured parent-associated gene, never a progress/state token.
  if len(chosen)<2:
   options=[g for comp in p['groups'] for g in comp if g not in chosen and g not in {'HLA-DRA','CD74','GZMB','PTPRC','ACTA2'} and g in e.index and e.loc[g,'coverage_fraction']>=.8 and e.loc[g,'pct_expressing_measured']>=.05]
   options.sort(key=lambda g:float(e.loc[g,'pct_expressing_measured']),reverse=True);chosen+=options[:2-len(chosen)]
  name=p['parent']+'_'+'_'.join(chosen) if len(chosen)==2 else ''
  competitors=[]
  for k,other in enumerate(rules):
   if k==j or other['family']==p['family']:continue
   if {p['parent'],other['parent']} in [{'HSPC','Myeloid'},{'HSPC','Granulocyte'}]:continue
   frac=float(t['positive_'+str(k)].mean())
   if frac>=.05:competitors.append(other['parent']+':'+format(frac,'.4f'))
  aligned=t.current_parent.astype(str).isin(p['aliases']);already=t.already_reviewed_broad&aligned
  basic=n>=20 and row['n_sources']>=3 and bool(name)
  diagnostic=basic and support>=.10
  # This screen nominates clusters for manual review; it never auto-accepts a publication freeze.
  recommended=basic and bool(aligned.mean()>=.8) and support>=.20 and not competitors
  rows.append({'audit_group':group,'rule_id':j,'suggested_parent':p['parent'],'suggested_label':name,'representative_markers':';'.join(chosen),'n_QC_keep':n,'n_direct_positive':npositive,'fraction_direct_positive':support,'n_sources':row['n_sources'],'max_source_fraction':row['max_source_fraction'],'current_parent_mode':parent,'old_major_mode':oldparent,'n_prior_reviewed_broad_matching':int(already.sum()),'n_current_parent_matching':int(aligned.sum()),'parent_matching_fraction':float(aligned.mean()),'competing_parent_components_above_5pct':';'.join(competitors) or 'none','screen_nominates_manual_review':diagnostic,'screen_nominates_coherent_parent_review':recommended,'final_manual_acceptance':False,'naming_scope':'major-plus-marker; does not imply a new subtype, per-cell expression, normality or a functional state'})
  if diagnostic or already.any():decisions.append({'audit_group':group,'rule_id':j,'suggested_label':name,'n_cells':n,'n_direct_positive':npositive,'n_prior_broad_matching':int(already.sum()),'manual_decision':'not_yet_reviewed','manual_reason':'','final_accepted':False})
pd.DataFrame(rows).to_csv(o/'all_group_name_proposals.tsv',sep='\t',index=False);pd.DataFrame(decisions).to_csv(o/'manual_review_entry.tsv',sep='\t',index=False)
receipt={'status':'full_axis_marker_name_proposals_not_freeze_acceptance','UTC':datetime.now(timezone.utc).isoformat(),'n_groups':len(summary),'n_rule_group_records':len(rows),'n_manual_review_entries':len(decisions),'formal_freeze_created':False,'automatic_screen_is_not_manual_acceptance':True};(o/'complete.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
