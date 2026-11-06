import json
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';a=r/'stage_a_full_axis';o=r/'stage_b2_refined_components';o.mkdir(exist_ok=False);assert json.loads((a/'independent_profile_audit.json').read_text())['status']=='PASS'
p=json.loads((r/'identity_policy_v2.json').read_text());q=pd.read_pickle(a/'all_cells_identity_evidence_PRIVATE.pkl.gz');raw=np.load(a/'marker_counts_PRIVATE.npy',mmap_mode='r');genes=np.load(a/'gene_symbols.npy').tolist();gx={x:i for i,x in enumerate(genes)};av=np.load(a/'source_gene_measured.npy');source=pd.read_csv(a/'source_axis.tsv',sep='\t').source.astype(str);si=pd.Index(source).get_indexer(q.gse_id.astype(str))
for j,rule in enumerate(p['rules']):
 f=np.ones(len(q),bool);e=np.ones(len(q),bool)
 for comp in rule['groups']:
  ii=[gx[x] for x in comp];f&=(raw[:,ii]>=2).any(axis=1);e&=av[:,ii].any(axis=1)[si]
 if rule['parent']=='Tcell':f|=(raw[:,[gx[x] for x in ['CD3D','CD3E','CD3G']]]>=2).sum(axis=1)>=2
 q['positive_'+str(j)]=f;q['evaluable_'+str(j)]=e
q.to_pickle(o/'all_cells_refined_identity_evidence_PRIVATE.pkl.gz')
s=pd.read_csv(a/'all_group_component_summary.tsv',sep='\t').set_index('audit_group');take=q.current_QC_keep
for j,rule in enumerate(p['rules']):
 counts=q.loc[take].groupby('audit_group',observed=True)['positive_'+str(j)].sum();ev=q.loc[take].groupby('audit_group',observed=True)['evaluable_'+str(j)].sum();s['n_positive_'+str(j)]=counts.reindex(s.index).to_numpy();s['n_evaluable_'+str(j)]=ev.reindex(s.index).to_numpy()
s.reset_index().to_csv(o/'all_group_component_summary.tsv',sep='\t',index=False)
old=pd.read_pickle(a/'all_cells_identity_evidence_PRIVATE.pkl.gz');rows=[]
for j in [3,10,18]:
 rows.append({'rule_id':j,'parent':p['rules'][j]['parent'],'v1_positive_QC_keep':int((old['positive_'+str(j)]&take).sum()),'v2_positive_QC_keep':int((q['positive_'+str(j)]&take).sum()),'meaning':'Conservative component refinement after shared-marker review; not QC deletion or negative identity.'})
pd.DataFrame(rows).to_csv(o/'component_refinement_counts.tsv',sep='\t',index=False)
out={'status':'refined_components_complete','UTC':datetime.now(timezone.utc).isoformat(),'n_all_axis':len(q),'n_QC_keep':int(take.sum()),'n_rules':len(p['rules']),'raw_counts_and_original_marker_tables_reused_unchanged':True,'no_QC_or_canonical_labels_changed':True,'changed_rules':[3,10,18]};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
