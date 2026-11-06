import json
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';a=r/'stage_a_full_axis';b=r/'stage_b2_refined_components';q=pd.read_pickle(b/'all_cells_refined_identity_evidence_PRIVATE.pkl.gz');old=pd.read_pickle(a/'all_cells_identity_evidence_PRIVATE.pkl.gz');assert q.index.equals(old.index) and np.array_equal(q.current_QC_keep.to_numpy(),old.current_QC_keep.to_numpy())
p=json.loads((r/'identity_policy_v2.json').read_text());raw=np.load(a/'marker_counts_PRIVATE.npy',mmap_mode='r');genes=np.load(a/'gene_symbols.npy').tolist();gx={g:i for i,g in enumerate(genes)};av=np.load(a/'source_gene_measured.npy');src=pd.Index(pd.read_csv(a/'source_axis.tsv',sep='\t').source.astype(str));si=src.get_indexer(q.gse_id.astype(str))
for left in range(0,len(q),50000):
 right=min(left+50000,len(q));x=raw[left:right];measured=av[si[left:right]]
 for j,rule in enumerate(p['rules']):
  flags=np.logical_and.reduce([(x[:,[gx[g] for g in comp]]>=2).any(axis=1) for comp in rule['groups']]);ev=np.logical_and.reduce([measured[:,[gx[g] for g in comp]].any(axis=1) for comp in rule['groups']])
  if rule['parent']=='Tcell':flags|=(x[:,[gx[g] for g in ['CD3D','CD3E','CD3G']]]>=2).sum(axis=1)>=2
  assert np.array_equal(flags,q['positive_'+str(j)].iloc[left:right].to_numpy()) and np.array_equal(ev,q['evaluable_'+str(j)].iloc[left:right].to_numpy())
 for j in set(range(22))-{3,10,18}:assert np.array_equal(q['positive_'+str(j)].iloc[left:right].to_numpy(),old['positive_'+str(j)].iloc[left:right].to_numpy())
s=pd.read_csv(b/'all_group_component_summary.tsv',sep='\t').set_index('audit_group');t=q.loc[q.current_QC_keep]
for j in range(22):assert np.array_equal(t.groupby('audit_group',observed=True)['positive_'+str(j)].sum().reindex(s.index).to_numpy(),s['n_positive_'+str(j)].to_numpy())
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_all_cells':len(q),'n_rules':22,'only_three_registered_positive_rules_changed':True,'all_positive_and_measurable_components_recomputed':True,'all_group_component_counts_verified':True,'underlying_raw_profile_audit_retained':True,'original_QC_preserved':True};(b/'independent_refinement_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
