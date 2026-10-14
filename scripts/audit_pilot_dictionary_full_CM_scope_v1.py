"""Independent np.isin gene audit and source sample-pivot CM denominators."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);O=R/'pilot_dictionary_full_CM_scope_preflight_20261014_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();iv=json.loads((O/'inventory.json').read_text());assert iv['status'].startswith('COMPLETE_FROZEN_PILOT_GENES')
for name,h in iv['output_sha256'].items():assert sha(O/name)==h
I=R/'all_CM_state_native_input_preflight_20261014_v1';u=pd.read_csv(I/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False);p=json.loads((R/'state_programme_dictionary_20261014_v2/frozen_programme_dictionary.json').read_text());programmes={x['programme_id']:x for x in p['programmes']};graph=pd.read_csv(os.environ['PANCANCER_CM_CORE_GRAPH'],sep='\t');axis=pd.read_csv(I/'all_source_native_scope.tsv',sep='\t');c=pd.read_csv(O/'all52_source_pilot_programme_gene_capture.tsv',sep='\t',keep_default_na=False);m=pd.read_csv(O/'all9CM_by52source_native_capture.tsv',sep='\t');memb=pd.read_csv(O/'all_CM_core_member_source_capture.tsv',sep='\t');checks=0
for sid,rows in c.groupby('gse_id',sort=True):
 gp=AR/'20261010_uniform_expression_annotation_v1/by_gse'/sid/'genes.tsv';assert sha(gp)==axis.set_index('gse_id').loc[sid,'native_gene_axis_sha256'];g=pd.read_csv(gp,sep='\t').gene.to_numpy(str);su=u[u.gse_id.eq(sid)]
 for r in rows.itertuples():
  pp=programmes[r.programme_id];names=np.asarray(pp['primary_genes']);mask=np.isin(names,g);n=int(mask.sum());assert r.native_measurable_genes==n and r.signature_genes==len(names) and np.isclose(r.native_gene_fraction,n/len(names));assert r.native_gene_axis_gate==(n>=10 and n>=.8*len(names));assert set(filter(None,r.signature_missing_genes.split(';')))==set(names[~mask]);q=su[su.major.eq(r.major)];assert r.cells==int(q.n_cells.sum()) and r.libraries==q.sample_uid.nunique() and r.original_nodes==q.subCluster.nunique() and r.sample_node_units_at_least20==int(q.n_cells.ge(20).sum());checks+=1
 # Each source cell count table gets an independent sample-by-frozen-node pivot.
 pv=su.pivot(index='sample_uid',columns='subCluster',values='n_cells').fillna(0)
 for row in m[m.gse_id.eq(sid)].itertuples():
  gg=graph[graph.cm.eq(row.cm)];cols=gg.node.to_list();a=pv.reindex(columns=cols,fill_value=0).to_numpy();w=gg.core_weight.to_numpy();anycaptured=(a>0).any(axis=0);eligible=(a>=20).any(axis=0);same=(a>=20).sum(axis=1)>=2;assert row.core_members==len(cols) and row.members_captured_any_cells==int(anycaptured.sum()) and row.members_at_least20_in_any_library==int(eligible.sum()) and row.libraries_with_at_least2_members20==int(same.sum());assert np.isclose(row.captured_core_weight_any_cells,w[anycaptured].sum()) and np.isclose(row.captured_core_weight_at_least20,w[eligible].sum());pilot=gg.major.isin(['T cell','Macrophage']).to_numpy();assert row.pilot_T_Mac_members_at_least20==int((eligible&pilot).sum()) and np.isclose(row.pilot_T_Mac_core_weight_at_least20,w[eligible&pilot].sum())
  for j,node in enumerate(cols):
   mr=memb[(memb.gse_id==sid)&(memb.cm==row.cm)&(memb.node==node)];assert len(mr)==1;r=mr.iloc[0];assert r.cells==int(a[:,j].sum()) and r.libraries==int((a[:,j]>0).sum()) and r.libraries_at_least20==int((a[:,j]>=20).sum())
assert checks==len(c)==52*99 and len(m)==52*9 and len(memb)==52*82
receipt=dict(status='PASS_ALL5148_PILOT_GENE_CAPTURE_ALL468_CM_SOURCE_AND4264_MEMBER_DENOMINATORS',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_inventory_sha256=sha(O/'inventory.json'),audit_code_sha256=sha(__file__),gene_records=checks,CM_source_records=len(m),member_source_records=len(memb),independent_gene_isin_and_sample_node_pivot=True,expression_activity_not_tested=True,no_LR_or_clinical_test=True,formal_states=False)
A=R/'pilot_dictionary_full_CM_scope_audit_20261014_v1';assert not A.exists();A.mkdir();(A/'independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
