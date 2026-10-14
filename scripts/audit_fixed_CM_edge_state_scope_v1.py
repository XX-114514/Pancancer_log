"""Long-form dictionary independent check of all source/stratum edge denominators."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'fixed_CM_edge_state_scope_preflight_20261014_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();iv=json.loads((O/'inventory.json').read_text());assert iv['unique_pairs']==60 and iv['original_CM_edge_records']==120
for name,h in iv['output_sha256'].items():assert sha(O/name)==h
I=R/'all_CM_state_native_input_preflight_20261014_v1';u=pd.read_csv(I/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False);meta=u[['sample_uid','gse_id','stratum','effective_resampling_block','unit_kind','analysis_eligible']].drop_duplicates().set_index('sample_uid');assert len(meta)==1290 and meta.index.is_unique;major_counts=u.groupby(['sample_uid','major']).n_cells.sum().to_dict();node_counts=u.set_index(['sample_uid','subCluster']).n_cells.to_dict();pairs=pd.read_csv(O/'all60_pair_scope_summary.tsv',sep='\t');records=pd.read_csv(O/'all_unique_pair_context_capture.tsv',sep='\t');mapping=pd.read_csv(O/'original120_to_unique60_pair_mapping.tsv',sep='\t');assert len(records)==5820 and len(mapping)==120
for row in records.itertuples():
 m=meta[meta.stratum.eq(row.stratum)];eligible=m[m.analysis_eligible.astype(str).eq('True')];cov=[];state=[];fra=[];both_all=0;fine_all=0
 for sid,r in m.iterrows():
  da=major_counts.get((sid,row.major_A),0);db=major_counts.get((sid,row.major_B),0);na=node_counts.get((sid,row.endpoint_A),0);nb=node_counts.get((sid,row.endpoint_B),0);both=da>=20 and db>=20;fine=na>=20 and nb>=20;both_all+=both;fine_all+=fine
  if sid in eligible.index and both:cov.append(sid);fra.append((na/da,nb/db,na,nb))
  if sid in eligible.index and fine:state.append(sid)
 variation=bool(fra and max(x[0] for x in fra)>min(x[0] for x in fra) and max(x[1] for x in fra)>min(x[1] for x in fra));assert row.composition_primary_libraries==len(cov) and row.state_measurement_primary_libraries==len(state) and row.total_libraries==len(m) and row.parent_major_capture20_libraries_allclasses==both_all and row.fine_both20_libraries_allclasses==fine_all;assert row.composition_variable_proportions==variation and row.composition_potential_testable==(len(cov)>=10 and variation) and row.state_potential_atleast20_libraries==(len(state)>=20);assert row.composition_zero_A_libraries==sum(x[2]==0 for x in fra) and row.composition_zero_B_libraries==sum(x[3]==0 for x in fra)
 for samples,prefix in [(cov,'composition'),(state,'state')]:
  mm=meta.loc[samples];reviewed=mm[mm.unit_kind.eq('reviewed_patient')];assert getattr(row,prefix+'_effective_blocks')==mm.effective_resampling_block.nunique() and getattr(row,prefix+'_reviewed_libraries')==len(reviewed) and getattr(row,prefix+'_reviewed_blocks')==reviewed.effective_resampling_block.nunique()
for r in pairs.itertuples():
 q=records[records.unique_pair_id.eq(r.unique_pair_id)];mp=mapping[mapping.unique_pair_id.eq(r.unique_pair_id)];assert r.original_edge_records==len(mp) and set(r.CMs.split(';'))==set(mp.cm);assert r.primary_composition_potential_contexts==int(q.composition_potential_testable.sum()) and r.primary_state20_potential_sources==q.loc[q.state_potential_atleast20_libraries,'gse_id'].nunique() and r.primary_state20_potential_contexts==int(q.state_potential_atleast20_libraries.sum())
for r in mapping.itertuples():
 a,b=sorted([r.original_endpoint1,r.original_endpoint2]);pid='PAIR_'+hashlib.sha256((a+'\0'+b).encode()).hexdigest()[:12];assert a==r.canonical_endpoint_A and b==r.canonical_endpoint_B and pid==r.unique_pair_id
receipt=dict(status='PASS_ALL5820_PAIR_CONTEXT_DENOMINATORS_60SUMMARIES_AND120_MAPS',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_inventory_sha256=sha(O/'inventory.json'),audit_code_sha256=sha(__file__),pair_context_records=5820,unique_pairs=60,CM_edge_records=120,independent_long_form_count_lookup=True,zero_values_block_reviewed_counts_and_dedup_checked=True,score_or_design_estimability_not_tested=True,formal_states=False)
A=R/'fixed_CM_edge_state_scope_audit_20261014_v1';assert not A.exists();A.mkdir();(A/'independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
