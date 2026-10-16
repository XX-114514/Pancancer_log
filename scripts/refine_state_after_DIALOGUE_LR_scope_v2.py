"""Refine one broad missing-score label; preserve all v1 effects byte-for-byte."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time
import pandas as pd
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=V/'state_after_DIALOGUE_LR_response_20261016_v1';start=time.monotonic();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
while not (P/'inventory.json').exists():
 e=json.loads((V/'state_after_DIALOGUE_LR_response_20261016_v1.execution.json').read_text());p=Path('/proc',str(e['child_pid']),'cmdline');assert p.exists() and e['producer'].encode() in p.read_bytes(),e;assert time.monotonic()-start<36*3600;time.sleep(20)
iv=json.loads((P/'inventory.json').read_text())
for f,h in iv['files'].items():assert sha(P/f)==h
scope=pd.read_csv(V/'state_after_DIALOGUE_frozen_projection_20261016_v1/public_slot_scope.tsv',sep='\t');t=pd.read_csv(P/'all258944_hypothesis_scope_counts.tsv',sep='\t');changed=[]
for i,r in t[t.status.eq('NO_GENE_SUPPORTED_RECEIVER_SCORE')].iterrows():
 root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{int(r.context):04d}';ci=json.loads((root/'private_input_inventory.json').read_text());assert r.target in ci['endpoints'];typ='type'+['A','B'][ci['endpoints'].index(r.target)];q=scope[(scope.context==r.context)&scope.type.eq(typ)&scope.program.eq(r.MCP)];assert len(q)==1
 if q.iloc[0].status=='NO_GENE_SUPPORTED_SCORE':new='NO_GENE_SUPPORTED_RECEIVER_SCORE'
 elif q.iloc[0].status=='FAILED_PARITY_OR_DESIGN_REVIEW':new='RECEIVER_SCORE_RECONSTRUCTION_OR_DESIGN_FAILED'
 else:raise AssertionError((r.context,typ,r.MCP,q.iloc[0].status))
 if new!=r.status:changed.append(dict(context=int(r.context),source=r.source,target=r.target,MCP=r.MCP,records=int(r.records),old_status=r.status,new_status=new));t.loc[i,'status']=new
assert int(t.records.sum())==258944;O=V/'state_after_DIALOGUE_LR_scope_detail_20261016_v2';assert not O.exists();O.mkdir(mode=0o700);t.to_csv(O/'all258944_hypothesis_scope_counts.tsv',sep='\t',index=False)
for f in ['all_estimable_exploratory_LR_receiver_associations.tsv','independent_scalar_OLS_checks.tsv']:
 (O/f).write_bytes((P/f).read_bytes());assert sha(O/f)==iv['files'][f]
(O/'inventory.json').write_text(json.dumps(dict(status='EFFECTS_UNCHANGED_MISSING_RECEIVER_SCORE_REASON_REFINED',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),parent_inventory_sha256=sha(P/'inventory.json'),projection_slot_scope_sha256=sha(V/'state_after_DIALOGUE_frozen_projection_20261016_v1/public_slot_scope.tsv'),changed_scope_rows=changed,all_effects_p_CI_BH_byte_identical=True,files={p.name:sha(p) for p in O.iterdir() if p.is_file()}),indent=2)+'\n')
