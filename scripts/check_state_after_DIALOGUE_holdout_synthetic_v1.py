from pathlib import Path
import tempfile,subprocess,sys
import numpy as np,pandas as pd
w=Path(tempfile.mkdtemp(prefix='dialogue_fixed_holdout_synthetic_'));p=w/'projection';m=w/'models';p.mkdir();m.mkdir();scope=[];scores=[]
for c in range(14):
 for prog in ['MCP1','MCP2']:
  for side in ['typeA','typeB']:
   scope.append(dict(context=c,program=prog,type=side,status='NO_GENE_SUPPORTED_SCORE' if c==1 else 'FROZEN_RECONSTRUCTION_PASS'))
   for part,n in [('discovery',20),('holdout',12)]:
    for i in range(n):
     v=float(i) if side=='typeA' else float((n-1-i) if prog=='MCP2' else i)
     if c==0 and side=='typeB':v=1.
     scores.append(dict(context=c,type=side,program=prog,partition=part,sample='b'+str(i),score=v))
pd.DataFrame(scope).to_csv(p/'public_slot_scope.tsv',sep='\t',index=False);pd.DataFrame(scores).to_csv(p/'private_block_mean_scores.tsv',sep='\t',index=False);pd.DataFrame(dict(context=list(range(14)))).to_csv(m/'public_fit_status.tsv',sep='\t',index=False)
subprocess.run([sys.executable,'/tmp/test_state_after_DIALOGUE_holdout_v1.py',str(p),str(m),str(w/'result')],check=True)
t=pd.read_csv(w/'result/public_internal_holdout_results.tsv',sep='\t');assert len(t)==28;assert t[t.context.isin([0,1])].p.isna().all();v=t[t.context>=2];assert len(v)==24;assert np.allclose(v[v.program=='MCP1'].holdout_rho,1);assert np.allclose(v[v.program=='MCP2'].holdout_rho,-1);assert (v.p<.001).all();assert np.allclose(v.BH_q,v.p);assert (v.valid_bootstrap==2000).all();print('PASS_28_SLOTS_POSITIVE_NEGATIVE_CONSTANT_MISSING_BH_BLOCK_BOOTSTRAP')
