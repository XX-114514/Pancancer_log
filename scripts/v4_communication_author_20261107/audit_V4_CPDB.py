"""Read back every real p/mean cellpair output by interaction ID; no effect gating."""
import json, sys, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np, pandas as pd

r=Path(sys.argv[1]);tag=sys.argv[2];out=r/'results'/tag
proof=json.loads((out/'complete.json').read_text())
if proof['status'].startswith('UNASSESSABLE_'):
    (out/'independent_audit.json').open('x').write(json.dumps({'status':'UNASSESSABLE_PRESERVED_NOT_A_PASS','tag':tag}));sys.exit(0)
for name,digest in proof['files'].items():assert hashlib.sha256((out/name).read_bytes()).hexdigest()==digest
def read(kind):
    files=list(out.glob('statistical_analysis_'+kind+'_*.txt'));assert len(files)==1
    t=pd.read_csv(files[0],sep='\t');assert t.id_cp_interaction.is_unique
    return t.set_index('id_cp_interaction')
p=read('pvalues');m=read('means');s=read('significant_means')
assert p.index.equals(m.index) and p.index.equals(s.index)
cols=[c for c in p.columns if '|' in c];assert set(cols)=={a+'|'+b for a in proof['label_cell_counts'] for b in proof['label_cell_counts']}
pv=p[cols].to_numpy(float);means=m[cols].to_numpy(float);sig=s[cols].to_numpy(float)
assert np.isfinite(pv).all() and ((pv>=0)&(pv<=1)).all() and np.isfinite(means).all()
assert np.all(np.isnan(sig)|(np.abs(sig-means)<1e-6))
assert np.array_equal(np.isfinite(sig),pv<=.05)
assert proof['iterations']==(20 if tag=='pilot' else 1000)
counts=s[cols].notna().sum()
pd.DataFrame([{'source':c.split('|')[0],'target':c.split('|')[1],'n_candidate_CCI':int(v)} for c,v in counts.items()]).to_csv(out/'audited_pair_counts.tsv',sep='\t',index=False)
receipt={'status':'PASS_independent_CPDB_hash_ID_pair_scope_p_range_and_significant_mean_readback','UTC':datetime.now(timezone.utc).isoformat(),
    'tag':tag,'n_cellpairs':len(cols),'n_retained_database_interactions':len(p),'n_candidate_cellpair_LR_events':int(np.isfinite(sig).sum()),
    'iterations':proof['iterations'],'pilot_not_scientific_result':tag=='pilot','patient_or_contact_or_causal_claim':False}
(out/'independent_audit.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt))
