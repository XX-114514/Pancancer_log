"""Independent final gate for all 18 CPDB producers; no patient inference."""
from pathlib import Path
import json,hashlib,pandas as pd,numpy as np,datetime
r=Path(__file__).resolve().parents[1];rows=[]
for cm in range(1,10):
 for group in ['High','Low']:
  tag=f'CM{cm:02d}_{group}';d=r/'results'/tag;receipt=json.loads((d/'complete.json').read_text())
  for name,sha in receipt['sha256'].items():assert hashlib.sha256((d/name).read_bytes()).hexdigest()==sha
  p=pd.read_csv(next(d.glob('statistical_analysis_pvalues_*.txt')),sep='\t').set_index('id_cp_interaction');s=pd.read_csv(next(d.glob('statistical_analysis_significant_means_*.txt')),sep='\t').set_index('id_cp_interaction').reindex(p.index);cols=[x for x in p if '|' in x];assert p.index.is_unique;assert np.isfinite(p[cols]).all().all();assert ((p[cols]>=0)&(p[cols]<=1)).all().all();assert ((p[cols]<=.05)==s[cols].notna()).all().all();rows.append({'tag':tag,'cells':receipt['cells'],'LR_rows':len(p),'ordered_pairs':len(cols),'candidate_events':int(s[cols].notna().sum().sum()),'status':'producer_sha_and_significance_readback_pass'})
pd.DataFrame(rows).to_csv(r/'audit/B_all18_readback.tsv',sep='\t',index=False)
(r/'audit/B_all18_acceptance.json').write_text(json.dumps({'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_groups':len(rows),'status':'technical_and_descriptive_readback_pass','patient_or_causal_inference':False},indent=2))
