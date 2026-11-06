"""Read back completed B groups only; a partial receipt never accepts the full family."""
from pathlib import Path
import json,hashlib,pandas as pd,numpy as np,datetime
r=Path(__file__).resolve().parents[1];rows=[]
for receipt_path in sorted((r/'results').glob('CM*/complete.json')):
 d=receipt_path.parent;receipt=json.loads(receipt_path.read_text())
 for name,sha in receipt['sha256'].items():assert hashlib.sha256((d/name).read_bytes()).hexdigest()==sha
 p=pd.read_csv(next(d.glob('statistical_analysis_pvalues_*.txt')),sep='\t').set_index('id_cp_interaction');s=pd.read_csv(next(d.glob('statistical_analysis_significant_means_*.txt')),sep='\t').set_index('id_cp_interaction').reindex(p.index);cols=[c for c in p if '|' in c];assert np.isfinite(p[cols]).all().all();assert ((p[cols]<=.05)==s[cols].notna()).all().all();rows.append({'tag':receipt['tag'],'cells':receipt['cells'],'LR_rows':len(p),'ordered_pairs':len(cols),'candidate_events':int(s[cols].notna().sum().sum()),'seconds':receipt['seconds'],'status':'SHA_and_significance_readback_pass'})
result={'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'groups':rows,'expected_groups':18,'full_B_accepted':False,'not_an_environment_difference_test':True};(r/'audit/B_available_group_readback.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
