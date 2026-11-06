from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
r=Path(__file__).resolve().parents[1];o=r/'results/pilot';done=json.loads((o/'complete.json').read_text());assert done['iterations']==20
for name,sha in done['outputs'].items():assert hashlib.sha256((o/name).read_bytes()).hexdigest()==sha
p=pd.read_csv(next(o.glob('*_pvalues_*.txt')),sep='\t');s=pd.read_csv(next(o.glob('*_significant_means_*.txt')),sep='\t');p=p.set_index('id_cp_interaction');s=s.set_index('id_cp_interaction').reindex(p.index);cols=[c for c in p if '|' in c];v=p[cols].to_numpy();assert np.isfinite(v).all() and v.min()>=0 and v.max()<=1
assert s[cols].shape==v.shape;assert not (s[cols].notna().to_numpy()&(v>.05)).any()
res={'status':'technical_pass','rows':len(p),'ordered_pair_columns':len(cols),'all_output_sha_verified':True,'pvalue_range':[float(v.min()),float(v.max())],'pilot_not_scientific_result':True};(r/'audit/pilot_acceptance.json').write_text(json.dumps(res,indent=2));print(res)
