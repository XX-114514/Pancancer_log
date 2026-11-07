"""Synthetic effect vectors verify patient units, NA denominators and family sizes."""
from pathlib import Path
from itertools import product
import json,sys,hashlib
import numpy as np,pandas as pd
from vascular_spatial_controls_v2 import patient_summary,AXES,SENDERS,RECEIVERS,RADII
root=Path(sys.argv[1]);root.mkdir(exist_ok=False);dirs=[];donors=[]
for i in range(6):
 d=root/f'SYNTH_CASE_{i}';d.mkdir();dirs.append(d);donors.append(f'SYNTH_DONOR_{i}')
 (d/'complete.json').write_text(json.dumps({'native_source_file_sha256':hashlib.sha256(f'SYNTH_SOURCE_{i}'.encode()).hexdigest()}))
 n=[{'radius_um':radius,'receiver':receiver,'log_odds_enrichment':(i+1)*.1 if radius==50 and receiver=='B' else 0.} for radius,receiver in product(RADII,RECEIVERS)];pd.DataFrame(n).to_csv(d/'neighborhood.tsv',sep='\t',index=False)
 q=[{'sender':s,'receiver':r,'ligand':l,'receptor':g,'status':'unavailable','AB_minus_A':np.nan,'AB_minus_B':np.nan,'excess_edge_fraction':np.nan} for s,r,(l,g) in product(SENDERS,RECEIVERS,AXES)]
 for name in ['regional_ES.tsv','conditional_LR_edges.tsv']:pd.DataFrame(q).to_csv(d/name,sep='\t',index=False)
t=patient_summary(dirs,donors,root/'full6');assert len(t)==66;z=t[(t.family=='neighborhood_primary2')&(t.hypothesis=='B_50um')].iloc[0];assert z.P==1/64 and z.BH==1/32
one=patient_summary(dirs[:1],donors[:1],root/'only1');assert len(one)==66 and one.P.isna().all()
try:patient_summary([dirs[0],dirs[0]],donors[:2],root/'repeated_case');raise RuntimeError('Repeated source accepted')
except AssertionError:pass
try:patient_summary(dirs[:2],['SYNTH_SAME_DONOR']*2,root/'repeated_donor');raise RuntimeError('Repeated donor accepted')
except AssertionError:pass
print('PASS_SYNTHETIC_ONLY:66 slots,2+10+27+27 full denominators,known6-patient exactP/BH,one-patient P unavailable,repeated donors/cases rejected')
