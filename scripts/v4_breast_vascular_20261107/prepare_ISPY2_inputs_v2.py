"""Outcome-blind measurement and exact-key assembly for an I-SPY2 proxy analysis."""
from pathlib import Path
import os,json,hashlib,gzip,collections
from datetime import datetime,timezone
import pandas as pd,numpy as np
w=Path(os.environ['BREAST_WORK']);source=Path(os.environ['ISPY2_SOURCE']);signatures=Path(os.environ['V4_SIGNATURES']);out=w/'clinical_v2';out.mkdir(exist_ok=False)
receipt=json.loads((source/'open_source_download_receipt_PRIVATE.json').read_text())
for z in receipt['files']:
 p=source/z['file'];assert p.stat().st_size==z['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256']
meta=json.loads((source/'parsed_metadata_PRIVATE.json').read_text());c=collections.Counter(x['patient id'] for x in meta);duplicates={k for k,n in c.items() if n>1}
assert all(x['tissue']=='Breast cancer biopsy (pre-treatment)' for x in meta)
m=pd.DataFrame([x for x in meta if x['patient id'] not in duplicates]);m.index=m['patient id'];assert m.index.is_unique
eligible=m['her2'].eq('0')&m['arm'].isin(['Paclitaxel + Pembrolizumab','Paclitaxel'])&m['hr'].isin(['0','1'])&m['pcr'].isin(['0','1']);m=m[eligible].copy()
f=source/'GSE194040_geneLevel_PRIVATE.txt.gz'
with gzip.open(f,'rt') as q:header=q.readline().rstrip('\n\r').split('\t')
expr=pd.read_csv(f,sep='\t',header=None,skiprows=1,index_col=0,float_precision='round_trip');assert expr.shape[1]==len(header);expr.columns=header
exprdup={k for k,n in collections.Counter(header).items() if n>1};m=m[~m.index.isin(exprdup)];missing=set(m.index)-set(expr.columns);assert not missing, 'Clinical keys absent from native expression'
x=expr.loc[:,m.index].T
program=pd.read_csv(signatures,sep='\t');panels=json.loads((w/'frozen_breast_plan_v1.json').read_text())['clinical_marker_panels'];panels['panT']=['CD3D','CD3E'];groups={**{cm:program.loc[program.cm.eq(cm),'gene'].tolist() for cm in ['CM03','CM04','CM07','CM09','CM11']},**panels}
scores={};coverage=[]
for k,genes in groups.items():
 assert len(genes)==(10 if k.startswith('CM') else 2 if k=='panT' else 3)
 measured=[g for g in genes if sum(x.columns==g)==1];coverage.append({'program':k,'required_genes':len(genes),'unique_measured_genes':len(measured),'missing_or_ambiguous':';'.join(g for g in genes if g not in measured)})
 if len(measured)!=len(genes):scores[k]=pd.Series(np.nan,index=x.index);continue
 z=x.loc[:,genes];ranks=z.rank(method='average',pct=True,na_option='keep');scores[k]=ranks.mean(axis=1).where(z.notna().all(axis=1))
s=pd.DataFrame(scores);s['Vascular_B_NK_equal_mean']=s[['B','EC','NK']].mean(axis=1).where(s[['B','EC','NK']].notna().all(axis=1));s['Vascular_B_NK_min_component']=s[['B','EC','NK']].min(axis=1).where(s[['B','EC','NK']].notna().all(axis=1));s['mean_B_NK']=s[['B','NK']].mean(axis=1).where(s[['B','NK']].notna().all(axis=1))
s.to_csv(out/'unscaled_scores_PRIVATE.csv',index_label='patient_key');z=(s-s.mean())/s.std(ddof=1);assert np.isfinite(z.dropna(axis=1).to_numpy()).all()
z['HR']=m.hr.astype(int);z['pCR']=m.pcr.astype(int);z['treatment']=m.arm.eq('Paclitaxel + Pembrolizumab').astype(int);z.to_csv(out/'analysis_inputs_PRIVATE.csv',index_label='patient_key')
pd.DataFrame(coverage).to_csv(out/'program_measurement_coverage.tsv',sep='\t',index=False)
aggregate=[]
for a,t in m.groupby('arm'):
 aggregate.append({'arm':a,'n':len(t),'n_HR_positive':int(t.hr.eq('1').sum()),'n_TNBC':int(t.hr.eq('0').sum()),'n_pCR':int(t.pcr.eq('1').sum())})
pd.DataFrame(aggregate).to_csv(out/'cohort_eligibility_aggregate.tsv',sep='\t',index=False)
q={'UTC':datetime.now(timezone.utc).isoformat(),'status':'EXACT_KEYS_AND_OUTCOME_BLIND_SCORE_MEASUREMENT_COMPLETE','n_source_arrays':len(header),'n_unique_expression_keys':len(set(header)),'n_metadata_records':len(meta),'n_unique_metadata_patient_keys':len(c),'n_duplicated_metadata_keys':len(duplicates),'n_duplicate_expression_keys':len(exprdup),'n_eligible_patients':len(m),'duplicate_policy':'all duplicate keys excluded before arm selection, no averaging or suffix stripping','missing_eligible_expression_keys':len(missing),'n_native_genes':len(expr),'source_files':[{'source':v['source'],'bytes':v['bytes'],'sha256':v['sha256']} for v in receipt['files']],'plan_sha256':hashlib.sha256((w/'frozen_breast_plan_v1.json').read_bytes()).hexdigest(),'signature_sha256':hashlib.sha256(signatures.read_bytes()).hexdigest(),'private_analysis_input_sha256':hashlib.sha256((out/'analysis_inputs_PRIVATE.csv').read_bytes()).hexdigest(),'patient_scores_not_public':True,'measurement_addendum_sha256':hashlib.sha256((w/'measurement_addendum_before_estimable_effects_v2.json').read_bytes()).hexdigest()}
(out/'input_measurement_receipt.json').write_text(json.dumps(q,indent=2)+'\n');print(json.dumps(q),flush=True)
