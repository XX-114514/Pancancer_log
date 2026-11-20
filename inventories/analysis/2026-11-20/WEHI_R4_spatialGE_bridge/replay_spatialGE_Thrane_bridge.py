import pathlib,sys,zipfile,hashlib,json,io,datetime
import pandas as pd
import numpy as np
root,out=map(pathlib.Path,sys.argv[1:]);r=json.loads((out/'complete.json').read_text())
z=zipfile.ZipFile(root/'spatialGE_Thrane_v1/melanoma_thrane.zip');a=zipfile.ZipFile(root/'Thrane2018_histology_v1/Thrane_et_al_2018_CAN_RES.zip')
read=lambda b:pd.read_csv(io.BytesIO(b),sep='\t',index_col=0)
values=spots=0
for row in r['samples']:
 sample=row['sample'];nb=z.read('melanoma_thrane/'+sample+'_counts.tsv');ob=a.read('Thrane_et_al_2018_CAN_RES/Count_data/'+sample+'_counts.tsv')
 assert hashlib.sha256(nb).hexdigest()==row['derivative_count_sha256'] and hashlib.sha256(ob).hexdigest()==row['original_count_sha256']
 n,o=read(nb),read(ob)
 assert np.array_equal(n.to_numpy(),o.to_numpy()) and list(n.columns)==list(o.columns)
 assert list(n.index)==list(o.index.str.split().str[0])
 m=read(z.read('melanoma_thrane/'+sample+'_mapping.tsv')).loc[n.columns]
 expected=np.array([list(map(int,s.split('x'))) for s in n.columns])
 assert np.array_equal(m[['xpos','ypos']].to_numpy(),expected)
 assert n.size==row['numeric_values'] and n.shape==(row['genes'],row['spots'])
 assert int(n.index.duplicated().sum())==row['duplicate_symbols']
 values+=n.size;spots+=n.shape[1]
result=dict(status='PASS',scope='All four expression arrays, symbols, duplicate counts and barcode coordinate values; does not verify primary clinical identity or image registration',numeric_values=values,spots=spots,slices=len(r['samples']),complete_sha256=hashlib.sha256((out/'complete.json').read_bytes()).hexdigest(),time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(out/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
