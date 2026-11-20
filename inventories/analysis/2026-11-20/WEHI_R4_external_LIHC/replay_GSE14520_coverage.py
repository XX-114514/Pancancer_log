import sys,pathlib,gzip,io,json,hashlib,datetime
import pandas as pd
base,bulk,out=map(pathlib.Path,sys.argv[1:])
r=json.loads((out/'qualification.json').read_text())
for name,entry in r['inputs'].items():
 p=bulk/name if name=='candidate_signatures.tsv' else (base/'metadata/platform_annotations'/name if name.startswith('GPL') else base/'raw/cohorts/GSE14520'/name)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
text=gzip.open(base/'metadata/platform_annotations/GPL3921.annot.gz','rt').read()
a=pd.read_csv(io.StringIO(text.split('!platform_table_begin\n')[1].split('!platform_table_end')[0]),sep='\t',dtype=str,keep_default_na=False)
text=gzip.open(base/'raw/cohorts/GSE14520/GSE14520-GPL3921_series_matrix.txt.gz','rt').read()
x=pd.read_csv(io.StringIO(text.split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0]),sep='\t',index_col=0)
assert x.shape==(r['probe_rows'],r['samples'])
a=a[a.ID.isin(x.index)]
actual=pd.read_csv(out/'CM13_probe_coverage.tsv',sep='\t'); n=0
for row in actual.itertuples():
 symbols=a['Gene symbol'].map(lambda v:set(z.strip() for z in v.split('///') if z.strip() not in ('','---')))
 u=a.loc[symbols.map(lambda v:v=={row.gene}),'ID'].nunique()
 b=a.loc[symbols.map(lambda v:row.gene in v and len(v)>1),'ID'].nunique()
 assert (u,b)==(row.unambiguous_probes,row.ambiguous_probes)
 n+=1
result=dict(status='PASS',scope='Independent pandas probe membership and all 120 gene coverage rows; source hashes and expression axes. Not clinical validation.',rows=n,shape=list(x.shape),qualification_sha256=hashlib.sha256((out/'qualification.json').read_bytes()).hexdigest(),time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(out/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
