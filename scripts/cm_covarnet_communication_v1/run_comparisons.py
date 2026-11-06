from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,time,json,hashlib
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method as cpdb
r=Path(__file__).resolve().parents[1];a=ad.read_h5ad(r/'inputs/full_LR.h5ad');cms=[f'CM{i:02d}' for i in range(1,10)]
for cm in cms:
 for group in ['High','Low']:
  tag=cm+'_'+group;o=r/'results'/tag;o.mkdir(exist_ok=True)
  if (o/'complete.json').exists():continue
  mask=a.obs.CMT==cm if group=='High' else a.obs.CMT!=cm;b=a[mask].copy();meta=o/'meta_PRIVATE.tsv';b.obs[['identity']].reset_index(names='Cell').rename(columns={'identity':'cell_type'}).to_csv(meta,sep='\t',index=False)
  t=time.time();print('START',tag,b.shape,flush=True)
  res=cpdb.call(cpdb_file_path=str(r/'sources/cellphonedb_v4.1.0.zip'),meta_file_path=str(meta),counts_file_path=b,counts_data='hgnc_symbol',output_path=str(o),iterations=1000,threshold=.1,threads=4,debug_seed=1337,pvalue=.05,subsampling=False,output_suffix=tag)
  receipt={'tag':tag,'cells':len(b),'seconds':time.time()-t,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in o.glob('*.txt')}};(o/'complete.json').write_text(json.dumps(receipt,indent=2));print('COMPLETE',tag,flush=True)
