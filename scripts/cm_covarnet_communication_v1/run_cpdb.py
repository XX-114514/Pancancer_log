from pathlib import Path
import json,sys,time,hashlib,importlib.metadata
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method as cpdb
r=Path(__file__).resolve().parents[1];mode=sys.argv[1];assert mode in ['pilot','global']
out=r/'results'/mode;out.mkdir(exist_ok=True);done=out/'complete.json'
if done.exists():print('already complete');sys.exit(0)
counts=r/'inputs/global_1000.h5ad';meta=r/'inputs/global_meta.tsv'
if mode=='pilot':
 import anndata as ad,pandas as pd,numpy as np
 a=ad.read_h5ad(counts);rng=np.random.default_rng(1337);ids=[]
 for _,g in a.obs.groupby('identity',observed=True):ids.extend(rng.choice(g.index,min(100,len(g)),replace=False))
 a=a[ids].copy();counts=r/'inputs/pilot.h5ad';a.write_h5ad(counts);meta=r/'inputs/pilot_meta.tsv';a.obs[['identity']].reset_index(names='Cell').rename(columns={'identity':'cell_type'}).to_csv(meta,sep='\t',index=False)
t=time.time();res=cpdb.call(cpdb_file_path=str(r/'sources/cellphonedb_v4.1.0.zip'),meta_file_path=str(meta),counts_file_path=str(counts),counts_data='hgnc_symbol',output_path=str(out),iterations=20 if mode=='pilot' else 1000,threshold=0.1,threads=4,debug_seed=1337,pvalue=0.05,subsampling=False,output_suffix=mode)
receipt={'mode':mode,'iterations':20 if mode=='pilot' else 1000,'seconds':time.time()-t,'package':importlib.metadata.version('cellphonedb'),'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.txt')},'scientific_acceptance':False}
assert receipt['outputs'];done.write_text(json.dumps(receipt,indent=2));print(receipt)
