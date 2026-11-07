"""Pinned CellPhoneDB reconstruction; every formal run uses 1000 permutations."""
import sys, json, time, hashlib
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version
import numpy as np, pandas as pd, anndata as ad
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method as cpdb

r=Path(sys.argv[1]);tag=sys.argv[2];start=time.time();out=r/'results'/tag;out.mkdir(exist_ok=False)
assert json.loads((r/'audit/input_independent_audit.json').read_text())['status'].startswith('PASS_')
if tag in ['pilot','global']:
    a=ad.read_h5ad(r/'inputs/global_1000_PRIVATE.h5ad')
    if tag=='pilot':
        rng=np.random.default_rng(2026110702);ids=[]
        for node in sorted(a.obs.subCluster.unique()):
            loc=np.flatnonzero(a.obs.subCluster.eq(node));ids.extend(rng.choice(loc,min(100,len(loc)),replace=False))
        a=a[np.sort(ids)].copy()
else:
    cm,group=tag.split('_');assert group in ['High','Low'];a=ad.read_h5ad(r/'inputs/full_LR_PRIVATE.h5ad')
    mask=a.obs.CMT.eq(cm).to_numpy();a=a[mask if group=='High' else ~mask].copy()
counts=a.obs.subCluster.value_counts()
if len(a)==0 or len(counts)<2:
    (out/'complete.json').open('x').write(json.dumps({'status':'UNASSESSABLE_FEWER_THAN_TWO_PRESENT_LABELS','n_cells':len(a),'n_labels':len(counts)}));sys.exit(0)
meta=out/'metadata_PRIVATE.tsv';a.obs[['subCluster']].reset_index(names='Cell').rename(columns={'subCluster':'cell_type'}).to_csv(meta,sep='\t',index=False)
n=20 if tag=='pilot' else 1000
result=cpdb.call(cpdb_file_path=str(r/'sources/cellphonedb_v4.1.0.zip'),meta_file_path=str(meta),
    counts_file_path=a,counts_data='hgnc_symbol',iterations=n,threshold=.1,threads=4,
    debug_seed=2026110702,pvalue=.05,subsampling=False,output_path=str(out),output_suffix=tag)
files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.txt')}
assert files
proof={'status':'CPDB_OUTPUT_WRITTEN_PENDING_INDEPENDENT_AUDIT','tag':tag,'UTC':datetime.now(timezone.utc).isoformat(),
    'n_cells':len(a),'n_labels':len(counts),'iterations':n,'threshold':.1,'threads':4,'seed':2026110702,
    'seconds':time.time()-start,'package':version('cellphonedb'),'database':'v4.1.0 reconstructed fixed snapshot',
    'label_cell_counts':counts.to_dict(),'files':files,'patient_or_causal_inference':False,
    'pilot_not_scientific_result':tag=='pilot'}
(out/'complete.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps({k:v for k,v in proof.items() if k!='files'}))
