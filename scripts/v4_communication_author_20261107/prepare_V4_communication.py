"""Same discovery cells, native counts and measured-gene mask; no CM re-fitting."""
import json, hashlib, time, shutil, zipfile
from pathlib import Path
from datetime import datetime, timezone
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io, anndata as ad, zarr, h5py
from anndata.io import read_elem

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1'
base=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2'
old=PROJECT_ROOT/'runs/20261106_V5_K9_CoVarNet_communication_v1'
start=time.time()
for folder in ['inputs','private','audit','results','sources','config','launch','scripts']:(r/folder).mkdir(exist_ok=True)
assert json.loads((base/'CM_results/FINAL_INDEPENDENT_AUDIT.json').read_text())['status'].startswith('PASS_')
q=pd.read_pickle(base/'inputs/selected_cell_annotation_and_sample_PRIVATE.pkl.gz')
H=pd.read_csv(base/'CM_results/H_PRIVATE.tsv',sep='\t',index_col=0)
W=pd.read_csv(base/'CM_results/W_PRIVATE.tsv',sep='\t',index_col=0)
assert len(q)==560871 and q.subCluster.nunique()==45 and q.sampleID.nunique()==164
assert q.index.is_unique and set(q.sampleID)==set(H.columns) and H.shape==(14,164)
assert np.isfinite(H.to_numpy()).all() and (H>=0).all().all() and H.sum().gt(0).all()
ties=np.isclose(H,H.max().to_numpy(),rtol=1e-10,atol=1e-12).sum(0)>1
assert not ties.any()
cmt=H.idxmax();q=q.copy();q['CMT']=q.sampleID.map(cmt)
nodes=pd.read_csv(base/'CM_results/each_node.tsv',sep='\t')
nodes.to_csv(r/'audit/network_members.tsv',sep='\t',index=False)
source=q.drop_duplicates('sampleID').set_index('sampleID').gse_id
pd.crosstab(cmt,source).rename_axis('cm').to_csv(r/'audit/CMT_source_sample_counts.tsv',sep='\t')
elig=[]
for cm in H.index:
    ids=set(nodes.loc[nodes.cm.eq(cm),'subCluster']);hi=q.CMT.eq(cm);nh=int(cmt.eq(cm).sum())
    sources=source[cmt.eq(cm)].value_counts()
    elig.append({'cm':cm,'n_network_members':len(ids),'n_high_samples':nh,'n_low_samples':164-nh,
                 'n_high_sources':len(sources),'largest_high_source_fraction':float(sources.max()/nh),
                 'n_high_cells':int(hi.sum()),'n_low_cells':int((~hi).sum()),
                 'status':'ORIGINAL_EMPTY_NETWORK_NO_WITHIN_CM_TEST' if not ids else 'EXPLORATORY_NO_SPATIAL_GATE_ACCEPTANCE_YET'})
pd.DataFrame(elig).to_csv(r/'audit/CM_comparison_eligibility.tsv',sep='\t',index=False)
q.groupby(['CMT','subCluster'],observed=True).size().reset_index(name='n_cells').to_csv(r/'audit/CMT_label_cell_counts.tsv',sep='\t',index=False)
root=PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects'
g=zarr.open_group(str(root/'pancancer_core_fullgene_v1.zarr'),mode='r')
axis=pd.Index(read_elem(g['obs/_index']));genes=pd.Index(read_elem(g['var/_index']))
q['raw_row']=axis.get_indexer(q.index);assert q.raw_row.ge(0).all()
q=q.sort_values('raw_row');positions=q.raw_row.to_numpy();assert np.all(np.diff(positions)>0)
available={}
for s in sorted(q.gse_id.unique()):
    with h5py.File(root/'harmonized_per_gse'/(s+'_symbol_counts.h5ad'),'r') as h:available[s]=set(read_elem(h['var']).index)
common=sorted(set.intersection(*available.values()))
assert len(common)>10000 and set(common)<=set(genes)
pd.DataFrame({'source':list(available),'n_native_genes':[len(v) for v in available.values()],
              'n_common_measured_genes':len(common)}).to_csv(r/'audit/native_measurement_scope.tsv',sep='\t',index=False)
for name,digest in [('cellphonedb_v4.1.0.zip','bac6dc33a959f39a79465e0d09b556094da4e154825e50f1470b909a6a62188c'),
                    ('signatures_one_to_one.tsv.gz','1317f2eaede1afcdc466dddce5091c97384835e1e6a5cafeeef8fcf3e2985671')]:
    p=old/'sources'/name;assert hashlib.sha256(p.read_bytes()).hexdigest()==digest
    shutil.copyfile(p,r/'sources'/name)
with zipfile.ZipFile(r/'sources/cellphonedb_v4.1.0.zip') as db:dbgenes=pd.read_csv(db.open('gene_table.csv')).hgnc_symbol.dropna().unique()
pd.DataFrame({'gene':dbgenes,'measured_in_all_13_sources':pd.Index(dbgenes).isin(common)}).to_csv(r/'audit/CPDB_gene_measurement.tsv',sep='\t',index=False)
chosen_cols=genes.get_indexer(common);ptr=g['X/indptr'][:];blocks=[]
for left in range(0,len(axis),20000):
    right=min(left+20000,len(axis));lo=np.searchsorted(positions,left);hi=np.searchsorted(positions,right)
    if hi==lo:continue
    a,b=int(ptr[left]),int(ptr[right]);data=g['X/data'][a:b];ind=g['X/indices'][a:b]
    assert np.isfinite(data).all() and (data>=0).all() and np.allclose(data,np.rint(data))
    raw=sp.csr_matrix((data,ind,ptr[left:right+1]-a),shape=(right-left,len(genes)))[positions[lo:hi]-left]
    total=np.asarray(raw.sum(1)).ravel();assert (total>0).all()
    x=raw[:,chosen_cols].astype(np.float32)
    x.data*=np.repeat(10000/total,np.diff(x.indptr));np.log1p(x.data,out=x.data)
    blocks.append(x)
    (r/'inputs/progress.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),
        'raw_rows_scanned':right,'selected_rows_completed':int(hi),'selected_rows_total':len(q)}))
X=sp.vstack(blocks,format='csr');assert X.shape==(len(q),len(common))
obs=q[['sampleID','gse_id','subCluster','original_major','majorCluster','CMT']].copy()
a=ad.AnnData(X,obs=obs,var=pd.DataFrame(index=pd.Index(common,name='gene')));a.uns['log1p']={'base':None}
a.write_h5ad(r/'inputs/full_common_PRIVATE.h5ad',compression='gzip')
lr=ad.AnnData(X[:,pd.Index(common).isin(dbgenes)].copy(),obs=obs.copy(),var=a.var.loc[pd.Index(common).isin(dbgenes)].copy())
lr.uns['log1p']={'base':None};lr.write_h5ad(r/'inputs/full_LR_PRIVATE.h5ad',compression='gzip')
rng=np.random.default_rng(2026110702);selected=[]
for node in sorted(obs.subCluster.unique()):
    ids=np.flatnonzero(obs.subCluster.eq(node));selected.extend(rng.choice(ids,min(1000,len(ids)),replace=False))
selected=np.sort(selected);glob=a[selected].copy();glob.write_h5ad(r/'inputs/global_1000_PRIVATE.h5ad',compression='gzip')
glob.obs[['subCluster']].reset_index(names='Cell').rename(columns={'subCluster':'cell_type'}).to_csv(r/'inputs/global_meta_PRIVATE.tsv',sep='\t',index=False)
scipy.io.mmwrite(r/'inputs/CellChat_logexpression_PRIVATE.mtx',glob.X.T)
pd.Series(glob.var_names).to_csv(r/'inputs/CellChat_genes_PRIVATE.txt',index=False,header=False)
glob.obs[['subCluster']].to_csv(r/'inputs/CellChat_meta_PRIVATE.tsv',sep='\t')
q.to_pickle(r/'private/cell_bridge_PRIVATE.pkl.gz',compression={'method':'gzip','compresslevel':1})
# Keep sample count and source context, without pretending samples certify patient identity.
cmt.rename('CMT').to_csv(r/'private/sample_CMT_PRIVATE.tsv',sep='\t')
proof={'status':'SAME_V4_COMMUNICATION_INPUTS_PREPARED_PENDING_INDEPENDENT_AUDIT',
    'UTC':datetime.now(timezone.utc).isoformat(),'seconds':time.time()-start,'n_cells':len(q),
    'n_samples':164,'n_sources':13,'n_original_CM':14,'n_labels':45,'n_common_measured_genes':len(common),
    'n_CPDB_measured_genes':lr.n_vars,'n_global_cells':glob.n_obs,'global_cap':1000,
    'normalization':'full native library total 1e4 then log1p; only genes measured in all retained 13 sources admitted',
    'source_or_label_or_CM_changed':False,'patient_independence_certified':False,
    'annotation_SHA256':'2a4bc82cb495d0fc2ffb40ecda6b1cbb53802eb2bf87efd2392746371653244c',
    'base_H_SHA256':hashlib.sha256((base/'CM_results/H_PRIVATE.tsv').read_bytes()).hexdigest(),
    'selected_rows_SHA256':hashlib.sha256(positions.tobytes()).hexdigest(),
    'global_input_SHA256':hashlib.sha256((r/'inputs/global_1000_PRIVATE.h5ad').read_bytes()).hexdigest()}
(r/'inputs/complete.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps(proof))
