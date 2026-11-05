"""Export a new V7-Mast-only raw-count object for lineage reanalysis.
No cell-level exports leave the server. No input/legacy annotation is modified.
"""
import os,json,hashlib,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,zarr,h5py,anndata as ad
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));assert ROOT.is_absolute()
RUN=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1'
OUT=RUN/'stage_b_mast_fullgene';OUT.mkdir(parents=True,exist_ok=False)
START=datetime.now(timezone.utc).isoformat()
B=ROOT/'runs/20260818_pancancer_5M_merge_prep'
df=pd.read_pickle(B/'11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz')
g=zarr.open_group(str(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r')
ix=pd.Index(read_elem(g['obs/_index']));assert ix.equals(df.index)
assert tuple(g['X'].attrs['shape'])==(len(df),95038)
rows=np.flatnonzero(df.major_celltype_v7.eq('Mast').to_numpy());assert len(rows)==45036
ptr=g['X/indptr'][:];starts=ptr[rows];ends=ptr[rows+1];lens=ends-starts;nnz=int(lens.sum());newptr=np.r_[0,np.cumsum(lens)].astype(np.int64)
data=np.empty(nnz,np.int32);indices=np.empty(nnz,np.int32)
# Coalesce requested raw CSR intervals. Read no full matrix or normalized layer.
block_start=0;blocks=[]
for i in range(1,len(rows)):
 if ends[i]-starts[block_start]>4000000 or starts[i]-ends[i-1]>226095:
  blocks.append((block_start,i));block_start=i
blocks.append((block_start,len(rows)))
for bi,(a,b) in enumerate(blocks):
 lo=int(starts[a]);hi=int(ends[b-1]);dd=g['X/data'][lo:hi];ii=g['X/indices'][lo:hi]
 for k in range(a,b):
  nlo=int(newptr[k]);nhi=int(newptr[k+1]);s=int(starts[k]-lo);e=int(ends[k]-lo)
  data[nlo:nhi]=dd[s:e];indices[nlo:nhi]=ii[s:e]
 if bi%50==0:
  (OUT/f'progress_{bi:05d}.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'blocks_done':bi+1,'blocks_total':len(blocks),'selected_cells':len(rows),'selected_nnz':nnz}))
assert (data>=0).all()
X=sp.csr_matrix((data,indices,newptr),shape=(len(rows),95038));X.check_format(full_check=True)
var=pd.DataFrame(index=pd.Index(read_elem(g['var/_index']),name='gene_symbol'))
for c in ['measured_in_n_gse','measured_in_gse_fraction']:var[c]=read_elem(g['var/'+c])
assert var.index.is_unique
obs=df.iloc[rows][['gse_id','sampleID','cancer_type','cohort','cluster_uid','major_celltype_v7','final_annotation_v7','annotation_confidence_v7','strict_major_eligible_v7','strict_subtype_eligible_v7','requires_recluster_for_subtype_v2','major_resolution_reason_v2','celltypist_major','celltypist_confidence']].copy()
for c in ['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_hb','doublet_method','doublet_score','doublet_status','decontx_contamination','decontx_status','ev_celltypist_label','ev_celltypist_conf','scatomic_full_label','scatomic_full_confidence']:
 if c in g['obs']:obs[c]=np.asarray(read_elem(g['obs/'+c]))[rows]
obs['global_position']=rows
# Use globally unique cell IDs as join keys; retain sample IDs privately.
assert obs.index.equals(ix[rows])
obj=ad.AnnData(X,obs=obs,var=var)
obj.uns['provenance']={'input_scope':'Frozen V7 major_celltype_v7 == Mast; candidate pool, not accepted identity','normalization':'X is source raw integer counts','later_core_annotation_labels_used':False,'axis_exact_with_frozen_V7':True,'UTC_start':START}
obj.write_h5ad(OUT/'mast_v7_fullgene_raw.h5ad',compression='gzip')
# 10X-compatible gene x cell CSC uses exactly the same buffers as cell x gene CSR.
with h5py.File(OUT/'mast_v7_fullgene_raw_10x.h5','x') as f:
 m=f.create_group('matrix');m.create_dataset('data',data=data,compression='gzip');m.create_dataset('indices',data=indices,compression='gzip');m.create_dataset('indptr',data=newptr,compression='gzip');m.create_dataset('shape',data=np.array([95038,len(rows)],dtype=np.int64))
 text_type=h5py.string_dtype('utf-8');m.create_dataset('barcodes',data=np.array(obs.index,dtype=object),dtype=text_type)
 ft=m.create_group('features')
 for k,v in [('id',var.index),('name',var.index),('feature_type',['Gene Expression']*95038),('genome',['GRCh38']*95038)]:ft.create_dataset(k,data=np.array(v,dtype=object),dtype=text_type)
obs.to_csv(OUT/'mast_v7_metadata_PRIVATE.tsv.gz',sep='\t',index_label='cell_id')
qc=[]
for label,s in obs.groupby('final_annotation_v7',observed=True):
 row={'label_v7':label,'n_cells':len(s),'n_gse':s.gse_id.nunique()}
 for c in ['n_genes_by_counts','total_counts','pct_counts_mt','pct_counts_hb','doublet_score','decontx_contamination']:
  if c in s:
   vals=pd.to_numeric(s[c],errors='coerce');row[c+'_n_available']=int(vals.notna().sum())
   for q in [.1,.5,.9]:row[c+'_q'+str(int(q*100))]=float(vals.quantile(q))
 for c in ['doublet_method','doublet_status','decontx_status']:
  if c in s:row[c+'_counts']=json.dumps({str(k):int(v) for k,v in s[c].value_counts(dropna=False).items()})
 qc.append(row)
pd.DataFrame(qc).to_csv(OUT/'mast_v7_QC_aggregate.tsv',sep='\t',index=False)
check=ad.read_h5ad(OUT/'mast_v7_fullgene_raw.h5ad',backed='r');assert check.shape==X.shape and check.obs_names.equals(obs.index);check.file.close()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024**2),b''):h.update(v)
 return h.hexdigest()
files=[{'name':p.name,'size':p.stat().st_size,'sha256':sha(p),'private_cell_or_matrix':p.suffix in ['.h5ad','.h5','.gz']} for p in OUT.iterdir() if p.is_file() and not p.name.startswith('progress')]
manifest={'stage':'B_Mast_fullgene_export_complete_not_annotation_complete','UTC_start':START,'UTC_end':datetime.now(timezone.utc).isoformat(),'shape':list(X.shape),'nnz':nnz,'n_read_blocks':len(blocks),'input_axis_exact':True,'n_annotation_changes':0,'scope':'Only V7 Mast; no recruitment outside this parent yet','outputs':files}
(OUT/'export_complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest))
