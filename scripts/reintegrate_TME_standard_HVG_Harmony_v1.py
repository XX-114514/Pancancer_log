"""Current-scope standard Seurat-v3 HVG/PCA/full Harmony, append-only.

Scientific recipe follows the existing HVG4000 pipeline and installed Scanpy.
Only row/feature/layer assembly is custom. No annotation or new QC is inferred.
"""
import argparse,collections,datetime,gc,hashlib,importlib.metadata,inspect,json,os,re,resource,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import psutil

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def checkpoint(R,stage,**kw):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw))
def csr32(x):
 x=x.tocsr().astype(np.float32);x.sum_duplicates();x.eliminate_zeros()
 if x.indices.dtype!=x.indptr.dtype:x.indices=x.indices.astype(np.int64);x.indptr=x.indptr.astype(np.int64)
 return x
def technical(g):
 return g.upper()=='MALAT1' or bool(re.search(r'^(MT-|RP[SL][0-9]|HB[ABDEGQMZ][0-9]|ERCC-)',g,re.I))

def main(R):
 started=time.monotonic();cfg=load(R/'config/spec_v1.json');J=Path(cfg['joined']);sys.path.insert(0,str(J/'scripts'))
 import anndata as ad,scanpy as sc,harmonypy as hm,numba
 from skmisc.loess import loess
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 sc.settings.n_jobs=cfg['threads'];numba.set_num_threads(cfg['threads']);storage_gate(R,cfg)
 ac=Path(cfg['working_acceptance']);assert sha(ac)==cfg['working_acceptance_SHA256'];a=load(ac);assert sha(a['producer_manifest']['path'])==a['producer_manifest']['SHA256'];w=load(a['producer_manifest']['path']);W=Path(w['index']['path']).parent;assert sha(w['index']['path'])==w['index']['SHA256'];wi=load(w['index']['path'])
 jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');native=NativeArchiveCollection(jm['archive_manifest']);sources=sorted(r['source'] for r in wi);lookup={r['source']:r for r in wi};common=None;metadata={};rows_total=0
 runtime={k:importlib.metadata.version(k) for k in ['scanpy','anndata','harmonypy','scikit-misc','numpy','scipy']};runtime['loess_extension']=inspect.getfile(loess);runtime['loess_extension_SHA256']=sha(runtime['loess_extension']);save(R/'provenance/runtime_v1.json',runtime)
 # Exercise the actual compiled LOESS-backed standard API before native reads.
 rng=np.random.default_rng(cfg['seed']);test=ad.AnnData(sp.csr_matrix(rng.poisson(np.geomspace(.1,10,500),(240,500)).astype(np.float32)));test.obs['source']=pd.Categorical(np.repeat(['a','b','c'],80));sc.pp.highly_variable_genes(test,flavor='seurat_v3',batch_key='source',n_top_genes=100,check_values=True);assert int(test.var.highly_variable.sum())==100 and np.isfinite(test.var.variances_norm).all();save(R/'provenance/standard_HVG_API_smoke_v1.json',dict(actual_compiled_LOESS_and_Scanpy_batch_API_passed=True,synthetic_not_biological_validation=True));del test
 for source in sources:
  wr=lookup[source];assert sha(W/wr['metadata'])==wr['metadata_SHA256'] and sha(W/wr['arrays'])==wr['arrays_SHA256']
  with np.load(W/wr['arrays'],allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','fit_pool_candidate_not_final_QC','source_qualification_code','native_total_RNA']}
  assert (q['source_qualification_code']==1).all();fit=np.flatnonzero(q['fit_pool_candidate_not_final_QC']);assert len(fit)>100;np.testing.assert_array_equal(native.rows(source)[q['source_native_position'][fit]],q['native_row'][fit]);genes,cols=native_measured_genes(native,source);common=set(genes) if common is None else common&set(genes)
  metadata[source]=dict(q=q,fit=fit,columns=dict(zip(genes,map(int,cols))),input_record=wr);rows_total+=len(fit)
  checkpoint(R,'scope_and_measured_features_verified',source=source,n_fit=len(fit));native._fields.clear()
 genes=sorted(common);assert len(genes)==cfg['expected_common_features'] and rows_total==cfg['expected_fit_records'] and len(sources)==cfg['expected_sources'];p=len(genes);assert p>=cfg['candidate_HVGs'];assert psutil.virtual_memory().available>cfg['required_available_RAM_bytes']
 save(R/'outputs/scope_and_common_genes_v1.json',dict(n_sources=len(sources),sources=sources,n_fit_records=rows_total,n_qualified_records=sum(len(v['q']['native_row']) for v in metadata.values()),genes=genes,all_features_uniquely_measured_in_every_source=True,missing_not_zero_imputed=True,working_acceptance_SHA256=sha(ac)))
 chunks=[];source_axis=[];position_axis=[];native_axis=[];RNA_axis=[];offset=0;scopes=[]
 for code,source in enumerate(sources):
  v=metadata[source];q=v['q'];fit=v['fit'];z=native.store(source);key=native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key'];csr=z[key];cols=np.asarray([v['columns'][g] for g in genes]);localRNA=[]
  for lo in range(0,len(fit),cfg['read_block_size']):
   take=fit[lo:lo+cfg['read_block_size']];x=gather_csr_rows(csr,q['native_row'][take],max_span=8192).tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all();tot=np.asarray(x.sum(axis=1,dtype=np.float64)).ravel();np.testing.assert_array_equal(tot,q['native_total_RNA'][take]);assert x.data.max(initial=0)<2**24;chunks.append(csr32(x[:,cols]));localRNA.append(tot)
  source_axis.append(np.full(len(fit),code,np.int16));position_axis.append(q['source_native_position'][fit]);native_axis.append(q['native_row'][fit]);RNA_axis.append(np.concatenate(localRNA));scopes.append(dict(source=source,offset=offset,n_fit=len(fit),n_qualified=len(q['native_row']),counts_key=key,store=str(z.store.path),has_adopted_ambient_delta='ambient_removed_counts' in z['layers'],working_arrays_SHA256=sha(W/v['input_record']['arrays'])));offset+=len(fit);checkpoint(R,'native_integer_HVG_input_source_complete',source=source,n_records=len(fit));native._fields.clear()
 X=csr32(sp.vstack(chunks,format='csr'));del chunks;gc.collect();source_code=np.concatenate(source_axis);positions=np.concatenate(position_axis);native_rows=np.concatenate(native_axis);native_RNA=np.concatenate(RNA_axis);del source_axis,position_axis,native_axis,RNA_axis
 axes=R/'outputs/private_cell_axes_v1.npz';np.savez_compressed(axes,source_code=source_code,source_native_position=positions,native_row=native_rows,native_total_RNA=native_RNA);save(R/'outputs/source_input_lineage_v1.json',scopes)
 A=ad.AnnData(X,obs=pd.DataFrame({'source':pd.Categorical.from_codes(source_code,sources)}),var=pd.DataFrame(index=pd.Index(genes)));del X;gc.collect();checkpoint(R,'official_seurat_v3_HVG_started',n_records=rows_total,n_common_features=p)
 sc.pp.highly_variable_genes(A,n_top_genes=cfg['candidate_HVGs'],flavor='seurat_v3',batch_key='source',inplace=True,check_values=True)
 table=A.var.copy();table['technical_gene']=[technical(g) for g in genes];eligible=table.loc[~table.technical_gene&table.highly_variable_rank.notna()].sort_values(['highly_variable_rank','highly_variable_nbatches'],ascending=[True,False],kind='mergesort');assert len(eligible)>=cfg['n_HVGs'];selected=list(eligible.index[:cfg['n_HVGs']]);table['selected_HVG']=table.index.isin(selected);table.index.name='gene';table.to_csv(R/'outputs/standard_batch_HVG_v1.tsv',sep='\t');save(R/'outputs/selected_HVG_names_v1.json',selected);del A;gc.collect();checkpoint(R,'official_seurat_v3_HVG_complete',n_selected=len(selected));storage_gate(R,cfg)
 # Normalize the full already-adopted analysis counts BEFORE feature subsetting.
 chunks=[];analysisRNA=[];layer_reports=[]
 for source in sources:
  v=metadata[source];q=v['q'];fit=v['fit'];z=native.store(source);key=native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key'];cols=np.asarray([v['columns'][g] for g in selected]);delta=z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None;n_changed=0;removed_total=0.
  for lo in range(0,len(fit),cfg['read_block_size']):
   take=fit[lo:lo+cfg['read_block_size']];raw=gather_csr_rows(z[key],q['native_row'][take],max_span=8192).astype(np.float64).tocsr();np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),q['native_total_RNA'][take]);x=raw
   if delta is not None:
    removed=gather_csr_rows(delta,q['native_row'][take],max_span=8192).astype(np.float64).tocsr();assert np.isfinite(removed.data).all() and (removed.data>=0).all();x=(raw-removed).tocsr();x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();loss=np.asarray(removed.sum(axis=1)).ravel();n_changed+=int((loss>0).sum());removed_total+=float(loss.sum())
   totals=np.asarray(x.sum(axis=1)).ravel();assert (totals>0).all();analysisRNA.append(totals);block=ad.AnnData(csr32(x));sc.pp.normalize_total(block,target_sum=10000);sc.pp.log1p(block);y=csr32(block.X[:,cols]);expected=np.log1p(np.asarray(x[:,cols[:3]].toarray())*(10000/totals)[:,None]);np.testing.assert_allclose(y[:,:3].toarray(),expected,rtol=3e-6,atol=2e-6);chunks.append(y)
  layer_reports.append(dict(source=source,recipe='published_counts_minus_already_adopted_ambient_delta' if delta is not None else 'published_legacy_integer_counts',n_current_fit_records=len(fit),n_fit_records_with_nonzero_previously_adopted_delta=n_changed,adopted_removed_molecules=removed_total,new_QC_or_ambient_refit=False));checkpoint(R,'adopted_analysis_layer_normalized_source_complete',source=source,n_fit=len(fit),n_previously_corrected=n_changed);native._fields.clear()
 Y=csr32(sp.vstack(chunks,format='csr'));del chunks;gc.collect();np.save(R/'outputs/private_analysis_RNA_totals_v1.npy',np.concatenate(analysisRNA));save(R/'outputs/adopted_analysis_layer_lineage_v1.json',layer_reports)
 B=ad.AnnData(Y,obs=pd.DataFrame({'source':pd.Categorical.from_codes(source_code,sources)}),var=pd.DataFrame(index=selected));del Y;gc.collect();checkpoint(R,'standard_centered_ARPACK_PCA_started',n_records=rows_total,n_features=len(selected));sc.pp.pca(B,n_comps=cfg['n_PCs'],zero_center=True,svd_solver='arpack',random_state=cfg['seed'],dtype='float32',mask_var=None)
 raw=np.asarray(B.obsm['X_pca'],np.float32,order='C');assert raw.shape==(rows_total,cfg['n_PCs']) and np.isfinite(raw).all();np.save(R/'outputs/X_pca_v1.npy',raw);np.savez_compressed(R/'outputs/PCA_model_v1.npz',loadings=B.varm['PCs'],mean=np.asarray(B.X.mean(axis=0)).ravel(),variance=B.uns['pca']['variance'],variance_ratio=B.uns['pca']['variance_ratio']);pca_summary=dict(n_records=rows_total,n_HVGs=len(selected),n_PCs=cfg['n_PCs'],zero_center=True,per_gene_variance_scaling=False,normalization_on_full_adopted_RNA_before_HVG_subset=True);save(R/'outputs/PCA_summary_v1.json',pca_summary);del B;gc.collect();storage_gate(R,cfg);checkpoint(R,'full_current_TME_Harmony_started',n_records=rows_total,n_batches=len(sources))
 harmony=hm.run_harmony(data_mat=raw,meta_data=pd.DataFrame({'source':pd.Categorical.from_codes(source_code,sources)}),vars_use=['source'],theta=cfg['theta'],lamb=cfg['lambda'],max_iter_harmony=cfg['max_iter_harmony'],random_state=cfg['seed'],ncores=cfg['threads'],verbose=True);z=np.asarray(harmony.Z_corr)
 if z.shape==raw.T.shape:z=z.T
 assert z.shape==raw.shape and np.isfinite(z).all();z=np.asarray(z,np.float32,order='C');np.save(R/'outputs/X_pca_harmony_v1.npy',z);history={k:np.asarray(getattr(harmony,k)).tolist() for k in ['objective_harmony','objective_kmeans','kmeans_rounds'] if hasattr(harmony,k)};conv=getattr(harmony,'converged',None);history.update(converged=None if conv is None else bool(conv));save(R/'outputs/Harmony_optimization_history_v1.json',history)
 files=[p for p in (R/'outputs').iterdir() if p.is_file()];manifest=dict(files=[dict(path=p.name,bytes=p.stat().st_size,SHA256=sha(p)) for p in files],n_fit_records=rows_total,n_sources=len(sources),count_matrices_copied_to_disk=False,working_labels_and_previous_models_unchanged=True,features='official_Scanpy_seurat_v3_batch_source_no_custom_dispersion',new_clustering_annotation_and_UMAP_pending=True);save(R/'outputs/standard_TME_Harmony_manifest_v1.json',manifest);save(R/'provenance/terminal_v1.json',dict(status='CURRENT_TME_STANDARD_SEURAT_V3_PCA_FULL_HARMONY_COMPLETE_REQUIRES_READBACK_AND_SAME_GRAPH_CLUSTER_UMAP',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),manifest_SHA256=sha(R/'outputs/standard_TME_Harmony_manifest_v1.json'),n_fit_records=rows_total,n_sources=len(sources),wall_seconds=time.monotonic()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,previous_counts_labels_and_QC_unchanged=True))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
