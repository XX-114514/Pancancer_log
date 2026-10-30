"""All-cell native cluster profiles and official all-common-gene Wilcoxon.

No embedding is used as expression. Counts/QC/identity labels stay unchanged.
Cell-level marker screening is not patient-level differential inference.
"""
import argparse,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import psutil
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def check(R,stage,**kw):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw))
def csr32(x):
 x=x.tocsr().astype(np.float32);x.sum_duplicates();x.eliminate_zeros()
 if x.indices.dtype!=x.indptr.dtype:x.indices=x.indices.astype(np.int64);x.indptr=x.indptr.astype(np.int64)
 return x
def main(R):
 started=time.monotonic();cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);G=Path(cfg['graph']);A=Path(cfg['graph_auditor'])
 for D in [P,G,A]:assert load(D/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 audit=load(A/'outputs/accepted_same_graph_delivery_v1.json');assert audit['graph_terminal_SHA256']==sha(G/'provenance/terminal_v1.json');graph_files={r['name']:r for r in audit['all_output_file_SHA256']};assert sha(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz')==graph_files['private_full_graph_cluster_axes_v1.tsv.gz']['SHA256']
 import scanpy as sc,anndata as ad,numba
 from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
 sc.settings.n_jobs=cfg['threads'];numba.set_num_threads(cfg['threads']);pc=load(P/'config/spec_v1.json');J=Path(pc['joined']);sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);assert psutil.virtual_memory().available>=cfg['minimum_available_RAM_bytes']
 # Actual official implementation sanity check with known opposite markers.
 toy=ad.AnnData(sp.csr_matrix(np.vstack([np.tile([5,1,2,2],(30,1)),np.tile([1,5,2,2],(30,1))]).astype(np.float32)));toy.var_names=['A','B','C','D'];toy.obs['cluster']=pd.Categorical(np.repeat(['0','1'],30));sc.pp.normalize_total(toy,target_sum=10000);sc.pp.log1p(toy);sc.tl.rank_genes_groups(toy,'cluster',method='wilcoxon',tie_correct=True,use_raw=False,reference='rest',n_genes=4,pts=True,corr_method='benjamini-hochberg');assert toy.uns['rank_genes_groups']['names']['0'][0]=='A' and toy.uns['rank_genes_groups']['names']['1'][0]=='B';assert toy.uns['rank_genes_groups']['scores']['0'][0]>0 and toy.uns['rank_genes_groups']['pvals_adj']['0'][0]<.05;save(R/'provenance/official_Wilcoxon_API_smoke_v1.json',dict(actual_tie_corrected_Wilcoxon_opposite_marker_fixture_passed=True,synthetic_not_biological_validation=True,Scanpy=sc.__version__));del toy
 scope=load(P/'outputs/scope_and_common_genes_v1.json');genes=scope['genes'];N=scope['n_fit_records'];lineage=load(P/'outputs/source_input_lineage_v1.json');assert [r['source'] for r in lineage]==scope['sources']
 with np.load(P/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 totals_saved=np.load(P/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');q=pd.read_csv(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz',sep='\t',usecols=['leiden_0_5']);labels=q.leiden_0_5.to_numpy(np.int32);del q;ng=audit['clusters']['leiden_0_5'];assert len(labels)==N and np.array_equal(np.unique(labels),np.arange(ng));native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
 assert sha(pc['working_acceptance'])==pc['working_acceptance_SHA256'];wa=load(pc['working_acceptance']);wm=load(wa['producer_manifest']['path']);W=Path(wm['index']['path']).parent;idx={r['source']:r for r in load(wm['index']['path'])};crosswalk=[];chunks=[];records=[]
 groups=[]
 for k in range(ng):
  h=hashlib.sha256();n_sources=0
  for code,source in enumerate(scope['sources']):
   members=axes['native_row'][(labels==k)&(axes['source_code']==code)]
   if not len(members):continue
   b=source.encode();h.update(len(b).to_bytes(8,'little'));h.update(b);h.update(len(members).to_bytes(8,'little'));h.update(np.asarray(members,dtype='<i8').tobytes());n_sources+=1
  groups.append(dict(cluster=str(k),n_cells=int((labels==k).sum()),n_sources=n_sources,exact_membership_SHA256=h.hexdigest()))
 save(R/'outputs/cluster_membership_catalog_v1.json',groups)
 save(R/'outputs/marker_input_contract_v1.json',dict(n_fit=N,n_clusters=ng,common_gene_names=genes,n_common_genes=len(genes),all_current_fit_cells_in_Wilcoxon=True,cell_subsampling=False,rank_method='official_scanpy_wilcoxon_tie_correct_true_one_vs_rest_BH_per_cluster',layer='full_adopted_analysis_RNA_CP10K_log1p',not_HVG_restricted=True,native_noncommon_genes_retained_in_source_profiles=True,cell_level_screen_not_patient_inference=True,graph_readback_acceptance_SHA256=sha(A/'outputs/accepted_same_graph_delivery_v1.json')))
 for code,r in enumerate(lineage):
  source=r['source'];wh=np.flatnonzero(axes['source_code']==code);assert len(wh)==r['n_fit'];wr=idx[source];meta=load(W/wr['metadata']);assert sha(W/wr['metadata'])==wr['metadata_SHA256'];assert sha(W/wr['arrays'])==wr['arrays_SHA256']
  with np.load(W/wr['arrays'],allow_pickle=False) as f:
   fit=np.flatnonzero(f['fit_pool_candidate_not_final_QC']);np.testing.assert_array_equal(f['native_row'][fit],axes['native_row'][wh]);oldmajor=np.asarray(meta['candidate_dictionaries']['consensus_major'])[f['consensus_major_code'][fit]];oldfine=np.asarray(meta['candidate_dictionaries']['consensus_celltype'])[f['consensus_celltype_code'][fit]];oldgroup=f['annotation_group_code'][fit]
  oldview=pd.DataFrame(dict(cluster=labels[wh],previous_major=oldmajor,previous_celltype=oldfine,previous_source_group=oldgroup));tab=oldview.groupby(list(oldview.columns),observed=True).size().reset_index(name='n_cells');tab.insert(0,'source',source);crosswalk.extend(tab.to_dict('records'));del oldview,oldmajor,oldfine,oldgroup
  z=native.store(source);names,physical=native_measured_genes(native,source);name_index={g:i for i,g in enumerate(names)};common_local=np.asarray([name_index[g] for g in genes]);n=np.bincount(labels[wh],minlength=ng);stats={k:np.zeros((ng,len(names)),np.float64) for k in ['CP10K_sum','logCP10K_sum','detected']};delta=z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None;nnz=0
  for lo in range(0,len(wh),cfg['read_block_size']):
   take=wh[lo:lo+cfg['read_block_size']];raw=gather_csr_rows(z[r['counts_key']],axes['native_row'][take],max_span=8192).astype(np.float64).tocsr();np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),axes['native_total_RNA'][take]);assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();x=raw
   if delta is not None:x=(raw-gather_nonnegative_float_delta(delta,axes['native_row'][take])).tocsr()
   x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();RNA=np.asarray(x.sum(axis=1)).ravel();np.testing.assert_allclose(RNA,totals_saved[take],rtol=1e-10,atol=1e-7);assert (RNA>0).all();native_x=x[:,physical].tocsr();nnz+=native_x.nnz;cp=native_x.multiply((10000/RNA)[:,None]).tocsr();logged=cp.copy();logged.data=np.log1p(logged.data);link=sp.csr_matrix((np.ones(len(take)),(labels[take],np.arange(len(take)))),shape=(ng,len(take)))
   for key,values in [('CP10K_sum',cp),('logCP10K_sum',logged),('detected',native_x>0)]:stats[key]+=(link@values).toarray()
   chunks.append(csr32(logged[:,common_local]));del raw,x,native_x,cp,logged,link
  out=R/'outputs'/('native_cluster_profiles_'+source+'.npz');np.savez_compressed(out,cluster_numbers=np.arange(ng,dtype=np.int32),n_cells=n,native_gene_names=np.asarray(names,dtype=str),native_columns=physical,**stats);record=dict(source=source,n_cells=len(wh),n_native_unique_measured_genes=len(names),n_unique_native_nnz_read=nnz,path=out.name,SHA256=sha(out),zero_size_source_cluster_is_unobserved=True,source_input_lineage=r);records.append(record);save(R/'outputs'/('native_cluster_profiles_'+source+'.json'),record);del stats;native._fields.clear();gc.collect();storage_gate(R,cfg);check(R,'full_native_cluster_source_profiles_complete',source=source,n_cells=len(wh),n_genes=len(names))
 pd.DataFrame(crosswalk).to_csv(R/'outputs/previous_working_labels_to_new_clusters_v1.tsv',sep='\t',index=False);save(R/'outputs/native_profiles_index_v1.json',records);X=csr32(sp.vstack(chunks,format='csr'));del chunks;gc.collect();assert X.shape==(N,len(genes));np.save(R/'outputs/private_Wilcoxon_cluster_axis_v1.npy',labels);a=ad.AnnData(X,obs=pd.DataFrame({'cluster':pd.Categorical.from_codes(labels,[str(i) for i in range(ng)])}),var=pd.DataFrame(index=genes));a.uns['log1p']={'base':None};del X;check(R,'official_all_cell_common_gene_Wilcoxon_started',n_cells=N,n_genes=len(genes),n_clusters=ng)
 sc.tl.rank_genes_groups(a,'cluster',method='wilcoxon',tie_correct=True,use_raw=False,reference='rest',n_genes=len(genes),pts=True,corr_method='benjamini-hochberg');df=sc.get.rank_genes_groups_df(a,group=None);assert len(df)==ng*len(genes);assert np.isfinite(df['scores']).all() and df['pvals'].between(0,1).all() and df['pvals_adj'].between(0,1).all();df.to_csv(R/'outputs/all_cell_common_gene_Wilcoxon_v1.tsv.gz',sep='\t',index=False);top=df.groupby('group',observed=True,sort=False).head(30);top.to_csv(R/'outputs/all_cell_Wilcoxon_top30_per_cluster_v1.tsv',sep='\t',index=False);check(R,'official_all_cell_Wilcoxon_complete',n_marker_records=len(df));del a;gc.collect();storage_gate(R,cfg)
 save(R/'provenance/terminal_v1.json',dict(status='FULL_NATIVE_CLUSTER_PROFILES_AND_OFFICIAL_ALL_CELL_COMMON_GENE_WILCOXON_COMPLETE_REQUIRES_INDEPENDENT_AUDIT_AND_MARKER_ADJUDICATION',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,n_clusters=ng,n_common_genes=len(genes),cell_subsampling=False,original_counts_QC_and_labels_unchanged=True,no_new_final_celltype_assigned=True,wall_seconds=time.monotonic()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
