"""Full current-TME same-graph Leiden/UMAP; labels remain numeric candidates."""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def check(R,stage,**v):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**v))
def main(R):
 started=time.monotonic();cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor'])
 for D in [P,A]:
  while not (D/'provenance/supervisor_terminal_v1.json').exists():time.sleep(5)
  assert load(D/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json')
 import anndata as ad,scanpy as sc,numba
 sc.settings.n_jobs=cfg['threads'];numba.set_num_threads(cfg['threads']);J=Path(load(P/'config/spec_v1.json')['joined']);sys.path.insert(0,str(J/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);import psutil;assert psutil.virtual_memory().available>=cfg['minimum_available_RAM_bytes'];scope=load(P/'outputs/scope_and_common_genes_v1.json');z=np.load(P/'outputs/X_pca_v1.npy');assert np.isfinite(z).all()
 with np.load(P/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:source=f['source_code'];position=f['source_native_position'];native_row=f['native_row']
 a=ad.AnnData(sp.csr_matrix((len(z),0)),obs=pd.DataFrame({'source':pd.Categorical.from_codes(source,scope['sources']),'source_native_position':position,'native_row':native_row}));a.obsm['X_pca']=z;check(R,'full_raw_PCA_knn_started',n_records=len(z),n_neighbors=cfg['n_neighbors'])
 sc.pp.neighbors(a,use_rep='X_pca',n_neighbors=cfg['n_neighbors'],random_state=cfg['seed'],transformer='pynndescent',metric='euclidean');check(R,'full_raw_PCA_knn_complete',n_connectivity_edges=a.obsp['connectivities'].nnz);sp.save_npz(R/'outputs/full_graph_connectivities_v1.npz',a.obsp['connectivities'],compressed=False);sp.save_npz(R/'outputs/full_graph_distances_v1.npz',a.obsp['distances'],compressed=False)
 for resolution in cfg['resolutions']:
  key='leiden_'+str(resolution).replace('.','_');sc.tl.leiden(a,resolution=resolution,random_state=cfg['seed'],flavor='igraph',n_iterations=2,directed=False,key_added=key);check(R,'numeric_Leiden_complete',resolution=resolution,n_clusters=int(a.obs[key].nunique()))
 a.obs.to_csv(R/'outputs/private_full_graph_cluster_axes_v1.tsv.gz',sep='\t');check(R,'full_same_graph_UMAP_started',n_records=a.n_obs);sc.tl.umap(a,random_state=cfg['seed']);xy=np.asarray(a.obsm['X_umap'],np.float32);assert xy.shape==(len(z),2) and np.isfinite(xy).all();np.save(R/'outputs/full_same_graph_UMAP_v1.npy',xy);save(R/'outputs/UMAP_actual_parameters_v1.json',dict(umap=a.uns['umap'],neighbors=a.uns['neighbors'],same_connectivities_consumed_by_Leiden_and_UMAP=True))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 rng=np.random.default_rng(cfg['seed']);take=[]
 for code in range(len(scope['sources'])):
  q=np.flatnonzero(source==code);take.append(rng.choice(q,min(len(q),cfg['plot_per_source']),replace=False))
 primary='leiden_'+str(cfg['primary_resolution']).replace('.','_');labels=a.obs[primary].to_numpy().astype(str)
 for k in set(labels):
  q=np.flatnonzero(labels==k);take.append(rng.choice(q,min(len(q),5),replace=False))
 chosen=np.unique(np.concatenate(take));fig,axs=plt.subplots(1,2,figsize=(19,8));cmap=plt.get_cmap('tab20')
 for i,k in enumerate(sorted(set(labels),key=int)):
  q=chosen[labels[chosen]==k];axs[0].scatter(xy[q,0],xy[q,1],s=.8,c=[cmap(i%20)],linewidths=0,rasterized=True);center=np.median(xy[np.flatnonzero(labels==k)],axis=0);axs[0].text(center[0],center[1],k,fontsize=7,ha='center',bbox=dict(facecolor='white',alpha=.7,pad=.2,edgecolor='none'))
 axs[1].scatter(xy[chosen,0],xy[chosen,1],s=.8,c=source[chosen],cmap='turbo',linewidths=0,rasterized=True)
 for ax in axs:ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
 axs[0].set_title('Numeric Leiden clusters | full raw PCA graph');axs[1].set_title('Source context on identical full-graph coordinates');fig.suptitle('Current TME raw PCA baseline; same graph for clustering and UMAP; source-balanced plotted subset');fig.tight_layout();fig.savefig(R/'outputs/current_TME_same_graph_cluster_source_UMAP_v1.png',dpi=180,bbox_inches='tight');plt.close(fig)
 summary=a.obs.groupby(primary,observed=True).agg(n_cells=('source','size'),n_sources=('source','nunique'));summary.to_csv(R/'outputs/numeric_cluster_summary_v1.tsv',sep='\t');save(R/'provenance/terminal_v1.json',dict(status='FULL_CURRENT_TME_RAW_PCA_SAME_GRAPH_LEIDEN_UMAP_COMPLETE_REQUIRES_BIOLOGICAL_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit_records=len(z),n_primary_clusters=len(summary),primary_resolution=cfg['primary_resolution'],full_graph_UMAP_not_reference_overlay=True,n_plot_points=len(chosen),new_final_celltype_labels_assigned=False,wall_seconds=time.monotonic()-started))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
