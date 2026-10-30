"""Independent whole-axis graph/Leiden/UMAP delivery readback; no new labels."""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def main(R):
 cfg=load(R/'config/spec_v1.json');G=Path(cfg['graph']);P=Path(cfg['producer']);A=Path(cfg['producer_auditor'])
 for D in [G,P,A]:assert load(D/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 t=load(G/'provenance/terminal_v1.json');scope=load(P/'outputs/scope_and_common_genes_v1.json');N=scope['n_fit_records'];assert t['n_fit_records']==N
 pc=load(P/'config/spec_v1.json');sys.path.insert(0,str(Path(pc['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg)
 with np.load(P/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 q=pd.read_csv(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz',sep='\t',index_col=0);assert len(q)==N and not q.isna().any().any();np.testing.assert_array_equal(np.asarray(q.index,dtype=np.int64),np.arange(N));np.testing.assert_array_equal(q['source_native_position'],axes['source_native_position']);np.testing.assert_array_equal(q['native_row'],axes['native_row']);np.testing.assert_array_equal(q['source'],np.asarray(scope['sources'])[axes['source_code']])
 xy=np.load(G/'outputs/full_same_graph_UMAP_v1.npy',mmap_mode='r');assert xy.shape==(N,2)
 for lo in range(0,N,65536):assert np.isfinite(xy[lo:lo+65536]).all()
 counts={}
 for k in ['leiden_0_5','leiden_1_0','leiden_2_0']:
  assert k in q and (pd.to_numeric(q[k])>=0).all();v=q[k].value_counts();assert v.sum()==N;counts[k]=int(len(v))
 assert counts['leiden_0_5']==t['n_primary_clusters'];saved=pd.read_csv(G/'outputs/numeric_cluster_summary_v1.tsv',sep='\t',index_col=0);expected=q.groupby('leiden_0_5',observed=True).agg(n_cells=('source','size'),n_sources=('source','nunique'));saved.index=saved.index.astype(str);expected.index=expected.index.astype(str);pd.testing.assert_frame_equal(saved.sort_index(),expected.sort_index(),check_index_type=False)
 graph=sp.load_npz(G/'outputs/full_graph_connectivities_v1.npz');assert graph.shape==(N,N);assert np.isfinite(graph.data).all() and graph.data.min(initial=0)>=0 and graph.data.max(initial=0)<=1.00001;diff=(graph-graph.T).tocsr();assert np.abs(diff.data).max(initial=0)<1e-6;zero_rows=int((graph.getnnz(axis=1)==0).sum());del diff,graph
 distances=sp.load_npz(G/'outputs/full_graph_distances_v1.npz');assert distances.shape==(N,N) and np.isfinite(distances.data).all() and distances.data.min(initial=0)>=0;del distances
 outputs=[dict(name=p.name,bytes=p.stat().st_size,SHA256=sha(p)) for p in (G/'outputs').iterdir() if p.is_file()]
 export_names=['current_TME_same_graph_cluster_source_UMAP_v1.png','numeric_cluster_summary_v1.tsv','UMAP_actual_parameters_v1.json'];export=[r for r in outputs if r['name'] in export_names];assert len(export)==3
 save(R/'outputs/accepted_same_graph_delivery_v1.json',dict(n_fit=N,n_sources=len(scope['sources']),clusters=counts,all_row_axes_and_coordinates_checked=True,full_graph_numeric_symmetry_checked=True,graph_zero_degree_rows=zero_rows,primary_cluster_counts_and_source_coverage_independently_reproduced=True,all_output_file_SHA256=outputs,graph_terminal_SHA256=sha(G/'provenance/terminal_v1.json'),safe_static_group_exports=export,no_final_celltype_labels_assigned=True,biological_quality_not_proven_by_numeric_readback=True))
 save(R/'provenance/terminal_v1.json',dict(status='WHOLE_AXIS_SAME_GRAPH_LEIDEN_UMAP_DELIVERY_NUMERICALLY_ACCEPTED_REQUIRES_MARKER_AND_BIOLOGICAL_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(R/'outputs/accepted_same_graph_delivery_v1.json'),n_fit=N))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
