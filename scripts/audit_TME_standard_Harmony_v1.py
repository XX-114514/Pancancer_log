"""Read back full model axes and bounded independent native-count PCA probes."""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def main(R):
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer'])
 while not (P/'provenance/supervisor_terminal_v1.json').exists():time.sleep(5)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;pc=load(P/'config/spec_v1.json');t=load(P/'provenance/terminal_v1.json');m=load(P/'outputs/standard_TME_Harmony_manifest_v1.json');assert sha(P/'outputs/standard_TME_Harmony_manifest_v1.json')==t['manifest_SHA256']
 for r in m['files']:
  p=P/'outputs'/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['SHA256']
 scope=load(P/'outputs/scope_and_common_genes_v1.json');selected=load(P/'outputs/selected_HVG_names_v1.json');tab=pd.read_csv(P/'outputs/standard_batch_HVG_v1.tsv',sep='\t',index_col=0);assert len(selected)==len(set(selected))==pc['n_HVGs'];assert tab.selected_HVG.sum()==pc['n_HVGs'];assert not tab.loc[selected,'technical_gene'].any()
 with np.load(P/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 raw=np.load(P/'outputs/X_pca_v1.npy',mmap_mode='r');corrected=np.load(P/'outputs/X_pca_harmony_v1.npy',mmap_mode='r');totals=np.load(P/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');assert raw.shape==corrected.shape==(scope['n_fit_records'],pc['n_PCs']);assert len(totals)==len(raw)
 for lo in range(0,len(raw),65536):assert np.isfinite(raw[lo:lo+65536]).all() and np.isfinite(corrected[lo:lo+65536]).all() and (totals[lo:lo+65536]>0).all()
 with np.load(P/'outputs/PCA_model_v1.npz',allow_pickle=False) as f:loadings=f['loadings'];mean=f['mean'];variance=f['variance'];vr=f['variance_ratio']
 np.testing.assert_allclose(loadings.T@loadings,np.eye(pc['n_PCs']),rtol=3e-4,atol=3e-4);assert (variance>=0).all() and np.isfinite(vr).all() and 0<vr.sum()<=1.0001
 J=Path(pc['joined']);jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(Path(cfg['independent_native_reader'])/'scripts'));from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(jm['archive_manifest']);a=load(pc['working_acceptance']);w=load(a['producer_manifest']['path']);W=Path(w['index']['path']).parent;idx={r['source']:r for r in load(w['index']['path'])};checks=0
 for code,source in enumerate(scope['sources']):
  wh=np.flatnonzero(axes['source_code']==code);wr=idx[source]
  with np.load(W/wr['arrays'],allow_pickle=False) as f:fit=np.flatnonzero(f['fit_pool_candidate_not_final_QC']);np.testing.assert_array_equal(axes['source_native_position'][wh],f['source_native_position'][fit]);np.testing.assert_array_equal(axes['native_row'][wh],f['native_row'][fit]);np.testing.assert_array_equal(axes['native_total_RNA'][wh],f['native_total_RNA'][fit])
  q=wh[np.unique(np.linspace(0,len(wh)-1,min(len(wh),cfg['checks_per_source'])).astype(int))];z=native.store(source);key=native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key'];x=gather(z[key],axes['native_row'][q]).astype(np.float64).tocsr();np.testing.assert_array_equal(np.asarray(x.sum(axis=1)).ravel(),axes['native_total_RNA'][q])
  if 'ambient_removed_counts' in z['layers']:x=(x-gather(z['layers/ambient_removed_counts'],axes['native_row'][q])).tocsr()
  x.eliminate_zeros();assert x.data.min(initial=0)>=-1e-8;RNA=np.asarray(x.sum(axis=1)).ravel();np.testing.assert_allclose(RNA,totals[q],rtol=1e-10,atol=1e-7);genes,cols=native_measured_genes(native,source);by=dict(zip(genes,cols));v=x[:,[int(by[g]) for g in selected]].toarray();v=np.log1p(v*(10000/RNA)[:,None]);expected=(v-mean)@loadings;np.testing.assert_allclose(expected,raw[q],rtol=3e-4,atol=3e-3);checks+=len(q);native._fields.clear()
 save(R/'outputs/standard_Harmony_readback_acceptance_v1.json',dict(n_sources=len(scope['sources']),n_fit_records=len(raw),n_independent_native_PCA_probes=checks,full_row_axes_and_all_coordinate_finiteness_checked=True,model_files_full_SHA_verified=True,independent_native_count_and_adopted_layer_PCA_probe_pass=True,QC_biology_and_integration_quality_not_proven_by_numeric_checks=True));save(R/'provenance/terminal_v1.json',dict(status='STANDARD_TME_HARMONY_ROW_AXES_AND_BOUNDED_NATIVE_PCA_READBACK_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),n_fit_records=len(raw)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
