"""Independent model/latent/exposure audit; CPU inference of GPU-trained model."""
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
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer'])
 while not (P/'provenance/supervisor_terminal_v1.json').exists():time.sleep(5)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 pc=load(P/'config/spec_v1.json');terminal=load(P/'provenance/terminal_v1.json');manifest=load(P/'outputs/scVI_GPU_manifest_v1.json');assert sha(P/'outputs/scVI_GPU_manifest_v1.json')==terminal['manifest_SHA256']
 for r in manifest['files']:
  p=P/'outputs'/r['path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['SHA256']
 inp=load(P/'outputs/private_native_input_spec_v1.json');assert sha(inp['axes_path'])==inp['axes_SHA256']==manifest['axes_SHA256'];H=Path(pc['producer']);assert sha(H/'outputs/selected_HVG_names_v1.json')==inp['genes_SHA256']==manifest['genes_SHA256'];N=inp['n_fit']
 with np.load(inp['axes_path'],allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 latent=np.load(P/'outputs/X_scVI_v1.npy',mmap_mode='r');assert latent.shape==(N,pc['n_latent'])
 for lo in range(0,N,65536):assert np.isfinite(latent[lo:lo+65536]).all()
 epochs=[]
 for p in sorted((P/'checkpoints').glob('*.json')):
  r=load(p)
  if r['stage']=='official_GPU_scVI_epoch_complete':epochs.append(r)
 assert [r['epoch'] for r in epochs]==list(range(1,pc['max_epochs']+1));assert all(r['device'].startswith('cuda:') for r in epochs);assert all(np.isfinite(list(r['metrics'].values())).all() for r in epochs);assert epochs[-1]['kl_weight']>=.999;assert sum(r['kl_weight']>=.999 for r in epochs)>=pc['max_epochs']-pc['n_epochs_kl_warmup']
 for name in ['elbo_train','elbo_validation','train_loss','validation_loss']:
  q=pd.read_csv(P/'outputs'/('history_'+name+'.tsv'),sep='\t',index_col=0);assert len(q)==pc['max_epochs'] and np.isfinite(q.to_numpy(dtype=float)).all()
 probes=[]
 for code,r in enumerate(inp['records']):
  wh=np.flatnonzero(axes['source_code']==code);assert len(wh)==r['n_fit'];probes.extend(wh[np.unique(np.linspace(0,len(wh)-1,min(cfg['probes_per_source'],len(wh))).astype(int))])
 probes=np.asarray(probes,np.int64);X=np.zeros((len(probes),len(inp['genes'])),np.float32)
 import zarr,anndata as ad,torch,scvi
 torch.set_num_threads(cfg['threads']);sys.path.insert(0,str(R/'scripts'));from audit_TME_full_cluster_markers_v1 import read_rows
 for code,r in enumerate(inp['records']):
  where=np.flatnonzero(axes['source_code'][probes]==code);take=probes[where];z=zarr.open_group(r['store'],mode='r');var=ad.io.read_elem(z['var']);names=[str(var.iloc[c]['gene_symbol']) if 'gene_symbol' in var else str(var.index[c]) for c in r['selected_columns']];assert names==inp['genes'];raw=read_rows(z[r['counts_key']],axes['native_row'][take],True);np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),axes['native_total_RNA'][take]);X[where]=raw[:,r['selected_columns']].toarray().astype(np.float32)
 assert (X.sum(axis=1)>0).all();a=ad.AnnData(sp.csr_matrix(X),obs=pd.DataFrame({'source':pd.Categorical.from_codes(axes['source_code'][probes],[r['source'] for r in inp['records']])}),var=pd.DataFrame(index=inp['genes']));model=scvi.model.SCVI.load(str(P/'outputs/model'),adata=a,accelerator='cpu',device='auto');assert str(next(model.module.parameters()).device)=='cpu';check=np.asarray(model.get_latent_representation(batch_size=cfg['inference_batch_size']),np.float32);np.testing.assert_allclose(check,latent[probes],rtol=cfg['CPU_GPU_latent_rtol'],atol=cfg['CPU_GPU_latent_atol'])
 # Actual saved train/validation indices; not just intended fractions.
 exposure=np.full(N,-1,np.int8);counts={}
 for code,name in enumerate(['train_indices','validation_indices','test_indices']):
  v=getattr(model,name)
  if v is None:v=np.empty(0,np.int64)
  v=np.asarray(v,np.int64);assert len(v)==len(np.unique(v)) and ((v>=0)&(v<N)).all() and (exposure[v]==-1).all();exposure[v]=code;counts[name]=len(v)
 assert (exposure>=0).all() and counts['test_indices']==0 and counts['train_indices']>0 and counts['validation_indices']>0;np.save(R/'outputs/private_scVI_gradient_exposure_v1.npy',exposure)
 save(R/'outputs/accepted_scVI_saved_model_v1.json',dict(n_fit=N,n_sources=len(inp['records']),n_independent_native_count_CPU_model_probes=len(probes),all_latent_rows_and_file_SHA_verified=True,actual_GPU_epochs=len(epochs),actual_full_KL_weight_epochs=sum(r['kl_weight']>=.999 for r in epochs),training_and_validation_history_finite=True,saved_gradient_exposure_counts=counts,exposure_from_actual_saved_indices_not_reconstructed=True,internal_random_cell_validation_not_independent_patient_or_source_validation=True,independent_CPU_inference_matches_GPU_latent=True,maximum_absolute_CPU_GPU_latent_error=float(np.max(np.abs(check-latent[probes]))),CPU_GPU_precision_tolerances=dict(rtol=cfg['CPU_GPU_latent_rtol'],atol=cfg['CPU_GPU_latent_atol']),training_was_GPU_only=True,biological_integration_quality_and_final_identity_not_yet_accepted=True,producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json')))
 save(R/'provenance/terminal_v1.json',dict(status='WHOLE_SCVI_LATENT_HISTORIES_ACTUAL_EXPOSURE_AND_INDEPENDENT_NATIVE_CPU_RELOAD_ACCEPTED_REQUIRES_BIOLOGICAL_COMPARISON',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json')))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
