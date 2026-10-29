"""Official GPU scVI on fixed native counts/axes/HVG; no Harmony chaining."""
import argparse,datetime,gc,hashlib,importlib.metadata,inspect,json,os,re,subprocess,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def check(R,stage,**v):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**v))
def gather_counts(csr,rows):
 unique,inverse=np.unique(rows,return_inverse=True);parts=[];low=0
 while low<len(unique):
  high=int(np.searchsorted(unique,unique[low]+2048));first=int(unique[low]);last=int(unique[high-1])+1
  ptr=np.asarray(csr['indptr'][first:last+1],np.int64);start=int(ptr[0]);end=int(ptr[-1]);v=np.asarray(csr['data'][start:end],np.float64);ix=np.asarray(csr['indices'][start:end],np.int64)
  if not np.isfinite(v).all() or (v<0).any() or not np.equal(v,np.floor(v)).all() or v.max(initial=0)>=2**24:raise ValueError('Invalid native integer counts')
  x=sp.csr_matrix((v,ix,ptr-start),shape=(last-first,int(csr.attrs['shape'][1])));x.sum_duplicates();parts.append(x[unique[low:high]-first]);low=high
 return sp.vstack(parts,format='csr')[inverse]
def main(R):
 started=time.monotonic();cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor'])
 for D in [P,A]:
  while not (D/'provenance/supervisor_terminal_v1.json').exists():time.sleep(5)
  assert load(D/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 import torch,scvi,anndata as ad,zarr,psutil
 from lightning.pytorch.callbacks import Callback
 assert torch.cuda.is_available(),'GPU required; never silently fall back to CPU'
 assert os.environ.get('SLURM_JOB_ID'),'Scheduled GPU allocation required'
 torch.set_num_threads(cfg['threads']);scvi.settings.seed=cfg['seed'];scvi.settings.dl_num_workers=0
 torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
 free,total=torch.cuda.mem_get_info();assert free>=cfg['minimum_GPU_free_bytes']
 versions={k:importlib.metadata.version(k) for k in ['scvi-tools','torch','anndata','numpy','scipy','lightning','zarr']}
 save(R/'provenance/runtime_v1.json',dict(versions=versions,GPU_name=torch.cuda.get_device_name(0),GPU_free_bytes=free,CUDA=torch.version.cuda,SLURM_JOB_ID=os.environ['SLURM_JOB_ID'],SCVI_API=str(inspect.signature(scvi.model.SCVI)),train_API=str(inspect.signature(scvi.model.SCVI.train))))
 # Actual GPU training, serialization, and reload before reading full inputs.
 rng=np.random.default_rng(cfg['seed']);toy=ad.AnnData(sp.csr_matrix(rng.poisson(2,(96,40)).astype(np.float32)));toy.obs['source']=pd.Categorical(np.repeat(['a','b','c'],32));scvi.model.SCVI.setup_anndata(toy,batch_key='source');m=scvi.model.SCVI(toy,n_latent=3,n_layers=1,n_hidden=16,gene_likelihood=cfg['gene_likelihood'],dispersion='gene')
 m.train(max_epochs=2,accelerator='gpu',devices=1,batch_size=32,train_size=.8,validation_size=.2,enable_progress_bar=False,enable_checkpointing=False,default_root_dir=str(R/'provenance/toy_logs'))
 assert str(next(m.module.parameters()).device).startswith('cuda');z=m.get_latent_representation();assert np.isfinite(z).all();m.save(str(R/'provenance/toy_model'),save_anndata=False);reload=scvi.model.SCVI.load(str(R/'provenance/toy_model'),adata=toy,accelerator='gpu',device=0);np.testing.assert_allclose(reload.get_latent_representation(),z,rtol=2e-5,atol=2e-5);save(R/'provenance/actual_GPU_SCVI_smoke_v1.json',dict(actual_GPU_training_and_saved_model_reload_passed=True,synthetic_not_biological_validation=True));del m,reload,toy;gc.collect();torch.cuda.empty_cache()
 # Prepare only column metadata in the original native reader environment.
 spec=R/'outputs/private_native_input_spec_v1.json';subprocess.run([cfg['native_python'],'-B',str(R/'scripts/prepare_TME_scVI_native_feature_map_v1.py'),'--producer',str(P),'--output',str(spec)],check=True)
 inp=load(spec);assert sha(inp['axes_path'])==inp['axes_SHA256'];assert psutil.virtual_memory().available>=cfg['minimum_RAM_available_bytes']
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg)
 with np.load(inp['axes_path'],allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 chunks=[]
 for code,r in enumerate(inp['records']):
  wh=np.flatnonzero(axes['source_code']==code);assert len(wh)==r['n_fit'] and wh[0]==r['offset'];z=zarr.open_group(r['store'],mode='r');cols=r['selected_columns']
  for lo in range(0,len(wh),cfg['read_block_size']):
   take=wh[lo:lo+cfg['read_block_size']];raw=gather_counts(z[r['counts_key']],axes['native_row'][take]);np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),axes['native_total_RNA'][take]);block=raw[:,cols].astype(np.float32).tocsr();assert (np.asarray(block.sum(axis=1)).ravel()>0).all();chunks.append(block)
  check(R,'native_integer_scVI_source_complete',source=r['source'],n_fit=len(wh));gc.collect()
 X=sp.vstack(chunks,format='csr');del chunks;gc.collect();assert X.shape==(inp['n_fit'],len(inp['genes']))
 a=ad.AnnData(X,obs=pd.DataFrame({'source':pd.Categorical.from_codes(axes['source_code'],[r['source'] for r in inp['records']])}),var=pd.DataFrame(index=inp['genes']));del X;scvi.model.SCVI.setup_anndata(a,batch_key='source');scvi.settings.seed=cfg['seed'];model=scvi.model.SCVI(a,n_latent=cfg['n_latent'],n_layers=cfg['n_layers'],n_hidden=cfg['n_hidden'],gene_likelihood=cfg['gene_likelihood'],dispersion='gene')
 class EpochReceipt(Callback):
  def on_train_epoch_end(self,trainer,pl_module):
   metrics={k:float(v.detach().cpu()) for k,v in trainer.callback_metrics.items() if getattr(v,'numel',lambda:0)()==1 and bool(torch.isfinite(v).all())}
   check(R,'official_GPU_scVI_epoch_complete',epoch=int(trainer.current_epoch)+1,device=str(next(pl_module.parameters()).device),metrics=metrics,GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated(),kl_weight=float(pl_module.kl_weight))
 check(R,'full_GPU_SCVI_training_started',n_fit=a.n_obs,n_features=a.n_vars,max_epochs=cfg['max_epochs'])
 model.train(max_epochs=cfg['max_epochs'],accelerator='gpu',devices=1,batch_size=cfg['batch_size'],train_size=.9,validation_size=.1,shuffle_set_split=True,early_stopping=False,plan_kwargs={'n_epochs_kl_warmup':cfg['n_epochs_kl_warmup']},check_val_every_n_epoch=1,enable_progress_bar=False,enable_checkpointing=False,default_root_dir=str(R/'provenance/train_logs'),callbacks=[EpochReceipt()])
 device=str(next(model.module.parameters()).device);assert device.startswith('cuda');assert float(model.trainer.lightning_module.kl_weight)>=.999;model.save(str(R/'outputs/model'),save_anndata=False)
 history=model.history or {}
 for key,value in history.items():
  assert re.fullmatch(r'[A-Za-z0-9_]+',key);pd.DataFrame(value).to_csv(R/'outputs'/('history_'+key+'.tsv'),sep='\t')
 latent=np.asarray(model.get_latent_representation(batch_size=cfg['batch_size']),np.float32);assert latent.shape==(inp['n_fit'],cfg['n_latent']) and np.isfinite(latent).all();np.save(R/'outputs/X_scVI_v1.npy',latent)
 probe=np.unique(np.linspace(0,len(latent)-1,512).astype(int));loaded=scvi.model.SCVI.load(str(R/'outputs/model'),adata=a,accelerator='gpu',device=0);np.testing.assert_allclose(loaded.get_latent_representation(indices=probe,batch_size=cfg['batch_size']),latent[probe],rtol=3e-5,atol=3e-5)
 files=[dict(path=str(p.relative_to(R/'outputs')),bytes=p.stat().st_size,SHA256=sha(p)) for p in (R/'outputs').rglob('*') if p.is_file()];storage_gate(R,cfg);save(R/'outputs/scVI_GPU_manifest_v1.json',dict(files=files,n_fit=inp['n_fit'],n_latent=cfg['n_latent'],axes_SHA256=inp['axes_SHA256'],genes_SHA256=inp['genes_SHA256'],actual_device=device,expression_matrix_copied_to_disk=False,original_counts_labels_and_QC_unchanged=True))
 save(R/'provenance/terminal_v1.json',dict(status='FULL_CURRENT_TME_GPU_SCVI_TRAINED_AND_SAVED_MODEL_RELOAD_ACCEPTED_REQUIRES_GRAPH_AND_BIOLOGICAL_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=inp['n_fit'],actual_device=device,max_epochs=cfg['max_epochs'],GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated(),wall_seconds=time.monotonic()-started,manifest_SHA256=sha(R/'outputs/scVI_GPU_manifest_v1.json')))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
