"""Prespecified identity-aware Spectra blocks, official GPU algorithm unchanged."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc,random,re,subprocess
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,torch
from scipy.sparse import vstack
from Spectra.Spectra_gpu import est_spectra
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def free():
 v=os.statvfs(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);return v.f_bavail*v.f_frsize
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_ALL_LINEAGE_SPECTRA_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256']
A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';iv=json.loads(A.read_text());assert iv['status']=='ACCEPTED_ALL18_NATIVE_DISCOVERY_INPUTS_NOT_FORMAL_STATES' and iv['cells']==100766
assert sha(A)==plan['accepted_input_inventory_sha256'] and sha(A.parent/'independent_audit.json')==plan['accepted_input_audit_sha256']
assert json.loads((A.parent/'independent_audit.json').read_text())['status'].startswith('PASS_ALL100766')
assert free()>=(plan['reserve_GiB']+plan['model_budget_GiB'])*2**30
torch.set_num_threads(2);assert torch.cuda.is_available()
runtime=R/'state_model_runtime_20261014_v1';resources=runtime/'resources';receipt=json.loads((resources/'runtime_resource_receipts.json').read_text());assert receipt['Spectra_repository_commit']==plan['Spectra_commit']
for item in receipt['source_files']:assert sha(resources/'spectra_snapshot'/item['path'])==item['sha256']
ref=json.loads((resources/'spectra_snapshot/src/Spectra/data/default_gene_sets.json').read_text());lookup={i['major']:i for i in iv['files']};V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'all_lineage_Spectra_state_discovery_20261014_v2';ledger=json.loads((V/'recovery_input_inventory.json').read_text());assert ledger['scientific_parameters_unchanged'];start=datetime.now(timezone.utc).isoformat();summaries=json.loads((O/'imported_v1_progress.json').read_text())['summary'];done={(r['block'],r['seed']) for r in summaries}
for name,h in ledger['imported_file_sha256']['Spectra'].items():assert sha(O/name)==h
assert sorted(m for b in plan['blocks'].values() for m in b)==sorted(lookup)
for bi,(block,majors) in enumerate(plan['blocks'].items()):
 if all((block,seed) in done for seed in plan['seeds']):
  print(json.dumps(dict(block=block,phase='REUSED_SHA_VERIFIED_COMPLETE_V1_BLOCK')),flush=True);continue
 assert free()>=(plan['reserve_GiB']+1)*2**30
 objs=[]
 for major in majors:
  item=lookup[major];assert sha(item['path'])==item['sha256'];a=ad.read_h5ad(item['path']);assert set(a.obs.major.astype(str))=={major};objs.append(a)
 common=sorted(set.intersection(*[set(a.var_names) for a in objs]));assert len(common)>=plan['min_common_genes']
 mats=[];frames=[]
 for a in objs:
  x=a[:,common].X.astype(np.float32).multiply((10000/a.obs.full_RNA_total.to_numpy())[:,None]).tocsr();x.data=np.log1p(x.data);mats.append(x);frames.append(a.obs.copy())
 base=ad.AnnData(X=vstack(mats,format='csr'),obs=pd.concat(frames),var=pd.DataFrame(index=pd.Index(common,name='gene')));assert base.obs_names.is_unique
 dictionary={'global':ref['global']}
 for major in majors:dictionary[major]={scope+'__'+name:gs for scope in plan['prior_mapping'][major] for name,gs in ref[scope].items()}
 axis=set(common);dictionary={scope:{n:sorted(set(gs)&axis) for n,gs in sets.items() if len(set(gs)&axis)>=plan['min_reference_genes']} for scope,sets in dictionary.items()}
 hvgs=set();cap=int(np.floor(plan['HVG_budget_per_block']/len(majors)))
 for major in majors:
  part=base[base.obs.major.eq(major)].copy();sc.pp.highly_variable_genes(part,flavor='seurat',batch_key='gse_id',n_top_genes=cap);hvgs.update(part.var_names[part.var.highly_variable]);del part
 prior=set(gene for sets in dictionary.values() for gs in sets.values() for gene in gs);features=sorted(hvgs|prior)
 L={scope:len(sets)+(plan['global_free_factors'] if scope=='global' else plan['lineage_free_factors_with_prior'] if sets else plan['lineage_free_factors_without_prior']) for scope,sets in dictionary.items()}
 # Conservative GPU tensor allowance incl. dense adjacency/gradient temporaries.
 estimated=12*len(L)*len(features)**2*4+base.n_obs*len(features)*4*2+base.n_obs*sum(L.values())*4*8
 gpu_free,gpu_total=torch.cuda.mem_get_info();assert estimated<=gpu_free*plan['GPU_free_fraction_limit'],(block,estimated,gpu_free)
 data=base[:,features].copy();data.X=data.X.toarray().astype(np.float32);root=O/block;fixed=dict(block=block,majors=majors,common_genes=common,features=features,HVG=sorted(hvgs),reference=dictionary,L=L,estimated_GPU_bytes=estimated,free_GPU_bytes_at_start=gpu_free)
 if root.exists():
  oldfixed=json.loads((root/'fixed_feature_dictionary.json').read_text())
  for field in ['block','majors','common_genes','features','HVG','reference','L']:assert oldfixed[field]==fixed[field],field
 else:root.mkdir();(root/'fixed_feature_dictionary.json').write_text(json.dumps(fixed,indent=2)+'\n')
 del objs,base,mats,frames;gc.collect()
 for seed in plan['seeds']:
  if (block,seed) in done:
   print(json.dumps(dict(block=block,seed=seed,phase='REUSED_SHA_VERIFIED_COMPLETE_V1_SEED')),flush=True);continue
  assert free()>=(plan['reserve_GiB']+1)*2**30
  random.seed(seed+bi*100);np.random.seed(seed+bi*100);torch.manual_seed(seed+bi*100);torch.cuda.manual_seed_all(seed+bi*100)
  a=data.copy();beg=datetime.now(timezone.utc).isoformat();print(json.dumps(dict(block=block,seed=seed,phase='TRAINING',cells=a.n_obs,features=len(features),L=L,started_UTC=beg)),flush=True)
  model=est_spectra(a,dictionary,L=L,use_highly_variable=False,cell_type_key='major',filter_sets=False,n_top_vals=50,lam=.008,delta=.001,rho=.05,num_epochs=plan['num_epochs'],batch_size=plan['batch_size'],skip=15,verbose=False)
  scores=np.asarray(a.obsm['SPECTRA_cell_scores']);weights=np.asarray(a.uns['SPECTRA_factors']);mask=model.internal_model.alpha_mask.detach().cpu().numpy();assert scores.shape==(a.n_obs,sum(L.values())) and weights.shape==(sum(L.values()),len(features));assert np.isfinite(scores).all() and np.isfinite(weights).all() and (scores>=0).all() and (weights>=0).all()
  scopes=[];order=[]
  for scope in model.internal_model.ct_order:order.extend([scope]*L[scope])
  for j,scope in enumerate(order):
   expected=np.ones(a.n_obs,bool) if scope=='global' else a.obs.major.eq(scope).to_numpy();assert np.array_equal(mask[:,j]>0,expected)
   scopes.append('block_global' if scope=='global' else scope)
  dest=root/('seed'+str(seed));dest.mkdir();np.savez_compressed(dest/'factors_and_discovery_scores.npz',factors=weights,cell_scores=scores,cell_ids=np.asarray(a.obs_names,dtype=str),genes=np.asarray(features,dtype=str),factor_scope=np.asarray(scopes,dtype=str));torch.save(model.internal_model.state_dict(),dest/'model_state_dict.pt');a.obs.to_csv(dest/'private_discovery_metadata.tsv.gz',sep='\t',index=True,compression='gzip')
  r=dict(block=block,seed=seed,actual_seed=seed+bi*100,cells=a.n_obs,features=len(features),factors=sum(L.values()),started_UTC=beg,finished_UTC=datetime.now(timezone.utc).isoformat(),block_shared_does_not_imply_all18_shared=True);summaries.append(r);(O/'progress.json').write_text(json.dumps(dict(status='ALL_LINEAGE_SPECTRA_FIT_PENDING_STABILITY',summary=summaries),indent=2)+'\n');print(json.dumps(r),flush=True);del a,model,scores,weights,mask;gc.collect();torch.cuda.empty_cache()
 del data;gc.collect()
out=dict(status='ALL6_BLOCKS_3SEEDS_SPECTRA_CANDIDATES_PENDING_STABILITY',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),accepted_inventory_sha256=sha(A),cells=100766,majors=18,summary=summaries,formal_states=False,all9CM_mapping_completed=False,DIALOGUE_completed=False,spatial_completed=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='outputs'}),flush=True)
