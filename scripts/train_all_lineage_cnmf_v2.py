"""Official cNMF all18 original-major contrasts; no prior or outcome selection."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc,re,subprocess,importlib.metadata
import numpy as np,anndata as ad
from cnmf import cNMF
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def free():
 v=os.statvfs(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);return v.f_bavail*v.f_frsize
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(os.environ['PANCANCER_ALL_LINEAGE_CNMF_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'] and importlib.metadata.version('cnmf')=='1.7.0'
A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';iv=json.loads(A.read_text());assert iv['status']=='ACCEPTED_ALL18_NATIVE_DISCOVERY_INPUTS_NOT_FORMAL_STATES' and iv['cells']==100766
assert sha(A)==plan['accepted_input_inventory_sha256'] and sha(A.parent/'independent_audit.json')==plan['accepted_input_audit_sha256']
assert json.loads((A.parent/'independent_audit.json').read_text())['status'].startswith('PASS_ALL100766')
assert free()>=(plan['reserve_GiB']+plan['model_budget_GiB'])*2**30
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'all_lineage_cNMF_state_contrast_20261014_v2';ledger=json.loads((V/'recovery_input_inventory.json').read_text());assert ledger['scientific_parameters_unchanged'];start=datetime.now(timezone.utc).isoformat();summaries=json.loads((O/'imported_v1_progress.json').read_text())['summary'];done={r['major'] for r in summaries}
for name,h in ledger['imported_file_sha256']['cNMF'].items():assert sha(O/name)==h
for i,item in enumerate(sorted(iv['files'],key=lambda a:a['major'])):
 if item['major'] in done:
  print(json.dumps(dict(major=item['major'],phase='REUSED_SHA_VERIFIED_ALL_V1_K_REPEATS_CONSENSUS')),flush=True);continue
 assert free()>=(plan['reserve_GiB']+1)*2**30
 assert sha(item['path'])==item['sha256'];raw=ad.read_h5ad(item['path']);major=item['major'];assert set(raw.obs.major.astype(str))=={major};root=O/('major_%02d'%i);root.mkdir()
 # Reuse the accepted count file; only native-fullRNA TPM is an extra copy.
 tpm=raw.copy();tpm.X=raw.X.astype(np.float64).multiply((1000000/raw.obs.full_RNA_total.to_numpy())[:,None]).tocsr();tpm.write_h5ad(root/'private_native_fullRNA_TPM.h5ad',compression='gzip',compression_opts=4);del tpm;gc.collect()
 model=cNMF(output_dir=str(root),name='cnmf');model.prepare(counts_fn=item['path'],tpm_fn=str(root/'private_native_fullRNA_TPM.h5ad'),components=plan['K'],n_iter=plan['repeats_per_K'],seed=plan['seed']+i,num_highvar_genes=plan['HVG_per_major'],max_NMF_iter=plan['max_NMF_iter'],init='random',beta_loss='frobenius')
 print(json.dumps(dict(major=major,phase='FACTORIZING',K=plan['K'],repeats=plan['repeats_per_K'],started_UTC=datetime.now(timezone.utc).isoformat())),flush=True);model.factorize(worker_i=0,total_workers=1,skip_completed_runs=False);model.combine(components=plan['K'],skip_missing_files=False)
 for k in plan['K']:model.consensus(k=k,density_threshold=.5,local_neighborhood_size=.3,show_clustering=False)
 r=dict(major=major,cells=raw.n_obs,genes=raw.n_vars,K=plan['K'],repeats=plan['repeats_per_K'],all_consensus_completed=True,finished_UTC=datetime.now(timezone.utc).isoformat(),private_relative_root='major_%02d'%i);summaries.append(r);(O/'progress.json').write_text(json.dumps(dict(status='ALL18_CNMF_PENDING_STABILITY',summary=summaries),indent=2)+'\n');print(json.dumps(r),flush=True);del raw,model;gc.collect()
out=dict(status='ALL18_CNMFK_REPEATS_CONSENSUS_PENDING_COMPARISON',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),accepted_inventory_sha256=sha(A),cells=100766,majors=18,summary=summaries,no_prior_or_clinical_LR_used=True,formal_states=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='outputs'}),flush=True)
