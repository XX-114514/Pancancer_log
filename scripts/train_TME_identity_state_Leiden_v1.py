"""Balanced development reference for identity-conditioned expression patterns.

Official Scanpy Leiden partitions fixed programme coordinates, never UMAP.
Pattern integers are computational partitions; biological programmes may overlap.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
import joblib

def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def binding(p):return dict(path=str(p),SHA256=sha(p))
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def wait(root,seconds,start):
 p=Path(root)/'provenance/supervisor_terminal_v1.json'
 while not p.exists():
  assert time.monotonic()-start<seconds;time.sleep(5)
 assert load(p)['actual_child_exit_code']==0

def main(R):
 cfg=load(R/'config/spec_v1.json');start=time.monotonic()
 for root in [cfg['inputs'],cfg['sensitivity'],cfg['math_test_stage']]:wait(root,cfg['dependency_wait_seconds'],start)
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 import scanpy as sc,anndata,sklearn
 storage_gate(R,cfg);I=Path(cfg['inputs']);ip=I/'outputs/accepted_development_state_multiview_inputs_v1.json';assert sha(ip)==load(I/'provenance/terminal_v1.json')['acceptance_SHA256'];mi=load(ip)
 candidate=load(bound(mi['candidate_acceptance']));readouts=load(bound(candidate['readout_acceptance']));app=load(bound(mi['applicability']));roles=load(bound(mi['source_roles']));audit=load(bound(readouts['score_numerical_audit']));idx=load(bound(audit['producer_index']));names=idx['program_ids'];index={n:j for j,n in enumerate(names)}
 role={r['source']:r['role'] for r in roles['source_roles']};ctx=np.load(bound(app['context_codes']));allowed=np.load(bound(roles['development_mask']));sources=[]
 for rec in idx['sources']:
  if role[rec['source']]!='development':continue
  with np.load(bound(rec['cell_axis']),allow_pickle=False) as z:rows=z['fit_row'];units=z['qualification_unit_code']
  assert allowed[rows].all()
  sources.append(dict(record=rec,rows=rows,units=units,context=ctx[rows],coverage=np.array([r['native_coverage_pass'] for r in rec['coverage_and_interpretation']],bool)))
 models=[];held=[];rng=np.random.default_rng(cfg['sampling_seed'])
 for view in mi['views']:
  c=view['context']['code'];features=view['program_ids'];columns=np.array([index[n] for n in features],int)
  with np.load(bound(view['values']),allow_pickle=False) as zz:eligible_units=zz['unit_codes'][zz['n_context_cells']>=cfg['minimum_context_cells']]
  context_masks=[(s['context']==c)&np.isin(s['units'],eligible_units) for s in sources]
  adequate=[s for s,mk in zip(sources,context_masks) if np.count_nonzero(mk)>=cfg['minimum_source_context_cells']]
  if len(adequate)<cfg['minimum_training_sources']:
   held.append(dict(context_code=c,reason='insufficient_development_sources'));continue
  complete=np.logical_and.reduce([s['coverage'][columns] for s in adequate]);omitted=[features[i] for i in np.flatnonzero(~complete)];columns=columns[complete];features=[features[i] for i in np.flatnonzero(complete)]
  if len(features)<cfg['minimum_programmes']:
   held.append(dict(context_code=c,reason='insufficient_jointly_measured_programmes',n_programmes=len(features)));continue
  row=[];sourcecodes=[];unit=[];values=[]
  for sn,s in enumerate(sources):
   if not s['coverage'][columns].all():continue
   choose=[]
   for u in np.unique(s['units'][context_masks[sn]]):
    positions=np.flatnonzero(context_masks[sn]&(s['units']==u));choose.extend(rng.choice(positions,min(len(positions),cfg['maximum_cells_per_unit']),replace=False).tolist())
   if len(choose)>cfg['maximum_cells_per_source']:choose=rng.choice(choose,cfg['maximum_cells_per_source'],replace=False).tolist()
   if not choose:continue
   choose=np.sort(choose);score=np.load(bound(s['record']['scores']),mmap_mode='r');yy=np.asarray(score[choose][:,columns]);assert np.isfinite(yy).all()
   values.append(yy);row.extend(s['rows'][choose].tolist());unit.extend(s['units'][choose].tolist());sourcecodes.extend([sn]*len(choose))
  if len(row)<cfg['minimum_training_cells']:
   held.append(dict(context_code=c,reason='insufficient_training_cells',n_cells=len(row)));continue
  row=np.array(row,np.int64);unit=np.array(unit,np.int32);sourcecodes=np.array(sourcecodes,np.int32);y=np.vstack(values)
  if len(row)>cfg['maximum_cells_per_context']:
   # Round-robin balanced over sources; no score-dependent sampling.
   strata=[rng.permutation(np.flatnonzero(sourcecodes==s)) for s in np.unique(sourcecodes)];chosen=np.array([a[k] for k in range(max(map(len,strata))) for a in strata if len(a)>k][:cfg['maximum_cells_per_context']]);row=row[chosen];unit=unit[chosen];sourcecodes=sourcecodes[chosen];y=y[chosen]
  keep=np.ptp(y,axis=0)>0
  if keep.sum()<cfg['minimum_programmes']:
   held.append(dict(context_code=c,reason='insufficient_nonconstant_programmes'));continue
  constant=[features[i] for i in np.flatnonzero(~keep)];features=[features[i] for i in np.flatnonzero(keep)];columns=columns[keep];y=y[:,keep]
  assert max(np.unique(unit,return_counts=True)[1])<=cfg['maximum_cells_per_unit'];assert max(np.unique(sourcecodes,return_counts=True)[1])<=cfg['maximum_cells_per_source'];assert len(row)<=cfg['maximum_cells_per_context']
  scale=StandardScaler();z=scale.fit_transform(y);pca=PCA(n_components=min(cfg['maximum_PCs'],z.shape[1],z.shape[0]-1),svd_solver='full');pc=pca.fit_transform(z)
  a=anndata.AnnData(z.astype(np.float32),var=pd.DataFrame(index=features));a.obsm['X_state_PCA']=pc.astype(np.float32);sc.pp.neighbors(a,n_neighbors=cfg['n_neighbors'],use_rep='X_state_PCA',random_state=cfg['graph_seed'])
  labels={};classifiers={};stability=[]
  for resolution in cfg['resolutions']:
   for seed in cfg['seeds']:
    key='r'+str(resolution).replace('.','p')+'_seed'+str(seed);sc.tl.leiden(a,resolution=resolution,random_state=seed,key_added=key,flavor='igraph',directed=False,n_iterations=2);labs=a.obs[key].astype(str).to_numpy();labels[key]=labs
    classifier=KNeighborsClassifier(n_neighbors=cfg['n_neighbors'],weights='distance',n_jobs=cfg['threads']);classifier.fit(pc,labs);classifiers[key]=classifier
    stability.append(dict(partition=key,n_patterns=len(set(labs)),source_NMI=float(normalized_mutual_info_score(sourcecodes,labs)),unit_NMI=float(normalized_mutual_info_score(unit,labs))))
  for res in cfg['resolutions']:
   ka='r'+str(res).replace('.','p')+'_seed'+str(cfg['seeds'][0]);kb='r'+str(res).replace('.','p')+'_seed'+str(cfg['seeds'][1]);stability.append(dict(resolution=res,seed_adjusted_Rand=float(adjusted_rand_score(labels[ka],labels[kb]))))
  out=R/'outputs'/('private_context_'+str(c)+'_state_partition_reference_v1.joblib');joblib.dump(dict(context_code=c,program_ids=features,program_columns=columns,scaler=scale,PCA=pca,classifiers=classifiers,primary_key='r'+str(cfg['primary_resolution']).replace('.','p')+'_seed'+str(cfg['seeds'][0]),reference_fit_rows=row,reference_labels=labels,reference_source_codes=sourcecodes,reference_unit_codes=unit),out,compress=3)
  graph=R/'outputs'/('private_context_'+str(c)+'_state_partition_graph_v1.h5ad');a.obs_names=[str(int(i)) for i in row];a.write_h5ad(graph,compression='gzip')
  profiles=[]
  for key,labs in labels.items():
   for lab in sorted(set(labs),key=int):
    mask=labs==lab;delta=(y[mask].mean(0)-y.mean(0))/np.maximum(y.std(0),1e-12)
    profiles.append(dict(partition=key,pattern_integer=int(lab),n_training_cells=int(mask.sum()),n_development_sources=len(np.unique(sourcecodes[mask])),n_candidate_work_units=len(np.unique(unit[mask])),mean_programme_scores=y[mask].mean(0).tolist(),standardized_programme_difference=delta.tolist(),functional_name_not_assigned=True))
  models.append(dict(context=view['context'],n_training_cells=len(row),n_training_sources=len(np.unique(sourcecodes)),n_programmes=len(features),program_ids=features,coverage_omitted_programmes=omitted,constant_programmes=constant,reference_model=binding(out),reference_graph=binding(graph),profiles=profiles,stability=stability,reference_source_labels=[s['record']['source'] for s in sources]))
  storage_gate(R,cfg);save(R/'checkpoints'/('context_'+str(c)+'_complete_v1.json'),dict(context_code=c,n_completed_contexts=len(models),n_training_cells=len(row)))
 op=R/'outputs/accepted_development_identity_state_partitions_v1.json';save(op,dict(schema='TME_development_identity_conditioned_state_partition_models_v1',input_acceptance=binding(ip),source_roles=mi['source_roles'],applicability=mi['applicability'],score_numerical_acceptance=readouts['score_numerical_audit'],n_context_models=len(models),models=models,held_contexts=held,software=dict(scanpy=sc.__version__,sklearn=sklearn.__version__),primary_resolution=cfg['primary_resolution'],resolutions=cfg['resolutions'],seeds=cfg['seeds'],method='StandardScaler_PCA_Scanpy_neighbors_Leiden_then_fixed_distance_weighted_kNN_projection',reference_sampling_independent_of_programme_scores=True,all_nonconstant_jointly_measured_features_no_variance_ranking=True,n_lockbox_score_arrays_opened=0,clinical_outcomes_used=False,patterns_are_computational_partitions_not_exclusive_biological_states=True,identity_and_continuous_vectors_preserved=True,projection_and_state_CoVarNet_not_yet_executed=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='DEVELOPMENT_IDENTITY_STATE_PATTERN_REFERENCES_FROZEN',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
