"""Independent native-score/axis replay of frozen developmental state graphs.

Audits numerics, sampling and source-role boundaries, not biological function.
The official stored graph partitions remain unchanged, including low support.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,joblib,anndata
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score

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
def main(R):
 cfg=load(R/'config/spec_v1.json');T=Path(cfg['training']);start=time.monotonic()
 while not (T/'provenance/supervisor_terminal_v1.json').exists():
  assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
 assert load(T/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 ap=T/'outputs/accepted_development_identity_state_partitions_v1.json';assert sha(ap)==load(T/'provenance/terminal_v1.json')['acceptance_SHA256'];m=load(ap);tc=load(T/'config/spec_v1.json');threadpool_limits(limits=tc['threads'])
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);roles=load(bound(m['source_roles']));development=np.load(bound(roles['development_mask']));app=load(bound(m['applicability']));context=np.load(bound(app['context_codes']));read_audit=load(bound(m['score_numerical_acceptance']));idx=load(bound(read_audit['producer_index']));sources={r['source']:r for r in idx['sources']};names=idx['program_ids'];columns={v:j for j,v in enumerate(names)}
 inp=load(bound(m['input_acceptance']));inputviews={r['context']['code']:r for r in inp['views']};role={r['source']:r['role'] for r in roles['source_roles']};reports=[];total=0
 for item in m['models']:
  c=item['context']['code'];model=joblib.load(bound(item['reference_model']));rows=model['reference_fit_rows'];u=model['reference_unit_codes'];sourcecodes=model['reference_source_codes'];n=len(rows);assert len(np.unique(rows))==n and development[rows].all() and (context[rows]==c).all()
  assert len(rows)<=tc['maximum_cells_per_context'] and np.unique(u,return_counts=True)[1].max()<=tc['maximum_cells_per_unit'] and np.unique(sourcecodes,return_counts=True)[1].max()<=tc['maximum_cells_per_source']
  yy=np.full((n,len(model['program_ids'])),np.nan,np.float32);source_labels=item['reference_source_labels']
  for s in np.unique(sourcecodes):
   source=source_labels[s];assert role[source]=='development';record=sources[source];at=np.flatnonzero(sourcecodes==s)
   with np.load(bound(record['cell_axis']),allow_pickle=False) as z:axis=z['fit_row'];uu=z['qualification_unit_code']
   positions=np.searchsorted(axis,rows[at]);np.testing.assert_array_equal(axis[positions],rows[at]);np.testing.assert_array_equal(uu[positions],u[at])
   y=np.load(bound(record['scores']),mmap_mode='r');js=[columns[k] for k in model['program_ids']];yy[at]=y[positions][:,js]
  yy=np.asfortranarray(yy);assert np.isfinite(yy).all();np.testing.assert_allclose(model['scaler'].mean_,yy.mean(0,dtype=np.float64),rtol=1e-10,atol=1e-10);np.testing.assert_allclose(model['scaler'].var_,yy.var(0,dtype=np.float64),rtol=1e-9,atol=1e-10)
  # StandardScaler transforms float32 in-place, rounding after each
  # operation. PCA.fit_transform stores U*S; it is not bitwise the later
  # (X-mean)@V projection. Re-execute the original full-SVD algorithm.
  z=yy.copy(order='F');z-=model['scaler'].mean_;z/=model['scaler'].scale_
  graph=anndata.read_h5ad(bound(item['reference_graph']));np.testing.assert_array_equal(graph.X,z)
  replay=PCA(n_components=model['PCA'].n_components_,svd_solver='full');pc=replay.fit_transform(z)
  np.testing.assert_array_equal(replay.components_,model['PCA'].components_)
  np.testing.assert_array_equal(replay.singular_values_,model['PCA'].singular_values_)
  assert graph.obs_names.tolist()==[str(int(v)) for v in rows];assert graph.var_names.tolist()==model['program_ids'];np.testing.assert_array_equal(graph.X,z);np.testing.assert_array_equal(graph.obsm['X_state_PCA'],pc)
  co=graph.obsp['connectivities'];assert co.shape==(n,n) and np.isfinite(co.data).all() and co.data.min(initial=0)>=0;assert abs(co-co.T).max()<1e-6
  with np.load(bound(inputviews[c]['values']),allow_pickle=False) as z:allowed=set(z['unit_codes'][z['n_context_cells']>=tc['minimum_context_cells']].tolist())
  assert set(u)<=allowed
  replay=[]
  for key,labs in model['reference_labels'].items():
   np.testing.assert_array_equal(graph.obs[key].astype(str).to_numpy(),labs);assert set(labs)==set(model['classifiers'][key].classes_)
   for p in [v for v in item['profiles'] if v['partition']==key]:
    use=labs==str(p['pattern_integer']);assert int(use.sum())==p['n_training_cells'];np.testing.assert_allclose(yy[use].mean(0),p['mean_programme_scores'],rtol=1e-6,atol=2e-7)
    expected=(yy[use].mean(0)-yy.mean(0))/np.maximum(yy.std(0),1e-12);np.testing.assert_allclose(expected,p['standardized_programme_difference'],rtol=1e-5,atol=2e-6)
   replay.append(dict(partition=key,n_patterns=len(set(labs)),source_NMI=float(normalized_mutual_info_score(sourcecodes,labs)),work_unit_NMI=float(normalized_mutual_info_score(u,labs))))
  for res in m['resolutions']:
   k=['r'+str(res).replace('.','p')+'_seed'+str(seed) for seed in m['seeds']];replay.append(dict(resolution=res,seed_adjusted_Rand=float(adjusted_rand_score(model['reference_labels'][k[0]],model['reference_labels'][k[1]]))))
  reports.append(dict(context_code=c,n_training_cells=n,n_programmes=len(model['program_ids']),n_sources=len(np.unique(sourcecodes)),source_budget_axis_PCA_graph_and_profile_replay_pass=True,exact_float32_Fortran_layout_scaler_original_thread_budget_and_full_SVD_fit_transform_replayed=True,diagnostics=replay));total+=n;storage_gate(R,cfg);save(R/'checkpoints'/('context_'+str(c)+'_audited_v1.json'),dict(n_contexts=len(reports),n_replayed_cells=total))
 op=R/'outputs/accepted_development_state_reference_numeric_audit_v1.json';save(op,dict(schema='TME_development_state_partition_numeric_audit_v1',producer_acceptance=binding(ap),n_contexts=len(reports),n_replayed_reference_cells=total,reports=reports,n_lockbox_arrays_opened=0,all_development_source_qualification_unit_and_cell_budgets_checked=True,all_original_native_UCell_feature_vectors_replayed=True,all_stored_profiles_and_graph_label_axes_checked=True,source_and_seed_diagnostics_descriptive_not_automatic_biological_validation=True,clinical_outcomes_used=False))
 save(R/'provenance/terminal_v1.json',dict(status='FROZEN_STATE_REFERENCE_NUMERIC_AND_SCOPE_AUDIT_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
