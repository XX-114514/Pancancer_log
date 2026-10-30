"""Freeze development-only identity-conditioned multiview programme inputs.

Raw means, source-centred and depth-adjusted views are different estimands.
Missing observation opportunity remains NaN; signed values never enter CoVarNet.
"""
import argparse,datetime,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd

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
def residual(y,covariates,minimum):
 y=np.asarray(y,float);x=np.asarray(covariates,float);assert np.isfinite(y).all() and np.isfinite(x).all()
 if len(y)<minimum:return None,'insufficient_units'
 scale=np.std(x,axis=0);xx=np.column_stack([np.ones(len(y)),np.divide(x-x.mean(0),scale,out=np.zeros_like(x),where=scale>0)])
 rank=np.linalg.matrix_rank(xx)
 if rank!=xx.shape[1]:return None,'collinear_covariates'
 if len(y)<=2*rank:return None,'insufficient_residual_degrees_of_freedom'
 return y-xx@np.linalg.lstsq(xx,y,rcond=None)[0],'estimable'
def transform(values):
 out=np.full(values.shape,np.nan,float);parameters=[]
 for j in range(values.shape[1]):
  good=np.isfinite(values[:,j]);yy=values[good,j];mu=float(yy.mean()) if len(yy) else None;sd=float(yy.std()) if len(yy) else None
  if sd is not None and sd>1e-12:out[good,j]=(yy-mu)/sd
  parameters.append(dict(n=int(good.sum()),mean=mu,standard_deviation=sd,constant=sd is None or sd<=1e-12))
 return out,parameters

def main(R):
 cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);C=Path(cfg['candidates']);assert load(C/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 cp=C/'outputs/accepted_development_state_measurement_candidates_v1.json';assert sha(cp)==load(C/'provenance/terminal_v1.json')['acceptance_SHA256'];cand=load(cp)
 reads=load(bound(cand['readout_acceptance']));app=load(bound(cand['applicability']));roles=load(bound(cand['source_roles']));assert reads['n_lockbox_score_arrays_read']==0
 pc=load(Path(cfg['score_producer'])/'config/spec_v1.json');catalog=load(bound(pc['catalog']));genes={p['program_id']:tuple(sorted(set(p['genes']))) for p in catalog['programs']};names=reads['program_ids'];index={p:i for i,p in enumerate(names)}
 role={x['source']:x['role'] for x in roles['source_roles']};assert all(role[x['source']]=='development' for x in reads['readout_sources'])
 contexts={x['code']:x for x in app['contexts']};loaded={};all_units=[]
 for r in reads['readout_sources']:
  with np.load(bound(r['readouts']),allow_pickle=False) as f:d={k:f[k] for k in ['unit_code','context_code','n_cells','unit_eligible_for_module_analysis','means','mean_log1p_RNA','mean_detected_native_genes']}
  loaded[r['source']]=d;all_units.extend(d['unit_code'].tolist())
 units=np.array(sorted(set(all_units)),np.int32);unitpos={int(u):i for i,u in enumerate(units)};source_axis=np.full(len(units),'',object)
 for source,d in loaded.items():
  for u in np.unique(d['unit_code']):
   i=unitpos[int(u)];assert source_axis[i] in ['',source];source_axis[i]=source
 assert all(source_axis)
 source_labels=sorted(set(source_axis));source_codes=np.array([source_labels.index(s) for s in source_axis],np.int32)
 views=[];held=[];counts=[];details=[]
 for c in sorted({x['context_code'] for x in cand['context_programme_candidates']}):
  candidates=[x for x in cand['context_programme_candidates'] if x['context_code']==c];eligible=[];by_genes={}
  for r in sorted(candidates,key=lambda x:x['program_id']):
   if r['n_supporting_development_sources']<cfg['minimum_sources']:
    held.append(dict(context_code=c,program_id=r['program_id'],reason='source_limited_retained_for_context_specific_analysis'));continue
   key=genes[r['program_id']]
   if key in by_genes:
    held.append(dict(context_code=c,program_id=r['program_id'],reason='identical_external_gene_set_alias',representative=by_genes[key]));continue
   by_genes[key]=r['program_id'];eligible.append(r)
  if not eligible:continue
  pids=[r['program_id'] for r in eligible];raw=np.full((len(units),len(pids)),np.nan,float);adjusted=np.full_like(raw,np.nan);centred=np.full_like(raw,np.nan);ncell=np.zeros(len(units),np.int32);cov=np.full((len(units),3),np.nan,float)
  for source,d in loaded.items():
   ix=np.flatnonzero((d['context_code']==c)&d['unit_eligible_for_module_analysis']&(d['n_cells']>=cfg['minimum_context_cells']))
   if not len(ix):continue
   rows=np.array([unitpos[int(u)] for u in d['unit_code'][ix]],int);ncell[rows]=d['n_cells'][ix]
   cov[rows]=np.column_stack([d['mean_log1p_RNA'][ix],d['mean_detected_native_genes'][ix],np.log1p(d['n_cells'][ix])])
   for j,r in enumerate(eligible):
    if source not in r['supporting_sources']:continue
    yy=d['means'][ix,index[r['program_id']]];assert np.isfinite(yy).all();raw[rows,j]=yy;centred[rows,j]=yy-yy.mean()
    rr,status=residual(yy,cov[rows],cfg['minimum_residual_units'])
    if rr is not None:adjusted[rows,j]=rr
    details.append(dict(source=source,context_code=c,program_id=r['program_id'],n_units=len(rows),depth_residual_status=status,required_covariates=['mean_log1p_RNA','mean_detected_native_genes','log1p_context_cells'],no_collinear_covariate_automatically_removed=True))
  keep=np.array([np.count_nonzero(np.isfinite(raw[:,j]))>=cfg['minimum_feature_units'] and np.nanstd(raw[:,j])>1e-12 for j in range(len(pids))])
  for j in np.flatnonzero(~keep):held.append(dict(context_code=c,program_id=pids[j],reason='insufficient_finite_units_or_constant_readout'))
  if keep.sum()<cfg['minimum_features_per_view']:
   held.append(dict(context_code=c,reason='insufficient_programme_features_for_view'));continue
  pids=[pids[j] for j in np.flatnonzero(keep)];raw=raw[:,keep];centred=centred[:,keep];adjusted=adjusted[:,keep]
  variants={};params={}
  for key,x in [('raw',raw),('source_centered',centred),('source_depth_residual',adjusted)]:variants[key],params[key]=transform(x)
  path=R/'outputs'/('private_context_'+str(c)+'_state_multiview_v1.npz');np.savez_compressed(path,unit_codes=units,source_codes=source_codes,n_context_cells=ncell,raw_original=raw,**variants)
  view=dict(context=contexts[c],program_ids=pids,n_features=len(pids),n_units_with_any_observed=int(np.any(np.isfinite(raw),axis=1).sum()),n_observed_entries=int(np.isfinite(raw).sum()),values=binding(path),training_transform= params,identity_is_frozen_reference_context_not_independent_truth=True)
  views.append(view);counts.append(ncell)
 assert len(views)>=2
 useful=np.zeros(len(units),bool)
 for v in views:
  with np.load(v['values']['path'],allow_pickle=False) as f:useful|=np.any(np.isfinite(f['raw']),axis=1)
 axis=R/'outputs/private_state_multiview_axis_v1.npz';np.savez_compressed(axis,unit_codes=units,source_codes=source_codes,any_observed_view=useful)
 diag=R/'outputs/private_state_depth_residual_status_v1.tsv';pd.DataFrame(details).to_csv(diag,sep='\t',index=False)
 op=R/'outputs/accepted_development_state_multiview_inputs_v1.json';save(op,dict(schema='TME_development_identity_conditioned_state_multiview_inputs_v1',candidate_acceptance=binding(cp),source_roles=cand['source_roles'],applicability=cand['applicability'],axis=binding(axis),source_labels=source_labels,n_units=len(units),n_units_with_observed_views=int(useful.sum()),n_views=len(views),n_features=sum(v['n_features'] for v in views),views=views,held_candidates=held,depth_diagnostics=binding(diag),raw_primary_preserves_source_and_cancer_context=True,source_centered_and_depth_residual_are_sensitivity_estimands=True,missing_unreliable_or_absent_context_remains_NaN=True,signed_residual_not_passed_to_CoVarNet=True,gene_identical_aliases_collapsed_without_effect_selection=True,no_variance_rank_feature_selection=True,n_lockbox_score_arrays_opened=0,clinical_outcomes_used=False,actual_module_model_not_fitted=True,native_universe_sensitivity_still_requires_separate_audit=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='DEVELOPMENT_MULTIVIEW_STATE_INPUTS_FROZEN_MODULE_FIT_PENDING',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
