"""Replay official MOFA saved data, factor order, reconstruction and variance.

The HDF5 writer sorts factors; private NPZ keeps training order. The audit
binds both explicitly instead of treating matching factor numbers as identity.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,h5py
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

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
def source_R2(z,groups):
 den=np.square(z-z.mean(0)).sum(0);res=z.copy()
 for g in np.unique(groups):res[groups==g]-=res[groups==g].mean(0)
 return np.divide(den-np.square(res).sum(0),den,out=np.full(z.shape[1],np.nan),where=den>0)
def signed_match(A,B,minimum):
 den=np.outer(np.linalg.norm(A,axis=0),np.linalg.norm(B,axis=0));cos=np.divide(A.T@B,den,out=np.zeros(den.shape),where=den>0);a,b=cos.shape;benefit=np.zeros((a+b,a+b));benefit[:a,:b]=np.abs(cos)-minimum;ri,ci=linear_sum_assignment(-benefit);matches=[]
 for i,j in zip(ri,ci):
  if i<a and j<b and abs(cos[i,j])>=minimum:matches.append((int(i),int(j),float(cos[i,j])))
 return matches,[i for i in range(a) if i not in {x[0] for x in matches}],[j for j in range(b) if j not in {x[1] for x in matches}]

def main(R):
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);start=time.monotonic();T=Path(cfg['math_test_stage'])
 while not (T/'provenance/supervisor_terminal_v1.json').exists():
  assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
 assert load(T/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 test=T/'outputs/accepted_signed_factor_matching_math_v1.json';assert load(test)['all_tests_passed']
 assert sha(T/'scripts/audit_TME_state_MOFA_fits_v1.py')==sha(__file__)
 while not (P/'provenance/supervisor_terminal_v1.json').exists():
  assert time.monotonic()-start<cfg['dependency_wait_seconds'];time.sleep(5)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 ap=P/'outputs/accepted_development_state_MOFA_fits_v1.json';assert sha(ap)==load(P/'provenance/terminal_v1.json')['acceptance_SHA256'];m=load(ap);pcfg=load(P/'config/spec_v1.json');inp=load(bound(m['input_acceptance']))
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);iview={v['context']['code']:v for v in inp['views']};data=[]
 with np.load(bound(inp['axis']),allow_pickle=False) as z:units=z['unit_codes'];source=z['source_codes']
 for view in m['views']:
  c=int(view['name'].removeprefix('context_'));original=iview[c];js=np.array([original['program_ids'].index(n) for n in view['features']]);arrays={}
  with np.load(bound(original['values']),allow_pickle=False) as z:
   mask=np.isfinite(z['raw'])&np.isfinite(z['source_centered'])&np.isfinite(z['source_depth_residual'])
   for variant in pcfg['variants']:arrays[variant]=np.where(mask[:,js],z[variant][:,js],np.nan)
  data.append(dict(view=view,context=c,arrays=arrays,parent=original['context']['compartment']))
 selected=np.zeros(len(units),bool)
 for d in data:selected|=np.isfinite(d['arrays']['raw']).any(1)
 units=units[selected];source=source[selected];processed={}
 for variant in pcfg['variants']:
  processed[variant]=[]
  for d in data:
   x=d['arrays'][variant][selected].astype(float);center=np.nanmean(x,axis=0);x-=center;scale=np.nanstd(x);assert scale>0;x/=scale;processed[variant].append(x)
 models=[];details=[]
 for item in m['models']:
  with np.load(bound(item['arrays']),allow_pickle=False) as z:
   np.testing.assert_array_equal(z['qualification_unit_code'],units);np.testing.assert_array_equal(z['source_codes'],source);Z=z['Z'].copy();W=[z['W_'+str(d['context'])].copy() for d in data]
  assert Z.shape==(len(units),item['rank']) and all(w.shape==(len(d['view']['features']),item['rank']) for w,d in zip(W,data))
  y=processed[item['variant']];per=np.zeros((len(data),item['rank']));total=np.zeros(len(data))
  for j,(x,w) in enumerate(zip(y,W)):
   good=np.isfinite(x);ss=np.square(x[good]).sum();assert ss>0
   total[j]=1-np.square((x-Z@w.T)[good]).sum()/ss
   for k in range(item['rank']):per[j,k]=1-np.square((x-np.outer(Z[:,k],w[:,k]))[good]).sum()/ss
  np.testing.assert_allclose(per,item['view_factor_variance_explained'],rtol=1e-8,atol=1e-10);np.testing.assert_allclose(total,item['view_total_variance_explained'],rtol=1e-8,atol=1e-10)
  with h5py.File(bound(item['model']),'r') as h:
   assert h['samples/development'].asstr()[:].tolist()==['qualification_work_unit_'+str(int(u)) for u in units]
   WW=[]
   for j,d in enumerate(data):
    view=d['view']['name'];feature_names=h['features/'+view].asstr()[:].tolist();assert feature_names==[view+'::'+n for n in d['view']['features']]
    x=h['data/'+view+'/development'][:];np.testing.assert_array_equal(np.isnan(x),np.isnan(y[j]));np.testing.assert_allclose(x,y[j],rtol=1e-12,atol=1e-12,equal_nan=True);WW.append(h['expectations/W/'+view][:].T)
   zz=h['expectations/Z/development'][:].T;full=np.vstack(W);saved=np.vstack(WW);order=[];unused=set(range(item['rank']))
   for k in range(item['rank']):
    possible=[i for i in sorted(unused) if np.array_equal(full[:,i],saved[:,k]) and np.array_equal(Z[:,i],zz[:,k])];assert possible,'Saved factor cannot be paired to original W/Z';i=possible[0];unused.remove(i);order.append(i)
   assert not unused;np.testing.assert_array_equal(full[:,order],saved);np.testing.assert_array_equal(Z[:,order],zz)
   np.testing.assert_allclose(h['variance_explained/r2_per_factor/development'][:]/100,per[:,order],rtol=1e-8,atol=1e-10);np.testing.assert_allclose(h['variance_explained/r2_total/development'][:]/100,total,rtol=1e-8,atol=1e-10)
   elbo=h['training_stats/elbo'][:];finite=elbo[np.isfinite(elbo)];assert len(finite)>=2 and finite[-1]>=finite[0]-1e-7*(1+abs(finite[0]));iterations=len(h['training_stats/time'])
  sr=source_R2(Z,source);factor_scope=[]
  for k in range(item['rank']):
   active=[d['parent'] for j,d in enumerate(data) if per[j,k]>=cfg['minimum_view_variance_fraction']]
   factor_scope.append(dict(original_factor_index=k+1,n_active_contexts=len(active),active_broad_parents=sorted(set(active)),n_active_broad_parents=len(set(active)),source_group_R2=float(sr[k]) if np.isfinite(sr[k]) else None,view_variance_fraction_threshold=cfg['minimum_view_variance_fraction'],not_automatic_functional_or_communication_validation=True))
  detail=dict(variant=item['variant'],rank=item['rank'],seed=item['seed'],n_units=len(units),n_iterations_recorded=iterations,at_configured_iteration_cap=iterations>=pcfg['iterations'],HDF5_sorted_factor_to_original_index=[i+1 for i in order],all_saved_data_masks_features_and_WZ_permutation_replayed=True,independent_variance_explained_pass=True,ELBO_improvement_verified=True,factor_scope=factor_scope)
  details.append(detail);models.append(dict(item=item,W=np.vstack(W),Z=Z));storage_gate(R,cfg);save(R/'checkpoints'/('model_'+str(len(details))+'_audited_v1.json'),dict(n_models=len(details)))
 comparisons=[]
 for i,left in enumerate(models):
  for right in models[i+1:]:
   a=left['item'];b=right['item']
   if a['variant']!=b['variant'] and (a['rank']!=b['rank'] or a['seed']!=b['seed']):continue
   if a['variant']==b['variant'] and a['rank']!=b['rank'] and a['seed']!=b['seed']:continue
   matches,ua,ub=signed_match(left['W'],right['W'],cfg['minimum_cosine']);rows=[]
   for l,r,cos in matches:
    x=left['Z'][:,l];yy=right['Z'][:,r]*(1 if cos>=0 else -1);rho=float(spearmanr(x,yy).statistic) if np.ptp(x)>0 and np.ptp(yy)>0 else None
    rows.append(dict(left_factor=l+1,right_factor=r+1,loading_cosine_signed=cos,sign_corrected_activity_Spearman=rho))
   comparisons.append(dict(left={k:a[k] for k in ['variant','rank','seed']},right={k:b[k] for k in ['variant','rank','seed']},matches=rows,unmatched_left=[x+1 for x in ua],unmatched_right=[x+1 for x in ub]))
 op=R/'outputs/accepted_state_MOFA_numeric_comparison_v1.json';save(op,dict(schema='TME_state_MOFA_numeric_replay_and_signed_factor_comparison_v1',producer_acceptance=binding(ap),signed_matching_math_acceptance=binding(test),n_models=len(models),n_units=len(units),n_views=m['n_views'],n_features=m['n_features'],models=details,comparisons=comparisons,matching='absolute_cosine_one_to_one_Hungarian_with_unmatched_and_activity_sign_correction',minimum_cosine=cfg['minimum_cosine'],all_training_data_missing_masks_and_axes_replayed=True,no_lockbox_arrays_opened=True,clinical_outcomes_used=False,source_adjustment_changes_estimand_not_automatic_failure=True,source_holdout_and_composition_increment_validation_pending=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='OFFICIAL_MOFA_MASK_FACTOR_ORDER_VARIANCE_AND_COMPARISON_AUDIT_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
