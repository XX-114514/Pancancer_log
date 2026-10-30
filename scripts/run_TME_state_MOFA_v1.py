"""Official MOFA fits for signed, masked identity-conditioned state readouts.

Models are development-only exploratory coordination patterns. They are not
cell types, communications or clinical mechanisms. Missing inputs stay missing.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def binding(p):return dict(path=str(p),SHA256=sha(p))
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def wait(root,seconds,start):
 p=Path(root)/'provenance/supervisor_terminal_v1.json'
 while not p.exists():
  assert time.monotonic()-start<seconds;time.sleep(5)
 assert load(p)['actual_child_exit_code']==0 and load(p).get('science_terminal_exists',True) and not (Path(root)/'provenance/failure_v1.json').exists(),'Dependency requires scientific acceptance, not shell exit zero alone'
def fit_mofa(data,names,features,samples,cfg,out):
 from mofapy2.run.entry_point import entry_point
 entry=entry_point();entry.set_data_options(scale_views=True,scale_groups=False,center_groups=True,use_float32=False)
 entry.set_data_matrix([[np.asarray(x,float).copy()] for x in data],likelihoods=['gaussian']*len(data),views_names=names,groups_names=['development'],samples_names=[samples],features_names=[[view+'::'+feature for feature in ff] for view,ff in zip(names,features)])
 entry.set_model_options(factors=cfg['factors'],spikeslab_weights=True,ard_weights=True,spikeslab_factors=False,ard_factors=False)
 entry.set_train_options(iter=cfg['iterations'],convergence_mode='medium',seed=cfg['seed'],quiet=True,verbose=False,gpu_mode=False,weight_views=False,dropR2=None,save_interrupted=False)
 entry.build();entry.run();entry.save(str(out),save_data=True,save_parameters=False)
 ex=entry.model.getExpectations(True,'Z','W');Z=np.asarray(ex['Z']);W=[np.asarray(x) for x in ex['W']]
 assert Z.shape==(len(samples),cfg['factors']) and all(w.shape==(len(f),cfg['factors']) for w,f in zip(W,features));assert np.isfinite(Z).all() and all(np.isfinite(w).all() for w in W)
 return entry,Z,W

def main(R):
 cfg=load(R/'config/spec_v1.json');start=time.monotonic()
 for root in [cfg['inputs'],cfg['vendor_stage'],cfg['math_test_stage'],cfg['sensitivity']]:wait(root,cfg['dependency_wait_seconds'],start)
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);V=Path(cfg['vendor_stage']);vp=V/'outputs/accepted_official_MOFA_vendor_v1.json';assert sha(vp)==load(V/'provenance/terminal_v1.json')['acceptance_SHA256'];vendor=load(vp);sys.path.insert(0,vendor['vendor_path'])
 assert vendor['version']==cfg['expected_version'];import mofapy2;assert Path(mofapy2.__file__).resolve().is_relative_to(Path(vendor['vendor_path']).resolve())
 T=Path(cfg['math_test_stage']);tp=T/'outputs/accepted_MOFA_math_and_API_synthetic_v1.json';assert sha(tp)==load(T/'provenance/terminal_v1.json')['acceptance_SHA256'];assert load(tp)['all_tests_passed']
 I=Path(cfg['inputs']);ip=I/'outputs/accepted_development_state_multiview_inputs_v1.json';assert sha(ip)==load(I/'provenance/terminal_v1.json')['acceptance_SHA256'];m=load(ip)
 S=Path(cfg['sensitivity']);sp=S/'outputs/accepted_development_state_native_sensitivity_v1.json';assert sha(sp)==load(S/'provenance/terminal_v1.json')['acceptance_SHA256'];assert load(sp)['n_lockbox_arrays_opened']==0
 with np.load(bound(m['axis']),allow_pickle=False) as z:unit=z['unit_codes'];source=z['source_codes']
 viewdata=[]
 for view in m['views']:
  with np.load(bound(view['values']),allow_pickle=False) as z:d={k:z[k].copy() for k in ['raw','source_centered','source_depth_residual']}
  mask=np.isfinite(d['raw'])&np.isfinite(d['source_centered'])&np.isfinite(d['source_depth_residual'])
  cols=np.flatnonzero(mask.sum(0)>=cfg['minimum_feature_units'])
  if len(cols)<cfg['minimum_features_per_view']:continue
  # Equal observed entries and feature units for each paired estimand.
  d={k:np.where(mask[:,cols],v[:,cols],np.nan) for k,v in d.items()};viewdata.append(dict(context_code=view['context']['code'],data=d,features=[view['program_ids'][j] for j in cols]))
 assert len(viewdata)>=2
 seen=np.zeros(len(unit),bool)
 for v in viewdata:seen|=np.any(np.isfinite(v['data']['raw']),axis=1)
 units=unit[seen];source=source[seen];samples=['qualification_work_unit_'+str(int(i)) for i in units];names=['context_'+str(v['context_code']) for v in viewdata];features=[v['features'] for v in viewdata]
 assert len(samples)>15
 models=[]
 for variant in cfg['variants']:
  values=[v['data'][variant][seen] for v in viewdata]
  for rank in cfg['ranks']:
   for seed in cfg['seeds']:
    label=variant+'_K'+str(rank)+'_seed'+str(seed);model_path=R/'outputs'/('private_'+label+'.hdf5');request=dict(factors=rank,iterations=cfg['iterations'],seed=seed)
    entry,Z,W=fit_mofa(values,names,features,samples,request,model_path)
    arrays=R/'outputs'/('private_'+label+'_factors_v1.npz');np.savez_compressed(arrays,qualification_unit_code=units,source_codes=source,Z=Z,**{'W_'+str(v['context_code']):w for v,w in zip(viewdata,W)})
    per=np.asarray(entry.model.calculate_variance_explained()[0]);total=np.asarray(entry.model.calculate_variance_explained(total=True)[0]);assert per.shape==(len(values),rank) and np.isfinite(per).all() and np.isfinite(total).all()
    item=dict(variant=variant,rank=rank,seed=seed,model=binding(model_path),arrays=binding(arrays),view_factor_variance_explained=per.tolist(),view_total_variance_explained=total.tolist());models.append(item)
    save(R/'checkpoints'/('fit_'+label+'_complete_v1.json'),dict(n_completed=len(models),n_planned=len(cfg['variants'])*len(cfg['ranks'])*len(cfg['seeds'])));storage_gate(R,cfg)
 op=R/'outputs/accepted_development_state_MOFA_fits_v1.json';save(op,dict(schema='TME_development_identity_conditioned_MOFA_fits_v1',input_acceptance=binding(ip),software_acceptance=binding(vp),math_API_acceptance=binding(tp),native_sensitivity_acceptance=binding(sp),n_units=len(samples),n_views=len(viewdata),n_features=sum(map(len,features)),views=[dict(name=n,features=f) for n,f in zip(names,features)],models=models,primary_rank=cfg['primary_rank'],primary_estimand='raw_on_paired_observation_domain',sensitivities=['source_centered','source_depth_residual'],all_variants_share_identical_observation_mask=True,missing_values_passed_as_NaN_to_official_model=True,signed_residuals_not_passed_to_CoVarNet=True,n_lockbox_arrays_opened=0,clinical_outcomes_used=False,independent_patient_or_end_to_end_external_validation_not_claimed=True,actual_cell_state_partition_and_incremental_validation_pending=True))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='OFFICIAL_MOFA_CONTINUOUS_STATE_COORDINATION_FITS_COMPLETE_VALIDATION_PENDING',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
