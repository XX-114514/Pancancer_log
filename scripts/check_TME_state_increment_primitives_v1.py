"""Counterexamples for target leakage, grouped scaling and paired block risk."""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from TME_state_increment_primitives_v1 import project_signed_sender,project_nonnegative_sender,make_predictor_frame,fit_grouped_ridge,paired_block_uncertainty

def main(R):
 rng=np.random.default_rng(42);W=rng.normal(size=(12,3));values=rng.normal(size=(8,12));parents=np.repeat([0,1],6);first=project_signed_sender(values,W,parents,[0],1);change=values.copy();change[:,6:]+=1000;second=project_signed_sender(change,W,parents,[0],1);np.testing.assert_array_equal(first[0],second[0]);missing=values.copy();missing[:,:6]=np.nan;assert np.isnan(project_signed_sender(missing,W,parents,[0],1)[0]).all()
 A=np.array([[1,.1],[.1,1],[.3,.5],[.6,.2]]);frequency=np.array([[.2,.4,.3,.5],[.1,.2,.3,.4],[0,0,.4,.5]]);p=np.array([0,0,1,1]);low=np.array([0,.3,0,0]);span=np.ones(4);h,status,it=project_nonnegative_sender(frequency,A,p,[0],1,low,span);assert status[1]=='below_frozen_training_origin' and np.isnan(h[1]).all();positive=np.array([[1,.2],[.2,1],[.1,.4],[.5,.2]]);hh,ss,ii=project_nonnegative_sender(np.array([[0.,0.,.2,.4]]),positive,p,[0],1,np.zeros(4),np.ones(4));assert ss[0]=='estimable' and (hh[0]==0).all()
 unid,st,_=project_nonnegative_sender(frequency,np.ones((4,2)),p,[0],1,np.zeros(4),np.ones(4));assert (st=='sender_basis_not_identifiable').all() and np.isnan(unid).all()
 n=120;x=rng.normal(size=(n,2));state=rng.normal(size=n);target=2*x[:,0]+3*state+rng.normal(0,.02,n);blocks=np.array(['ARTIFICIAL_BLOCK_'+str(i//2) for i in range(n)]);cat=pd.DataFrame({'source':np.repeat(['ARTIFICIAL_A','ARTIFICIAL_B','ARTIFICIAL_C','ARTIFICIAL_HELD'],30)});train=np.arange(90);held=np.arange(90,120)
 b=make_predictor_frame(x,cat,blocks);e=make_predictor_frame(np.column_stack([x,state]),cat,blocks);base,registry=fit_grouped_ridge(b.iloc[train],target[train],blocks[train],n_jobs=1);enhanced,_=fit_grouped_ridge(e.iloc[train],target[train],blocks[train],n_jobs=1);p0=base.predict(b.iloc[held]).ravel();p1=enhanced.predict(e.iloc[held]).ravel();assert np.mean((target[held]-p1)**2)<.02*np.mean((target[held]-p0)**2);assert all(z['candidate_blocks_disjoint'] for z in registry)
 pre=enhanced.best_estimator_.regressor_.named_steps['preprocess'];feature_names=pre.get_feature_names_out();assert not any('__candidate_block' in n for n in feature_names)
 np.testing.assert_allclose(pre.named_transformers_['numeric'].mean_,np.column_stack([x,state])[train].mean(0));np.testing.assert_allclose(enhanced.best_estimator_.transformer_.mean_,target[train].mean())
 null_target=2*x[:,0]+rng.normal(0,.02,n);nb,_=fit_grouped_ridge(b.iloc[train],null_target[train],blocks[train],n_jobs=1);ne,_=fit_grouped_ridge(e.iloc[train],null_target[train],blocks[train],n_jobs=1);err0=np.mean((null_target[held]-nb.predict(b.iloc[held]).ravel())**2);err1=np.mean((null_target[held]-ne.predict(e.iloc[held]).ravel())**2);assert abs(err0-err1)<.01
 base_loss=np.full((24,3),2.);enh_loss=np.full((24,3),1.);base_loss[:2,1]=np.nan;enh_loss[:2,1]=np.nan;bl=np.repeat(np.arange(12),2).astype(str);a=paired_block_uncertainty(base_loss,enh_loss,bl,repeats=256);assert abs(a['mean_programme_candidate_block_scaled_error_improvement']-1)<1e-12 and a['CI_low']==1 and a['paired_block_prediction_swap_p']<.02
 duplicate=paired_block_uncertainty(np.repeat(base_loss,2,0),np.repeat(enh_loss,2,0),np.repeat(bl,2),repeats=256);assert a==duplicate;null=paired_block_uncertainty(base_loss,base_loss,bl,repeats=256);assert null['paired_block_prediction_swap_p']==1 and null['CI_low']==0
 out=R/'outputs/accepted_state_increment_primitives_synthetic_v1.json';out.write_text(json.dumps(dict(schema='TME_state_increment_primitives_synthetic_v1',all_tests_passed=True,n_tests=9,sender_receiver_target_leakage_checked=True,missing_sender_out_of_origin_and_rank_failure_checked=True,real_sklearn_grouped_ridge_and_training_scaler_executed=True,block_not_a_predictor=True,composition_only_null_and_training_only_XY_scaling_checked=True,paired_risk_missingness_and_duplicate_material_invariance=True,no_human_science_result=True),indent=2));(R/'provenance/terminal_v1.json').write_text(json.dumps(dict(status='STATE_INCREMENT_SENDER_AND_GROUPED_PREDICTION_PRIMITIVES_SYNTHETIC_PASS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=hashlib.sha256(out.read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  (R/'provenance/failure_v1.json').write_text(json.dumps(dict(traceback=traceback.format_exc()),indent=2));raise
