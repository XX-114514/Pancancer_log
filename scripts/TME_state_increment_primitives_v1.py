"""Standard grouped-ridge and sender-only state increment primitives.

Candidate blocks are not newly certified patients. Signed-basis coordinates
are a regularized projection, not exact MOFA posterior inference.
"""
import inspect
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.compose import ColumnTransformer,TransformedTargetRegressor
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold,GridSearchCV

def project_signed_sender(values,loadings,feature_parents,sender_parents,receiver_parent,alpha=1.):
 values=np.asarray(values,float);loadings=np.asarray(loadings,float);parent=np.asarray(feature_parents)
 assert values.ndim==2 and loadings.shape[0]==values.shape[1] and parent.shape==(values.shape[1],)
 assert np.isfinite(loadings).all() and alpha>0
 use=np.isin(parent,sender_parents)&(parent!=receiver_parent);observed=np.isfinite(values)&use[None,:];patterns,inverse=np.unique(observed,axis=0,return_inverse=True);k=loadings.shape[1];coordinates=np.full((len(values),k),np.nan);nobs=observed.sum(1);rank=np.zeros(len(values),int)
 for j,mask in enumerate(patterns):
  rows=np.flatnonzero(inverse==j)
  if mask.sum()<k or not np.any(loadings[mask]):continue
  x=loadings[mask];rank[rows]=np.linalg.matrix_rank(x)
  regression=Ridge(alpha=alpha,fit_intercept=False,solver='svd').fit(x,values[np.ix_(rows,mask)].T)
  coordinates[rows]=regression.coef_
 return coordinates,nobs,rank

def project_nonnegative_sender(frequency,effective_basis,node_parents,sender_parents,receiver_parent,minimum,span,maximum_condition=1e8,max_iter=2000):
 frequency=np.asarray(frequency,float);A=np.asarray(effective_basis,float);parent=np.asarray(node_parents);lo=np.asarray(minimum,float);scale=np.asarray(span,float)
 assert frequency.ndim==2 and A.shape[0]==frequency.shape[1] and parent.shape==lo.shape==scale.shape==(A.shape[0],) and (scale>0).all() and np.isfinite(A).all() and (A>=0).all()
 keep=np.isin(parent,sender_parents)&(parent!=receiver_parent);a=A[keep];k=A.shape[1];output=np.full((len(frequency),k),np.nan);status=np.full(len(frequency),'missing_sender_observation',object);iterations=np.zeros(len(frequency),int)
 if a.shape[0]<k or np.linalg.matrix_rank(a)<k or np.linalg.cond(a)>maximum_condition:return output,np.full(len(frequency),'sender_basis_not_identifiable',object),iterations
 x=(frequency[:,keep]-lo[keep])/scale[keep];finite=np.isfinite(x).all(1);below=finite&(x<0).any(1);status[below]='below_frozen_training_origin';valid=finite&~below
 if not valid.any():return output,status,iterations
 target=x[valid].T;H=np.full((k,target.shape[1]),1/k);den=np.maximum(a.sum(0),1e-12)[:,None];last=np.full(target.shape[1],np.inf);converged=np.zeros(target.shape[1],bool);it=np.zeros(target.shape[1],int)
 for step in range(1,max_iter+1):
  active=~converged
  if not active.any():break
  hh=H[:,active];xx=target[:,active];hh*=((a.T@(xx/np.maximum(a@hh,1e-12)))/den);H[:,active]=hh
  if step%20==0:
   yy=np.maximum(a@hh,1e-12);loss=(yy-xx).sum(0);term=np.zeros_like(xx);positive=xx>0;term[positive]=xx[positive]*np.log(xx[positive]/yy[positive]);loss+=term.sum(0);ids=np.flatnonzero(active);done=np.isfinite(last[ids])&(np.abs(last[ids]-loss)<1e-8*(1+np.abs(last[ids])));converged[ids[done]]=True;it[ids[done]]=step;last[ids]=loss
 remaining=~converged;it[remaining]=max_iter;output[valid]=H.T;iterations[valid]=it;status[valid]='estimable';status[np.flatnonzero(valid)[remaining]]='projection_iteration_cap'
 return output,status,iterations

def candidate_block_weights(blocks):
 b=np.asarray(blocks).astype(str);u,inverse,counts=np.unique(b,return_inverse=True,return_counts=True);w=1/counts[inverse];return w/w.mean()

def make_predictor_frame(numeric,categorical,blocks):
 numeric=np.asarray(numeric,float);assert numeric.ndim==2 and np.isfinite(numeric).all();c=pd.DataFrame(categorical).reset_index(drop=True).fillna('unknown').astype(str)
 assert len(c)==len(numeric)==len(blocks);n=pd.DataFrame(numeric,columns=['numeric_'+str(i) for i in range(numeric.shape[1])]);c.columns=['category_'+str(i) for i in range(c.shape[1])];out=pd.concat([n,c],axis=1);out['__candidate_block']=np.asarray(blocks).astype(str);return out

def normalized_block_MSE(estimator,X,y):
 y=np.asarray(y,float);y=y[:,None] if y.ndim==1 else y;p=np.asarray(estimator.predict(X));p=p[:,None] if p.ndim==1 else p;w=candidate_block_weights(X['__candidate_block']);e=((y-p)/estimator.transformer_.scale_)**2
 return -float(np.sum(w[:,None]*e)/(w.sum()*e.shape[1]))

def fit_grouped_ridge(frame,target,blocks,alphas=(.01,.1,1.,10.,100.,1000.),n_splits=5,n_jobs=1):
 y=np.asarray(target,float);y=y[:,None] if y.ndim==1 else y;blocks=np.asarray(blocks).astype(str);assert len(frame)==len(y)==len(blocks) and np.isfinite(y).all();assert np.array_equal(frame['__candidate_block'].to_numpy(),blocks)
 numeric=[c for c in frame if c.startswith('numeric_')];categorical=[c for c in frame if c.startswith('category_')];ohe_opts={'handle_unknown':'ignore'};ohe_opts['sparse_output' if 'sparse_output' in inspect.signature(OneHotEncoder).parameters else 'sparse']=False
 pre=ColumnTransformer([('numeric',StandardScaler(),numeric),('categorical',OneHotEncoder(**ohe_opts),categorical)],remainder='drop',sparse_threshold=0)
 estimator=TransformedTargetRegressor(regressor=Pipeline([('preprocess',pre),('ridge',Ridge(solver='cholesky'))]),transformer=StandardScaler())
 folds=list(GroupKFold(n_splits=min(n_splits,len(set(blocks)))).split(frame,y,blocks));registry=[]
 for tr,va in folds:
  assert not set(blocks[tr])&set(blocks[va]);registry.append(dict(n_train_units=len(tr),n_validation_units=len(va),n_train_blocks=len(set(blocks[tr])),n_validation_blocks=len(set(blocks[va])),candidate_blocks_disjoint=True))
 search=GridSearchCV(estimator,{'regressor__ridge__alpha':list(alphas)},scoring=normalized_block_MSE,cv=folds,refit=True,error_score='raise',n_jobs=n_jobs,pre_dispatch=max(1,n_jobs))
 search.fit(frame,y,ridge__sample_weight=candidate_block_weights(blocks))
 return search,registry

def block_loss_matrix(loss,blocks):
 loss=np.asarray(loss,float);loss=loss[:,None] if loss.ndim==1 else loss;b=np.asarray(blocks).astype(str);unique,inverse=np.unique(b,return_inverse=True);finite=np.isfinite(loss);total=np.zeros((len(unique),loss.shape[1]));count=np.zeros_like(total)
 np.add.at(total,inverse,np.where(finite,loss,0));np.add.at(count,inverse,finite);means=np.divide(total,count,out=np.full(total.shape,np.nan),where=count>0);return unique,means

def paired_block_uncertainty(baseline_loss,enhanced_loss,blocks,repeats=2000,seed=20261031):
 b,base=block_loss_matrix(baseline_loss,blocks);bb,enhanced=block_loss_matrix(enhanced_loss,blocks);assert np.array_equal(b,bb);assert np.array_equal(np.isfinite(base),np.isfinite(enhanced));delta=base-enhanced;valid=np.isfinite(delta);n=valid.sum(0);assert (n>0).all();point=float(np.mean(np.nansum(delta,0)/n));rng=np.random.default_rng(seed);draw=rng.integers(0,len(b),(repeats,len(b)));multiplicity=np.zeros((repeats,len(b)),int);np.add.at(multiplicity,(np.arange(repeats)[:,None],draw),1);den=multiplicity@valid.astype(float);num=multiplicity@np.where(valid,delta,0);boot=np.divide(num,den,out=np.full(num.shape,np.nan),where=den>0);bs=np.nanmean(boot,1);finite=bs[np.isfinite(bs)]
 signs=rng.choice([-1.,1.],size=(repeats,len(b)));permutation=np.mean((signs@np.where(valid,delta,0))/n,axis=1);p=(1+np.count_nonzero(np.abs(permutation)>=abs(point)))/(repeats+1)
 return dict(mean_programme_candidate_block_scaled_error_improvement=point,CI_low=float(np.quantile(finite,.025)),CI_high=float(np.quantile(finite,.975)),n_candidate_blocks=len(b),n_finite_bootstraps=len(finite),paired_block_prediction_swap_p=float(p),prediction_swap_requires_exchangeable_model_labels_within_block=True,not_an_independent_patient_or_clinical_validation=True)
