"""Official single-view RF backend; custom outer blocked evaluation and fusion."""
import hashlib
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def seed(*parts):
    return int.from_bytes(hashlib.sha256(('20261020_native_spatial_MV|'+ '|'.join(map(str,parts))).encode()).digest()[:8],'little')%(2**32-1)


def predict_views(train_base,train_near,train_far,y_train,test_base,test_near,test_far,
                  random_seed,trees=64,leaf=5):
    """No heldout outcome argument; both comparisons have three calibrated RFs."""
    from liana.method.sp import RandomForestModel
    assert all(np.isfinite(x).all() for x in [train_base,train_near,train_far,y_train,test_base,test_near,test_far])
    assert len(y_train)==len(train_base) and len(y_train)>=90
    assert np.var(y_train)>1e-12
    baseline_oob=[];baseline_test=[];full_oob=[];full_test=[];importance=[]
    for j in range(3):
        model=RandomForestModel(seed(random_seed,j),n_estimators=trees,n_jobs=1,max_features='sqrt',min_samples_leaf=leaf)
        names=['B%03d'%k for k in range(train_base.shape[1])]
        model.fit(y_train,train_base,names)
        assert np.isfinite(model.predictions).all()
        baseline_oob.append(model.predictions.copy());baseline_test.append(model.model.predict(test_base))
        if j==0:full_oob.append(model.predictions.copy());full_test.append(model.model.predict(test_base))
    for j,(tr,te) in enumerate([(train_near,test_near),(train_far,test_far)],start=1):
        X=np.column_stack([train_base,tr]);T=np.column_stack([test_base,te])
        names=['B%03d'%k for k in range(train_base.shape[1])]+['S%03d'%k for k in range(tr.shape[1])]
        model=RandomForestModel(seed(random_seed,j),n_estimators=trees,n_jobs=1,max_features='sqrt',min_samples_leaf=leaf)
        model.fit(y_train,X,names)
        assert np.isfinite(model.predictions).all()
        full_oob.append(model.predictions.copy());full_test.append(model.model.predict(T))
        importance.append(np.asarray(model.model.feature_importances_[train_base.shape[1]:],float))
    base=make_pipeline(StandardScaler(),Ridge(alpha=1.,fit_intercept=True))
    full=make_pipeline(StandardScaler(),Ridge(alpha=1.,fit_intercept=True))
    base.fit(np.column_stack(baseline_oob),y_train)
    full.fit(np.column_stack(full_oob),y_train)
    result=(base.predict(np.column_stack(baseline_test)),full.predict(np.column_stack(full_test)))
    assert all(np.isfinite(x).all() for x in result)
    return result[0],result[1],importance


def field_metrics(y,base,full,field):
    """Equal field weight; all returned values are aggregate moments/errors."""
    records=[]
    for f in sorted(set(field)):
        at=field==f;v=y[at]
        records.append(dict(field=str(f),cells=int(at.sum()),mean_y=float(v.mean()),mean_y2=float(np.mean(v*v)),
              mean_squared_error_base=float(np.mean((v-base[at])**2)),mean_squared_error_full=float(np.mean((v-full[at])**2))))
    return records


def aggregate(records):
    if not records:return dict(status='NO_VALID_TEST_FIELDS')
    mean=np.mean([r['mean_y'] for r in records]);var=np.mean([r['mean_y2'] for r in records])-mean**2
    if var<=1e-12:return dict(status='CONSTANT_HELDOUT_TARGET',test_fields=len(records))
    b=np.mean([r['mean_squared_error_base'] for r in records]);m=np.mean([r['mean_squared_error_full'] for r in records])
    return dict(status='ESTIMABLE_BLOCKED_PREDICTION',test_fields=len(records),test_cells=sum(r['cells'] for r in records),
                equal_field_target_variance=float(var),baseline_R2=float(1-b/var),multiview_R2=float(1-m/var),gain_R2=float((b-m)/var),
                equal_field_baseline_MSE=float(b),equal_field_multiview_MSE=float(m))


def fixture():
    rng=np.random.default_rng(2026102005);n=120
    B=rng.normal(size=(n,4));N=rng.normal(size=(n,3));F=rng.normal(size=(n,3));y=.8*B[:,0]+N[:,0]+rng.normal(size=n)*.2
    a,b,imp=predict_views(B[:90],N[:90],F[:90],y[:90],B[90:],N[90:],F[90:],1321,trees=32,leaf=5)
    a2,b2,_=predict_views(B[:90],N[:90],F[:90],y[:90],B[90:],N[90:],F[90:],1321,trees=32,leaf=5)
    np.testing.assert_array_equal(a,a2);np.testing.assert_array_equal(b,b2)
    records=field_metrics(y[90:],a,b,np.repeat(['A','B','C'],10));r=aggregate(records)
    assert r['status']=='ESTIMABLE_BLOCKED_PREDICTION' and np.isfinite(r['gain_R2'])
    altered=[dict(x,mean_squared_error_full=x['mean_squared_error_base']) for x in records]
    assert aggregate(altered)['gain_R2']==0
    assert len(imp)==2 and all(x.shape==(3,) for x in imp)
    return dict(status='PASS_OFFICIAL_RF_CUSTOM_TRAIN_ONLY_FUSION_DETERMINISM_AND_EQUAL_FIELD_ERROR_AGGREGATION',
                no_heldout_outcome_argument_in_fit_function=True,capacity_matched_three_RF_comparisons=True,synthetic_only=True)
