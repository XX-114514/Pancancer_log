"""Independent aggregate reconstruction and preselected direct sklearn refit."""
import hashlib,importlib.util,json,os,sys
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestRegressor


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def seed(*parts):
    return int.from_bytes(hashlib.sha256(('20261020_native_spatial_MV|'+'|'.join(map(str,parts))).encode()).digest()[:8],'little')%(2**32-1)


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def aggregate_direct(table):
    # Separate pandas groupby/moment implementation; never call production aggregate.
    assert table.field.is_unique and np.isfinite(table[['mean_y','mean_y2','mean_squared_error_base','mean_squared_error_full']]).all().all()
    mean=float(table.mean_y.sum()/len(table));var=float(table.mean_y2.sum()/len(table)-mean*mean)
    if var<=1e-12:return dict(status='CONSTANT_HELDOUT_TARGET',test_fields=len(table))
    b=float(table.mean_squared_error_base.sum()/len(table));f=float(table.mean_squared_error_full.sum()/len(table))
    return dict(status='ESTIMABLE_BLOCKED_PREDICTION',test_fields=len(table),test_cells=int(table.cells.sum()),equal_field_target_variance=var,
                equal_field_baseline_MSE=b,equal_field_multiview_MSE=f,baseline_R2=1-b/var,multiview_R2=1-f/var,gain_R2=(b-f)/var)


def close(expected,row):
    assert row['status']==expected['status']
    error=0.
    for k,v in expected.items():
        if k=='status':continue
        np.testing.assert_allclose(float(row[k]),float(v),atol=2e-10,rtol=1e-10)
        error=max(error,abs(float(row[k])-float(v)))
    return error


def direct_fusion(oob,test,y):
    # Explicit centered ridge normal equations; no production pipeline or Ridge call.
    X=np.column_stack(oob);T=np.column_stack(test);mu=X.mean(0);scale=X.std(0);scale[scale==0]=1.
    X=(X-mu)/scale;T=(T-mu)/scale;xm=X.mean(0);ym=y.mean();Xc=X-xm
    beta=np.linalg.solve(Xc.T@Xc+np.eye(X.shape[1]),Xc.T@(y-ym))
    return (T-xm)@beta+ym


def refit_fixed_case(control,plan,contract,output):
    # Before observing effects: image ordinal 0, first frozen target, outer fold 0,
    # both raw and background readouts. Reuse already-audited input transformation;
    # independently rebuild geometry split/buffer, RF estimator and fusion algebra.
    sys.path.insert(0,str(control));prod=load(control/'run_native_spatial_blocked_multiview_v1.py','production_input_transform')
    graph=load(control/'author_spatial_inflow_math_v1.py','audited_graph').graph
    cap=pd.read_csv('/dev/shm/pancancer-native-state-capture-20261020-v2/private_output/private_all_field_proxy_capture_and_quality.tsv.gz',sep='\t',dtype={'sample_id':str,'FOV':str})
    target=contract['targets'][0];fields,scope,_=prod.frames_for_image(target,cap,contract,plan,graph)
    assert len(fields)>=6
    assignments=KMeans(n_clusters=3,n_init=10,random_state=seed(0,'FOV_GEOMETRY')).fit_predict(np.stack([x['center'] for x in fields]))
    guard=max(x['wide_radius'] for x in fields);test=[f for i,f in enumerate(fields) if assignments[i]==0]
    tree=cKDTree(np.concatenate([f['xy'] for f in test]));m=target['proxy_mappings'][0];receiver=m['major']
    candidates=[]
    for i,f in enumerate(fields):
        if assignments[i]==0:continue
        rows=np.flatnonzero(f['pop'][receiver]);d=tree.query(f['xy'][rows])[0];eligible=rows[d>guard]
        if len(eligible)>=30:candidates.append((f,eligible))
    assert len(candidates)>=3
    per_field=min(80,min(len(x) for _,x in candidates));rng=np.random.default_rng(seed(0,receiver,0,'CELL_SELECTION'))
    train=[(f,np.sort(rng.choice(rows,per_field,replace=False))) for f,rows in candidates]
    assert min(float(tree.query(f['xy'][rows])[0].min()) for f,rows in train)>guard
    train_indices=np.concatenate([f['native_rows'][r] for f,r in train]);train_hash=hashlib.sha256(train_indices.astype('<i8').tobytes()).hexdigest()
    other=[i for i,x in enumerate(target['proxy_mappings']) if x['major']==receiver and i!=0]
    checks=[]
    for readout in ['RAW','NEGATIVE_EXCESS']:
        def arrays(parts):
            B=[];N=[];F=[];Y=[];V=[];I=[]
            for f,rows in parts:
                s=f['readouts'][readout]['scores'];near,far=f['readouts'][readout]['views'][receiver]
                B.append(np.column_stack([f['backgrounds'][receiver][rows],s[rows][:,other]]));N.append(near[rows]);F.append(far[rows]);Y.append(s[rows,0]);V.extend([f['FOV']]*len(rows));I.extend(f['native_rows'][rows])
            return np.vstack(B),np.vstack(N),np.vstack(F),np.concatenate(Y),np.asarray(V),np.asarray(I)
        tr=arrays(train);te=arrays([(f,np.flatnonzero(f['pop'][receiver])) for f in test]);assert not np.intersect1d(tr[5],te[5]).size
        base_oob=[];base_test=[];full_oob=[];full_test=[]
        for j in range(3):
            rf=RandomForestRegressor(n_estimators=64,n_jobs=1,max_features='sqrt',min_samples_leaf=5,oob_score=True,random_state=seed(seed(0,receiver,0,'RF'),j))
            rf.fit(tr[0],tr[3]);base_oob.append(rf.oob_prediction_);base_test.append(rf.predict(te[0]))
            if j==0:full_oob.append(rf.oob_prediction_);full_test.append(rf.predict(te[0]))
        for j in [1,2]:
            rf=RandomForestRegressor(n_estimators=64,n_jobs=1,max_features='sqrt',min_samples_leaf=5,oob_score=True,random_state=seed(seed(0,receiver,0,'RF'),j))
            rf.fit(np.column_stack([tr[0],tr[j]]),tr[3]);full_oob.append(rf.oob_prediction_);full_test.append(rf.predict(np.column_stack([te[0],te[j]])))
        b=direct_fusion(base_oob,base_test,tr[3]);f=direct_fusion(full_oob,full_test,tr[3])
        saved=pd.read_csv(output/'private_image_00_test_fields.tsv.gz',sep='\t',dtype={'field':str})
        saved=saved[saved.target_proxy_group_id.eq(m['proxy_group_id'])&saved.readout.eq(readout)&saved.fold.eq(0)].set_index('field')
        max_error=0.
        for fv in sorted(set(te[4])):
            at=te[4]==fv;y=te[3][at];got=[np.mean(y),np.mean(y*y),np.mean((y-b[at])**2),np.mean((y-f[at])**2)]
            expected=saved.loc[fv,['mean_y','mean_y2','mean_squared_error_base','mean_squared_error_full']].to_numpy(float)
            np.testing.assert_allclose(got,expected,atol=2e-11,rtol=1e-10);max_error=max(max_error,float(np.max(np.abs(np.array(got)-expected))))
        folds=pd.read_csv(output/'private_image_00_folds.tsv.gz',sep='\t');row=folds[folds.target_proxy_group_id.eq(m['proxy_group_id'])&folds.readout.eq(readout)&folds.fold.eq(0)].iloc[0]
        assert row.train_selection_sha256==train_hash and row.train_cells==len(tr[3]) and row.test_cells==len(te[3])
        checks.append(dict(readout=readout,status='PRESELECTED_DIRECT_SKLEARN_AND_EXPLICIT_RIDGE_REFIT_ACCEPTED',field_errors_rechecked=len(saved),maximum_aggregate_error=max_error,buffer_and_native_row_disjointness=True))
    return checks


def main():
    control=Path(os.environ['PANCANCER_MV_CONTROL']);out=Path(os.environ['PANCANCER_AUDIT_OUTPUT']);prodout=control/'private_output'
    assert not out.exists() and out.parent.parent==Path('/dev/shm');out.mkdir(mode=0o700)
    plan=json.loads((control/'plan.json').read_text());term=json.loads((prodout/'private_terminal.json').read_text())
    assert term['status']=='COMPLETE_ALL126_REAL_BLOCKED_MULTIVIEW_CASES_PENDING_AUDIT' and term['images_completed']==4
    assert sha(control/'plan.json')==term['plan_sha256'] and sha(control/'run_native_spatial_blocked_multiview_v1.py')==plan['producer_sha256']
    for name,entry in term['files'].items():assert sha(prodout/name)==entry['sha256']
    for name,digest in plan['helper_sha256'].items():assert sha(control/name)==digest
    assert sha(plan['private_contract'])==plan['private_contract_sha256'];contract=json.loads(Path(plan['private_contract']).read_text())
    family=json.loads((prodout/'private_frozen_full_family_before_effects.json').read_text());assert family['plan_sha256']==sha(control/'plan.json')
    expected={(i,m['proxy_group_id'],r) for i,t in enumerate(contract['targets']) for m in t['proxy_mappings'] for r in ['RAW','NEGATIVE_EXCESS']}
    assert len(expected)==126 and expected=={(x['image_ordinal'],x['target_proxy_group_id'],x['readout']) for x in family['family']}
    model=pd.read_csv(prodout/'private_complete126_blocked_model_cases.tsv',sep='\t')
    assert len(model)==126 and set(zip(model.image_ordinal,model.target_proxy_group_id,model.readout))==expected
    max_error=0.;field_positions=0;fold_positions=0;unestimable=0;structural_checks=[]
    for i in range(4):
        folds=pd.read_csv(prodout/('private_image_%02d_folds.tsv.gz'%i),sep='\t')
        fields=pd.read_csv(prodout/('private_image_%02d_test_fields.tsv.gz'%i),sep='\t',dtype={'field':str})
        scope=pd.read_csv(prodout/('private_image_%02d_field_scope.tsv.gz'%i),sep='\t',dtype={'FOV':str});assert scope.FOV.is_unique
        included=set(scope.loc[scope.status.eq('FIELD_INCLUDED'),'FOV']);submodel=model[model.image_ordinal.eq(i)]
        if len(folds) and 'train_selection_sha256' in folds:
            for _,block in folds.dropna(subset=['train_selection_sha256']).groupby(['response_lineage','fold']):assert block.train_selection_sha256.nunique()==1
        for row in submodel.to_dict('records'):
            fr=folds[folds.target_proxy_group_id.eq(row['target_proxy_group_id'])&folds.readout.eq(row['readout'])] if len(folds) else folds
            vr=fields[fields.target_proxy_group_id.eq(row['target_proxy_group_id'])&fields.readout.eq(row['readout'])] if len(fields) else fields
            good=[]
            if len(fr):
                assert fr.fold.is_unique and set(fr.fold)=={0,1,2};fold_positions+=len(fr)
                for f in fr.to_dict('records'):
                    if f['status'] in ['ESTIMABLE_BLOCKED_PREDICTION','CONSTANT_HELDOUT_TARGET']:
                        x=vr[vr.fold.eq(f['fold'])];max_error=max(max_error,close(aggregate_direct(x),f));field_positions+=len(x)
                        assert set(x.field)<=included and f['train_FOVs']>=3 and f['train_cells_per_field']>=30
                        assert f['train_cells']==f['train_FOVs']*f['train_cells_per_field'] and f['train_cells_per_field']<=80
                    if f['status']=='ESTIMABLE_BLOCKED_PREDICTION':good.append(f['fold'])
            assert row['valid_folds']==len(good)
            if len(good)==3:
                assert set(vr.field)==included and vr.field.is_unique
                max_error=max(max_error,close(aggregate_direct(vr),row))
            else:
                assert row['status']=='INCOMPLETE_OR_INSUFFICIENT_SPATIAL_FOLDS'
                assert all(pd.isna(row[k]) for k in ['gain_R2','baseline_R2','multiview_R2']);unestimable+=1
        structural_checks.append(dict(image_ordinal=i,all_family_preserved=True,common_fields_and_train_selection_verified=True))
    checks=refit_fixed_case(control,plan,contract,prodout)
    receipt=dict(status='ACCEPTED_FULL126_AGGREGATION_AND_PRESELECTED_TWO_READOUT_DIRECT_REFIT',completed_UTC=datetime.now(timezone.utc).isoformat(),
                 audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(prodout/'private_terminal.json'),full_family_positions=126,
                 aggregate_fold_positions_rechecked=fold_positions,aggregate_test_field_positions_rechecked=field_positions,
                 maximum_summary_error=max_error,unestimable_preserved_as_missing=True,all_image_checks=structural_checks,direct_refit_checks=checks,
                 input_transformation_reused_from_hash_bound_capture_accepted_production=True,independent_raw_RNA_archive_replay=False,
                 all_forest_fits_independently_retrained=False,scientific_acceptance_pending=True,no_arrays_or_models_persisted=True)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','full_family_positions','aggregate_fold_positions_rechecked','aggregate_test_field_positions_rechecked','maximum_summary_error']}))


if __name__=='__main__':main()
