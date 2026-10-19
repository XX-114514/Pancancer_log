"""Outcome-free nested case/whole-study validation; no target spatial effects."""
import argparse
import hashlib
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from prepare_CM02_panel_bridge_reference_scope_v1 import budget

CLASSES=['B_Plasma','Myeloid','NK','Nonimmune','T']


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def weights(groups):
    _,codes,n=np.unique(groups,return_inverse=True,return_counts=True)
    w=1/n[codes]
    return w/w.mean()


def wrho(y,p,w):
    if len(y)<3 or np.ptp(y)==0 or np.ptp(p)==0:return np.nan
    x=rankdata(y);z=rankdata(p);x-=np.average(x,weights=w);z-=np.average(z,weights=w)
    d=np.sqrt(np.sum(w*x*x)*np.sum(w*z*z))
    return float(np.sum(w*x*z)/d) if d>0 else np.nan


def prepare(counts,cols):
    x=np.asarray(counts[:,cols],np.float64)
    denom=x.sum(1)
    good=denom>0
    x[good]=np.log1p(x[good]/denom[good,None]*10000)
    return x,good


def fit(x,y,groups,param,kind):
    keep=np.var(x,axis=0)>1e-12
    assert keep.any()
    scaler=StandardScaler().fit(x[:,keep])
    z=scaler.transform(x[:,keep]);w=weights(groups)
    if kind=='lineage':
        model=LogisticRegression(C=param,solver='lbfgs',max_iter=1000,tol=1e-4,
                                 class_weight='balanced',random_state=20261020)
    else:model=Ridge(alpha=param,solver='cholesky')
    with warnings.catch_warnings(record=True) as messages:
        warnings.simplefilter('always');model.fit(z,y,sample_weight=w)
    converged=not any(issubclass(m.category,ConvergenceWarning) for m in messages)
    return model,scaler,keep,converged


def predict(bundle,x,kind):
    model,scaler,keep,_=bundle;z=scaler.transform(x[:,keep])
    if kind=='lineage':
        p=np.zeros((len(x),len(CLASSES)))
        native=model.predict_proba(z)
        for j,label in enumerate(model.classes_):p[:,CLASSES.index(label)]=native[:,j]
        return p
    return model.predict(z)


def choose(x,y,groups,kind):
    grid=[.1,1.,10.] if kind=='lineage' else [100.,10.,1.,.1]
    splitter=GroupKFold(n_splits=3)
    rows=[]
    for param in grid:
        scores=[]
        for k,(tr,va) in enumerate(splitter.split(x,y,groups)):
            assert not set(groups[tr]) & set(groups[va])
            bundle=fit(x[tr],y[tr],groups[tr],param,kind)
            pred=predict(bundle,x[va],kind);w=weights(groups[va])
            if kind=='lineage':score=f1_score(y[va],np.array(CLASSES)[np.argmax(pred,axis=1)],average='macro',sample_weight=w,zero_division=0)
            else:score=float(np.nanmean([wrho(y[va,j],pred[:,j],w) for j in range(2)]))
            if not bundle[3]:score=np.nan
            rows.append(dict(parameter=param,inner_fold=k,score=score,train_cases=len(set(groups[tr])),validation_cases=len(set(groups[va])),converged=bundle[3]))
            scores.append(score)
        value=float(np.mean(scores)) if np.isfinite(scores).all() else -np.inf
        rows.append(dict(parameter=param,inner_fold='mean',score=value,converged=bool(np.isfinite(value))))
    means=[r for r in rows if r['inner_fold']=='mean']
    best=max(means,key=lambda r:r['score'])
    assert np.isfinite(best['score']),'No converged inner validation model'
    return best['parameter'],rows


def save_bundle(path,bundle,cols,tr,te):
    m,s,k,converged=bundle
    np.savez_compressed(path,coef=m.coef_,intercept=m.intercept_,mean=s.mean_,scale=s.scale_,
                        selected_feature_positions=np.asarray(cols)[k],denominator_positions=np.asarray(cols),training_rows=tr,test_rows=te,
                        classes=np.asarray(getattr(m,'classes_',[]),str),converged=np.array(converged))


def main(a):
    budget(a.root,2*1024**3)
    original=a.root/'reference_native_inputs_v2'
    audit=a.root/'reference_native_audit_v1/independent_native_audit_technical_public.json'
    assert json.loads(audit.read_text())['producer_terminal_sha256']==sha(original/'native_input_technical_public.json')
    cpath=a.root/'reference_scope_v1/private_frozen_reference_contract.json';c=json.loads(cpath.read_text())
    meta=pd.read_csv(original/'private_selected_reference_metadata.tsv.gz',sep='\t',keep_default_na=False)
    counts=np.load(original/'private_panel_counts.npy',mmap_mode='r')
    native=np.load(original/'private_full_native_target_arrays.npz')
    depth=json.loads(a.depth.read_text());target_depth=int(round(depth['panel_RNA_median']))
    assert target_depth>0
    out=a.root/'model_validation_v1';assert not out.exists();out.mkdir(mode=0o700);(out/'models').mkdir(mode=0o700)
    contract=dict(created_utc=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),
                  scope_sha256=sha(cpath),native_audit_sha256=sha(audit),depth_contract_sha256=sha(a.depth),
                  preprocessing='Each model uses its own permitted feature-only biological CP10K/log1p denominator; no excluded LR/core in denominator. Training-fold-only variance/scaler.',
                  target_depth=target_depth,stress='native-trained model evaluated at fixed multinomial matched target panel depth; no stress-driven hyperparameter tuning',
                  metrics='Patient-balanced confusion/F1/precision; Spearman ranks weighted to give every held-out case equal total weight; plain rho also retained.',
                  no_outcomes_geometry_or_spatial_effect_loaded=True,
                  all_study_folds=c['sources'],frozen_engineering_thresholds=c['eligibility'],
                  gate_requires_native_and_matched_depth_stress=True,
                  evaluation_class_cell_case_minima=[100,5],score_label='new_RNA_measurement_bridge_not_original_CM')
    (out/'plan_before_any_bridge_fit.json').write_text(json.dumps(contract,indent=2)+'\n')
    all_genes=c['shared_panel_genes'];lc=[all_genes.index(g) for g in c['lineage_feature_genes']];pc=[all_genes.index(g) for g in c['programme_feature_genes']]
    x_lineage,good_lineage=prepare(counts,lc);x_programme,good_programme=prepare(counts,pc)
    rng=np.random.default_rng(20261020);thin=np.zeros(counts.shape,np.int32)
    for k in range(len(meta)):
        row=np.asarray(counts[k],np.int64);total=int(row.sum())
        if total>target_depth:thin[k]=rng.multinomial(target_depth,row/total)
        else:thin[k]=row
    # Saved technical realization allows an independent numerical replay.
    np.save(out/'private_fixed_capture_stress_counts.npy',thin)
    sx_lineage,stress_good_lineage=prepare(thin,lc);sx_programme,stress_good_programme=prepare(thin,pc)
    rows=[];classrows=[];innerrows=[];predictions=[]
    groups=meta.patient_block.astype(str).to_numpy();sources=meta.dataset.astype(str).to_numpy()
    lineage_mask=meta.lineage_evaluation_row.astype(str).str.lower().isin(['true','1']).to_numpy() & good_lineage
    targets=native['original_full_signature_UCell1500']
    programme_mask=meta.lineage.eq('T').to_numpy() & np.isfinite(targets).all(1) & good_programme
    for holdout in c['sources']:
        for kind,valid,x,sx,cols,stressgood in [
            ('lineage',lineage_mask,x_lineage,sx_lineage,lc,stress_good_lineage),
            ('programme',programme_mask,x_programme,sx_programme,pc,stress_good_programme)]:
            tr=np.flatnonzero(valid & (sources!=holdout));te=np.flatnonzero(valid & (sources==holdout))
            assert len(tr)>0 and len(te)>0 and not (set(groups[tr]) & set(groups[te]))
            y=meta.lineage.astype(str).to_numpy() if kind=='lineage' else targets
            print(json.dumps(dict(phase='nested_fit_started',holdout=holdout,kind=kind,train=len(tr),test=len(te))),flush=True)
            param,inner=choose(x[tr],y[tr],groups[tr],kind)
            for r in inner:r.update(holdout=holdout,kind=kind)
            innerrows.extend(inner)
            bundle=fit(x[tr],y[tr],groups[tr],param,kind)
            assert bundle[3], 'Outer model did not converge'
            save_bundle(out/'models'/(holdout+'_'+kind+'.npz'),bundle,cols,tr,te)
            for mode,xx,ok in [('native',x,np.ones(len(meta),bool)),('matched_depth_stress',sx,stressgood)]:
                p=predict(bundle,xx[te],kind);w=weights(groups[te]);ncase=len(set(groups[te]))
                if kind=='lineage':
                    truth=y[te];label=np.array(CLASSES)[p.argmax(1)]
                    pred=pd.DataFrame(dict(reference_row=te,holdout=holdout,kind=kind,mode=mode,truth=truth,predicted=label,
                                           pressure_positive_input=ok[te],max_probability=p.max(1)))
                    for j,g in enumerate(CLASSES):pred['prob_'+g]=p[:,j]
                    predictions.append(pred)
                    matrix=confusion_matrix(truth,label,labels=CLASSES,sample_weight=w)
                    np.save(out/(holdout+'_'+mode+'_private_confusion.npy'),matrix)
                    macro=f1_score(truth,label,labels=CLASSES,average='macro',sample_weight=w,zero_division=0)
                    precision=precision_score(truth,label,labels=CLASSES,average=None,sample_weight=w,zero_division=0)
                    recall=recall_score(truth,label,labels=CLASSES,average=None,sample_weight=w,zero_division=0)
                    rows.append(dict(holdout=holdout,kind=kind,mode=mode,parameter=param,test_cells=len(te),test_cases=ncase,
                                     macroF1=macro,converged=True,zero_input_cells=int((~ok[te]).sum())))
                    for j,g in enumerate(CLASSES):
                        at=truth==g;nc=int(at.sum());ng=len(set(groups[te][at]));supported=nc>=100 and ng>=5
                        classrows.append(dict(holdout=holdout,mode=mode,lineage=g,cells=nc,cases=ng,
                                              precision=float(precision[j]),recall=float(recall[j]),coverage_eligible=supported))
                else:
                    for j,name in enumerate(['regulatory','cytotoxic']):
                        rho=wrho(y[te,j],p[:,j],w)
                        rows.append(dict(holdout=holdout,kind=kind,mode=mode,target=name,parameter=param,test_cells=len(te),test_cases=ncase,
                                         rho=rho,plain_cell_spearman=float(spearmanr(y[te,j],p[:,j]).statistic),
                                         RMSE=float(np.sqrt(np.average((y[te,j]-p[:,j])**2,weights=w))),
                                         coverage_eligible=len(te)>=100 and ncase>=5,converged=True,
                                         zero_input_cells=int((~ok[te]).sum())))
                        predictions.append(pd.DataFrame(dict(reference_row=te,holdout=holdout,kind=kind,mode=mode,target=name,
                                                              truth=y[te,j],predicted=p[:,j],pressure_positive_input=ok[te])))
            print(json.dumps(dict(phase='whole_study_fold_finished',holdout=holdout,kind=kind)),flush=True)
    metrics=pd.DataFrame(rows);classes=pd.DataFrame(classrows)
    metrics.to_csv(out/'private_all_model_metrics.tsv',sep='\t',index=False)
    classes.to_csv(out/'private_all_class_precision_recall.tsv',sep='\t',index=False)
    pd.DataFrame(innerrows).to_csv(out/'private_all_inner_case_validation.tsv',sep='\t',index=False)
    pd.concat(predictions,ignore_index=True).to_csv(out/'private_all_heldout_predictions.tsv.gz',sep='\t',index=False,compression='gzip')
    gates=[]
    for mode in ['native','matched_depth_stress']:
        line=metrics[(metrics.kind=='lineage') & (metrics['mode']==mode)]
        global_ok=len(line)>=3 and line.macroF1.ge(.8).all()
        for label in CLASSES:
            subset=classes[(classes['mode']==mode) & classes.lineage.eq(label) & classes.coverage_eligible]
            gates.append(dict(mode=mode,target=label,kind='lineage',passed=bool(global_ok and len(subset)>=3 and subset.precision.ge(.9).all()),
                              eligible_sources=len(subset),reason='Complete macroF1>=.80 and >=3 held studies with >=100cells/5cases and classprecision>=.90'))
        for target in ['regulatory','cytotoxic']:
            subset=metrics[(metrics.kind=='programme') & (metrics['mode']==mode) & metrics.target.eq(target) & metrics.coverage_eligible.fillna(False)]
            gates.append(dict(mode=mode,target=target,kind='programme',passed=bool(len(subset)>=3 and subset.rho.notna().all() and subset.rho.median()>=.5 and subset.rho.ge(.3).all()),
                              eligible_sources=len(subset),reason='>=3 held studies with >=100cells/5cases; medianrho>=.50 and eachrho>=.30'))
    pd.DataFrame(gates).to_csv(out/'private_all_frozen_engineering_gates.tsv',sep='\t',index=False)
    spatial_ok=all(g['passed'] for g in gates if g['target'] in ['T','regulatory','cytotoxic'])
    (out/'private_spatial_admission_decision.json').write_text(json.dumps(dict(
        full_dual_arm_bridge_passed=spatial_ok,gate_rows=gates,
        no_relaxed_thresholds_or_target_effect_selection=True,
        spatial_effects_not_run_by_this_script=True),indent=2)+'\n')
    result=dict(status='ALL_FROZEN_NESTED_STUDY_CASE_VALIDATION_MODELS_COMPLETE_INDEPENDENT_AUDIT_PENDING',
                finished_utc=datetime.now(timezone.utc).isoformat(),outer_studies=len(c['sources']),
                outer_models=len(list((out/'models').glob('*.npz'))),
                producer_sha256=sha(Path(__file__)),fit_contract_sha256=sha(out/'plan_before_any_bridge_fit.json'),
                all_model_files_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},
                effects_predictions_coefficients_and_gate_values_private=True,no_spatial_P=True)
    (out/'model_validation_technical_public.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='all_model_files_sha256'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for field in ['root','depth']:p.add_argument('--'+field,type=Path,required=True)
    main(p.parse_args())
