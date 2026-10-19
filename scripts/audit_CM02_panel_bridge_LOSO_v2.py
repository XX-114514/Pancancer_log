"""Replay all outer predictions from coefficients; independent metrics/gates."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def rank(x):
    _,inverse,n=np.unique(x,return_inverse=True,return_counts=True)
    return (np.cumsum(n)-(n-1)/2)[inverse].astype(float)


def case_weights(groups):
    g=pd.Series(groups);n=g.map(g.value_counts()).to_numpy(float);w=1/n
    return w/w.mean()


def rho(y,p,w):
    x=rank(y);z=rank(p);x-=np.sum(w*x)/w.sum();z-=np.sum(w*z)/w.sum()
    den=np.sqrt(np.dot(w*x,x)*np.dot(w*z,z))
    return float(np.dot(w*x,z)/den) if den>0 else np.nan


def equal(a,b,tol=2e-10):
    if pd.isna(a) and pd.isna(b):return
    assert np.isfinite(a) and np.isfinite(b) and abs(a-b)<tol,(a,b)


def replay(model,counts,kind,classes):
    total=np.asarray(counts[:,model['denominator_positions']],float).sum(1)
    xx=np.asarray(counts[:,model['selected_feature_positions']],float)
    nz=total>0;xx[nz]=np.log1p(xx[nz]/total[nz,None]*10000)
    xx=(xx-model['mean'])/model['scale']
    values=xx@model['coef'].T+model['intercept']
    if kind=='programme':return values
    sourceclasses=model['classes'].astype(str).tolist()
    if len(sourceclasses)==2 and values.shape[1]==1:
        positive=1/(1+np.exp(-values[:,0]));prob=np.column_stack([1-positive,positive])
    else:
        values-=values.max(1)[:,None];prob=np.exp(values);prob/=prob.sum(1)[:,None]
    full=np.zeros((len(xx),len(classes)))
    for k,c in enumerate(sourceclasses):full[:,classes.index(c)]=prob[:,k]
    return full


def main(a):
    out=a.root/'model_validation_v1';proof=json.loads((out/'model_validation_technical_public.json').read_text())
    for name,digest in proof['all_model_files_sha256'].items():assert sha(out/name)==digest
    meta=pd.read_csv(a.root/'reference_native_inputs_v2/private_selected_reference_metadata.tsv.gz',sep='\t',keep_default_na=False)
    counts=np.load(a.root/'reference_native_inputs_v2/private_panel_counts.npy',mmap_mode='r')
    stress=np.load(out/'private_fixed_capture_stress_counts.npy',mmap_mode='r')
    archive=np.load(a.root/'reference_native_inputs_v2/private_full_native_target_arrays.npz')
    target={k:archive[k] for k in archive.files}
    pred=pd.read_csv(out/'private_all_heldout_predictions.tsv.gz',sep='\t',float_precision='round_trip')
    metrics=pd.read_csv(out/'private_all_model_metrics.tsv',sep='\t')
    classmetrics=pd.read_csv(out/'private_all_class_precision_recall.tsv',sep='\t')
    inner=pd.read_csv(out/'private_all_inner_case_validation.tsv',sep='\t')
    groups=meta.patient_block.astype(str).to_numpy();classes=['B_Plasma','Myeloid','NK','Nonimmune','T']
    max_prediction_error=0.;all_prediction_values=0;models=0
    contract=json.loads((a.root/'reference_scope_v1/private_frozen_reference_contract.json').read_text())
    diagnostics=pd.read_csv(out/'private_missing_spatial_marker_recovery_diagnostics.tsv',sep='\t')
    for modelpath in sorted((out/'models').glob('*.npz')):
        holdout,kind=modelpath.stem.rsplit('_',1);m=np.load(modelpath,allow_pickle=False)
        tr,te=m['training_rows'],m['test_rows'];assert len(set(tr)&set(te))==0
        assert not set(groups[tr])&set(groups[te])
        assert meta.dataset.iloc[te].eq(holdout).all() and not meta.dataset.iloc[tr].eq(holdout).any()
        forbidden=set(contract['LR_excluded_genes']+(contract['core14'] if kind=='programme' else []))
        assert not set(np.asarray(contract['shared_panel_genes'])[m['denominator_positions']])&forbidden
        # Re-evaluate the complete grid's means and the predeclared tie order.
        tab=inner[(inner.holdout==holdout)&(inner.kind==kind)]
        grid=[.1,1.,10.] if kind=='lineage' else [100.,10.,1.,.1]
        scored=[]
        for value in grid:
            chunk=tab[tab.parameter.eq(value)];parts=chunk[chunk.inner_fold.astype(str).ne('mean')]
            assert len(parts)==3
            mean=float(parts.score.mean()) if parts.score.notna().all() else -np.inf
            recorded=float(chunk[chunk.inner_fold.astype(str).eq('mean')].score.iloc[0]);equal(mean,recorded) if np.isfinite(mean) else None
            scored.append(mean)
        expected_parameter=grid[int(np.argmax(scored))]
        assert metrics[(metrics.holdout==holdout)&(metrics.kind==kind)].parameter.eq(expected_parameter).all()
        for mode,source in [('native',counts),('matched_depth_stress',stress)]:
            expected=replay(m,source[te],kind,classes);w=case_weights(groups[te])
            subset=pred[(pred.holdout==holdout)&(pred.kind==kind)&(pred['mode']==mode)]
            if kind=='lineage':
                d=subset.set_index('reference_row').loc[te]
                stored=d[['prob_'+x for x in classes]].to_numpy()
                max_prediction_error=max(max_prediction_error,float(np.max(abs(expected-stored))));all_prediction_values+=expected.size
                truth=meta.lineage.iloc[te].to_numpy();labels=np.asarray(classes)[expected.argmax(1)]
                assert np.array_equal(labels,d.predicted.to_numpy()) and np.array_equal(truth,d.truth.to_numpy())
                cm=np.zeros((len(classes),len(classes)))
                for i,c in enumerate(classes):
                    for j,k in enumerate(classes):cm[i,j]=w[(truth==c)&(labels==k)].sum()
                assert np.max(abs(cm-np.load(out/(holdout+'_'+mode+'_private_confusion.npy'))))<1e-8
                f=[]
                for j,c in enumerate(classes):
                    precision=cm[j,j]/cm[:,j].sum() if cm[:,j].sum() else 0.
                    recall=cm[j,j]/cm[j].sum() if cm[j].sum() else 0.
                    f.append(2*cm[j,j]/(cm[j].sum()+cm[:,j].sum()) if cm[j].sum()+cm[:,j].sum() else 0.)
                    row=classmetrics[(classmetrics.holdout==holdout)&(classmetrics['mode']==mode)&(classmetrics.lineage==c)].iloc[0]
                    equal(precision,row.precision);equal(recall,row.recall)
                    assert int(row.cells)==int((truth==c).sum()) and int(row.cases)==len(set(groups[te][truth==c]))
                row=metrics[(metrics.holdout==holdout)&(metrics.kind==kind)&(metrics['mode']==mode)].iloc[0]
                equal(np.mean(f),row.macroF1)
            else:
                for j,name in enumerate(['regulatory','cytotoxic']):
                    d=subset[subset.target==name].set_index('reference_row').loc[te]
                    max_prediction_error=max(max_prediction_error,float(np.max(abs(expected[:,j]-d.predicted.to_numpy(float)))));all_prediction_values+=len(te)
                    # Preserve exact serialized ties for independent metric ranks;
                    # coefficient reconstruction above has already checked values.
                    actual_prediction=d.predicted.to_numpy(float)
                    truth=target['original_full_signature_UCell1500'][te,j];np.testing.assert_allclose(d.truth.to_numpy(float),truth,atol=1e-14,rtol=0)
                    row=metrics[(metrics.holdout==holdout)&(metrics.kind==kind)&(metrics['mode']==mode)&(metrics.target==name)].iloc[0]
                    equal(rho(truth,actual_prediction,w),row.rho)
                    equal(rho(truth,actual_prediction,np.ones(len(te))),row.plain_cell_spearman)
                    equal(np.sqrt(np.dot(w,(truth-actual_prediction)**2)/w.sum()),row.RMSE)
                    for gene in (['TNFRSF18'] if name=='regulatory' else ['CTSW','GZMB','NKG7']):
                        raw=target['core14_native_counts'][te,contract['core14'].index(gene)]
                        value=np.log1p(raw/target['full_native_total'][te]*10000)
                        dr=diagnostics[(diagnostics.holdout==holdout)&(diagnostics['mode']==mode)&(diagnostics.target==name)&(diagnostics.missing_spatial_gene==gene)].iloc[0]
                        equal(rho(value,actual_prediction,w),dr.programme_prediction_vs_native_gene_rho)
        models+=1
    assert max_prediction_error<2e-10
    gates=pd.read_csv(out/'private_all_frozen_engineering_gates.tsv',sep='\t')
    for g in gates.itertuples():
        if g.kind=='lineage':
            allmodels=metrics[(metrics.kind=='lineage')&(metrics['mode']==g.mode)]
            sub=classmetrics[(classmetrics['mode']==g.mode)&classmetrics.lineage.eq(g.target)&classmetrics.coverage_eligible]
            passed=len(allmodels)>=3 and allmodels.macroF1.ge(.8).all() and len(sub)>=3 and sub.precision.ge(.9).all()
        else:
            sub=metrics[(metrics.kind=='programme')&(metrics['mode']==g.mode)&metrics.target.eq(g.target)&metrics.coverage_eligible.fillna(False)]
            passed=len(sub)>=3 and sub.rho.notna().all() and sub.rho.median()>=.5 and sub.rho.ge(.3).all()
        assert bool(g.passed)==bool(passed)
    dual=bool(gates[gates.target.isin(['T','regulatory','cytotoxic'])].passed.all())
    assert dual==json.loads((out/'private_spatial_admission_decision.json').read_text())['full_dual_arm_bridge_passed']
    destination=a.root/'model_independent_audit_v2';assert not destination.exists();destination.mkdir(mode=0o700)
    result=dict(status='ALL_OUTER_MODEL_PREDICTIONS_PATIENT_SOURCE_SPLITS_METRICS_AND_ENGINEERING_GATES_INDEPENDENTLY_ACCEPTED',
                finished_utc=datetime.now(timezone.utc).isoformat(),outer_models=models,
                replayed_prediction_values=all_prediction_values,maximum_prediction_difference=max_prediction_error,
                inner_parameter_selection_replayed_from_all_frozen_scores=True,
                no_claim_all_inner_optimizers_independently_refit=True,
                producer_terminal_sha256=sha(out/'model_validation_technical_public.json'),auditor_sha256=sha(Path(__file__)),
                scientific_metrics_gates_predictions_not_in_public_receipt=True,no_spatial_inference_run=True)
    (destination/'independent_model_audit_technical_public.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args())
