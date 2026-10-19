"""Source-heldout recovery of operational RNA classes; all learned arrays RAM only."""
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
import hashlib
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix, vstack
from scipy.special import softmax
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from sklearn.utils.extmath import randomized_svd


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def digest_rows(rows):
    return hashlib.sha256(np.asarray(rows, dtype='<i8').tobytes()).hexdigest()


def seed(text):
    return int(hashlib.sha256(text.encode()).hexdigest()[:8], 16)


def select_rows(frame, task, cfg):
    selected = []
    for (source, label), g in frame.groupby(['source', task], sort=True):
        groups = []
        for library, h in g.groupby('library', sort=True):
            prefix = str(cfg['seed']) + '|' + task + '|' + source + '|' + label + '|'
            rows = sorted(h.row.to_numpy(), key=lambda row: hashlib.sha256((prefix + str(int(row))).encode()).digest())
            groups.append((seed(prefix + library), rows[:cfg['library_class_cap']]))
        groups.sort(key=lambda x: x[0])
        chosen = []
        for position in range(cfg['library_class_cap']):
            for _, rows in groups:
                if position < len(rows): chosen.append(int(rows[position]))
                if len(chosen) == cfg['source_class_cap']: break
            if len(chosen) == cfg['source_class_cap']: break
        selected.extend(chosen)
    out = np.sort(np.array(selected, np.int64))
    assert len(out) == len(set(out))
    return out


def balanced_weights(y, source, library):
    f = pd.DataFrame(dict(label=y, source=source, library=library))
    size = f.groupby(['label', 'source', 'library']).size()
    units = f.groupby(['label', 'source']).library.nunique()
    studies = f.groupby('label').source.nunique()
    w = np.array([1. / (studies[label] * units[label, src] * size[label, src, lib]) for label, src, lib in zip(y, source, library)])
    w *= len(w) / w.sum()
    assert np.isfinite(w).all() and (w > 0).all()
    return w


def scale(train, query, w):
    mu = np.average(train, axis=0, weights=w)
    sd = np.sqrt(np.average((train - mu)**2, axis=0, weights=w))
    if not np.isfinite(sd).all() or np.any(sd < 1e-10): return None
    return (train - mu) / sd, (query - mu) / sd


def fit_model(train, test, y, w, classes):
    if not np.isfinite(train).all() or not np.isfinite(test).all(): return None, 'NONFINITE_DESIGN'
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        model = LogisticRegression(C=1., solver='lbfgs', max_iter=2000, tol=1e-8)
        model.fit(train, y, sample_weight=w)
    if any(type(x.message).__name__ == 'ConvergenceWarning' for x in ws): return None, 'OPTIMIZER_NOT_CONVERGED'
    assert model.classes_.tolist() == classes
    logits = train @ model.coef_.T + model.intercept_
    p = softmax(logits, axis=1)
    np.testing.assert_allclose(p, model.predict_proba(train), atol=1e-11, rtol=1e-11)
    target = np.eye(len(classes))[pd.Index(classes).get_indexer(y)]
    residual = w[:, None] * (p - target)
    gradient = residual.T @ train + model.coef_
    gradient_intercept = residual.sum(axis=0)
    error = max(float(np.abs(gradient).max()), float(np.abs(gradient_intercept).max())) / w.sum()
    assert error < 1e-5
    direct_loss = float(np.average(-np.log(np.maximum(p, 1e-300))[np.arange(len(y)), target.argmax(axis=1)], weights=w))
    np.testing.assert_allclose(direct_loss, log_loss(y, p, sample_weight=w, labels=classes), atol=1e-10, rtol=1e-10)
    probabilities = softmax(test @ model.coef_.T + model.intercept_, axis=1)
    np.testing.assert_allclose(probabilities, model.predict_proba(test), atol=1e-11, rtol=1e-11)
    np.testing.assert_allclose(probabilities.sum(axis=1), 1, atol=1e-12)
    return dict(probabilities=probabilities, training_scaled_KKT_error=error, n_iter=int(model.n_iter_.max())), ''


def metrics(y, probabilities, w, classes):
    truth = pd.Index(classes).get_indexer(y)
    predicted = probabilities.argmax(axis=1)
    loss = -np.log(np.clip(probabilities[np.arange(len(y)), truth], 1e-15, 1))
    np.testing.assert_allclose(np.average(loss, weights=w), log_loss(y, probabilities, sample_weight=w, labels=classes), atol=1e-10)
    k = len(classes)
    confusion = np.bincount(truth * k + predicted, weights=w, minlength=k * k).reshape(k, k)
    raw = np.bincount(truth * k + predicted, minlength=k * k).reshape(k, k)
    loop = np.zeros((k, k))
    for a, b, v in zip(truth, predicted, w): loop[a, b] += v
    np.testing.assert_allclose(confusion, loop, rtol=1e-12, atol=1e-10)
    recall = np.diag(confusion) / confusion.sum(axis=1)
    precision = np.divide(np.diag(confusion), confusion.sum(axis=0), out=np.zeros(k), where=confusion.sum(axis=0) > 0)
    f1 = np.divide(2 * precision * recall, precision + recall, out=np.zeros(k), where=(precision + recall) > 0)
    return dict(balanced_material_class_logloss=float(np.average(loss, weights=w)), sampled_cell_logloss=float(loss.mean()), balanced_accuracy=float(recall.mean()), macroF1=float(f1.mean())), loss, raw


def bootstrap_delta(y, library, delta, classes, random_seed):
    libs = sorted(set(library)); table = np.full((len(libs), len(classes)), np.nan)
    for i, lib in enumerate(libs):
        for j, label in enumerate(classes):
            take = (library == lib) & (y == label)
            if take.any(): table[i, j] = delta[take].mean()
    rng = np.random.default_rng(random_seed); values = []
    for _ in range(2000):
        draw = rng.integers(0, len(libs), len(libs))
        count = np.bincount(draw, minlength=len(libs))
        denominator = count @ np.isfinite(table)
        if (denominator == 0).any(): values.append(np.nan); continue
        value = float(np.mean((count @ np.nan_to_num(table)) / denominator))
        alternate = float(np.mean([np.mean([table[i, j] for i in draw if np.isfinite(table[i, j])]) for j in range(len(classes))]))
        assert abs(value - alternate) < 1e-12
        values.append(value)
    a = np.asarray(values); valid = int(np.isfinite(a).sum())
    result = dict(declared_draws=2000, valid_draws=valid, failed_missing_class_draws=2000-valid, conditional_fixed_prediction_library_resampling=True)
    if valid >= 1800:
        lo, hi = np.quantile(a[np.isfinite(a)], [.025, .975]); result.update(CI_low=float(lo), CI_high=float(hi))
    else: result['interval_status'] = 'INSUFFICIENT_COMPLETE_CLASS_RESAMPLES'
    return result


def RNA_coordinates(train, test, w, cfg, random_seed):
    first = np.asarray(train.T @ w).ravel() / w.sum()
    second = np.asarray(train.power(2).T @ w).ravel() / w.sum()
    frequency = np.asarray((train > 0).T @ w).ravel() / w.sum()
    variance = second - first**2
    eligible = np.flatnonzero((frequency >= .01) & (variance > 1e-10))
    if len(eligible) < cfg['HVG']: return None, 'INSUFFICIENT_TRAIN_ONLY_GENE_FEATURES'
    order = np.argsort(-variance[eligible], kind='stable')[:cfg['HVG']]
    columns = eligible[order]
    a = train[:, columns].toarray(); b = test[:, columns].toarray()
    mu = np.average(a, axis=0, weights=w)
    np.testing.assert_allclose(mu, first[columns], atol=1e-10, rtol=1e-10)
    centered = (a - mu) * np.sqrt(w[:, None])
    _, singular, loading = randomized_svd(centered, n_components=cfg['PCA_components'], n_iter=4, random_state=random_seed)
    if singular.min() <= 1e-10: return None, 'TRAIN_COMPONENT_RANK_FAILURE'
    pc_train = (a - mu) @ loading.T; pc_test = (b - mu) @ loading.T
    return dict(train=pc_train, test=pc_test, feature_indices_SHA256=digest_rows(columns), train_only_feature_count=len(columns), input_genes=len(first)), ''


def fixture():
    rng = np.random.default_rng(2026102006)
    y = np.repeat(['a', 'b', 'c'], 180)
    x = rng.uniform(0, .3, (len(y), 90))
    for k in range(3): x[k*180:(k+1)*180, k*20:(k+1)*20] += 2
    train = np.tile(np.r_[np.ones(120, bool), np.zeros(60, bool)], 3)
    w = np.ones(int(train.sum()))
    pc, reason = RNA_coordinates(csr_matrix(x[train]), csr_matrix(x[~train]), w, {'HVG':60, 'PCA_components':5}, 123)
    assert not reason
    scaled = scale(pc['train'], pc['test'], w); assert scaled is not None
    model, reason = fit_model(*scaled, y[train], w, ['a','b','c']); assert not reason
    m, _, _ = metrics(y[~train], model['probabilities'], np.ones((~train).sum()), ['a','b','c'])
    assert m['balanced_accuracy'] > .95
    assert scale(np.ones((5,2)), np.ones((3,2)), np.ones(5)) is None
    no, reason = RNA_coordinates(csr_matrix(np.zeros((10,90))), csr_matrix(np.zeros((4,90))), np.ones(10), {'HVG':60,'PCA_components':5}, 1)
    assert no is None and reason == 'INSUFFICIENT_TRAIN_ONLY_GENE_FEATURES'
    return dict(status='SYNTHETIC_SEPARABLE_RECOVERY_CONSTANT_QC_AND_UNESTIMABLE_RNA_ACCEPTED', real_inputs_used=False)


def main():
    cfg = json.load(sys.stdin)
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    assert sha(Path(os.environ['PANCANCER_REFERENCE_QUALIFICATION'])) == cfg['qualification_SHA256']
    from qualify_CM_lineage_state_reference_inputs_v1 import label_T
    old = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    assert sha(old / 'private_fixed34_gene_T_counts_and_depth.npz') == cfg['T_cache_SHA256']
    cache = dict(np.load(old / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False))
    fixed = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed) == cfg['fixed_panel_contract_SHA256']; panels = json.loads(fixed.read_text())['panels']
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    assert sha(atlas / 'manifests/immune_authority_release_public.json') == cfg['release_SHA256']
    cm_path = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    assert sha(cm_path) == cfg['CM_material_SHA256']
    cm = pd.read_csv(cm_path, sep='\t', keep_default_na=False, usecols=['sampleID'])
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    md = pd.DataFrame({k:np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset','library_id','major']})
    md['row'] = np.arange(len(md)); md['source'] = md.dataset; md['library'] = md.dataset + '::' + md.library_id
    material = md.library.isin(cm.sampleID).to_numpy(); assert sorted(set(md.loc[material,'source'])) == cfg['source_scope']
    tr = cache['global_original_T_row_index']; np.testing.assert_array_equal(tr, np.flatnonzero(md.major.eq('T')))
    md['coarse_lineage'] = 'Other_assigned'
    md.loc[md.major.isin(['T','NK']), 'coarse_lineage'] = md.loc[md.major.isin(['T','NK']), 'major']
    md.loc[md.major.isin(cfg['boundary_majors']), 'coarse_lineage'] = 'Boundary_unresolved'
    md['T_capture_state'] = ''
    md.loc[tr,'T_capture_state'] = label_T(cache['raw_gene_counts'], cache['genes'].tolist(), panels, cfg['anchor_gene'], cfg['effector_genes'], 1)
    selections = {}
    for task in cfg['tasks']:
        keep = material & (md.major.eq('T').to_numpy() if task == 'T_capture_state' else True)
        selections[task] = select_rows(md.loc[keep], task, cfg)
    rows = np.unique(np.concatenate(list(selections.values())))
    assert len(rows) < 40000
    source_order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured = np.asarray(read_elem(z['varm/measured_in_source']),bool)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    status = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    maps = json.loads((old / 'private_source_gene_measurement_contract.json').read_text())
    excluded = {int(v) for item in maps for v in item['physical_columns'] if v>=0}
    good = measured[:,[source_order.index(s) for s in cfg['source_scope']]].all(axis=1)
    good &= ~np.char.startswith(canonical,'AMBIGUOUS::') & (status!='ambiguous_source_feature')
    good &= ~pd.Index(canonical).duplicated(keep=False); good[list(excluded)] = False
    common = np.flatnonzero(good); assert len(common) == cfg['expected_source_common_disjoint_genes'] and not set(common)&excluded
    out.mkdir(mode=0o700)
    selection_proof = {task:dict(cells=len(rr), native_rows_SHA256=digest_rows(rr)) for task,rr in selections.items()}
    (out/'private_plan_before_RNA_and_model_fitting.json').write_text(json.dumps(dict(plan=cfg, selection=selection_proof, union_native_rows_SHA256=digest_rows(rows), common_physical_gene_indices_SHA256=digest_rows(common), source_SHA256=sha(__file__)),indent=2)+'\n')
    counts=z['layers/counts']; nc,nf=map(int,counts.attrs['shape']); pieces=[]; total=[]; detected=[]; seen=0
    for start in range(0,nc,cfg['chunk_rows']):
        stop=min(nc,start+cfg['chunk_rows']); lo=int(np.searchsorted(rows,start)); hi=int(np.searchsorted(rows,stop))
        if hi==lo: continue
        ptr=np.asarray(counts['indptr'][start:stop+1],np.int64); a,b=int(ptr[0]),int(ptr[-1])
        x=csr_matrix((np.asarray(counts['data'][a:b]),np.asarray(counts['indices'][a:b]),ptr-ptr[0]),shape=(stop-start,nf))[rows[lo:hi]-start]
        x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
        assert np.isfinite(x.data).all() and (x.data>0).all() and (x.data==np.floor(x.data)).all()
        local_source=md.source.to_numpy()[rows[lo:hi]]
        for src in set(local_source): assert measured[x[local_source==src].indices,source_order.index(src)].all()
        tt=md.major.to_numpy()[rows[lo:hi]]=='T'
        totals=np.asarray(x.sum(axis=1)).ravel().astype(np.int64); det=np.diff(x.indptr)
        if tt.any():
            pos=np.searchsorted(tr,rows[lo:hi][tt]); np.testing.assert_array_equal(tr[pos],rows[lo:hi][tt])
            np.testing.assert_array_equal(totals[tt],cache['full_native_total_counts'][pos]); np.testing.assert_array_equal(det[tt],cache['full_native_detected_genes'][pos])
            for src in set(local_source[tt]):
                mp=next(item for item in maps if item['dataset']==src); cc=np.array(mp['physical_columns']); assert (cc>=0).all()
                take=tt&(local_source==src); old_pos=np.searchsorted(tr,rows[lo:hi][take]); np.testing.assert_array_equal(x[take][:,cc].toarray(),cache['raw_gene_counts'][old_pos])
        pieces.append(x[:,common].astype(np.float64)); total.extend(totals.tolist()); detected.extend(det.tolist()); seen+=hi-lo
        (out/'private_progress.json').write_text(json.dumps(dict(phase='EXTRACTING_EXISTING_RNA',original_rows_scanned=stop,selected_rows_read=seen,selected_rows_total=len(rows)))+'\n')
    X=vstack(pieces,format='csr'); del pieces
    total=np.asarray(total);detected=np.asarray(detected); assert X.shape==(len(rows),len(common)) and (total>0).all()
    logX=X.multiply(10000/total[:,None]).tocsr(); logX.data=np.log1p(logX.data); del X
    QC=np.c_[np.log1p(total),np.log1p(detected)]
    results=[];class_scores=[];confusions=[];numeric=[]
    for task in cfg['tasks']:
        rr=selections[task]; ix=np.searchsorted(rows,rr); frame=md.iloc[rr].reset_index(drop=True)
        classes=sorted(cfg['coarse_classes'] if task=='coarse_lineage' else cfg['T_classes'])
        assert set(frame[task])==set(classes)
        for heldout in cfg['source_scope']:
            train=frame.source.ne(heldout).to_numpy(); test=~train
            train_ix=ix[train];test_ix=ix[test]; yt=frame.loc[train,task].to_numpy(); yv=frame.loc[test,task].to_numpy()
            record=dict(task=task,heldout_source=heldout,training_cells=int(train.sum()),heldout_cells=int(test.sum()),heldout_material_units=int(frame.loc[test,'library'].nunique()),status='UNESTIMABLE')
            support=[]
            for label in classes:
                a=frame.loc[train&frame[task].eq(label).to_numpy()]; b=frame.loc[test&frame[task].eq(label).to_numpy()]
                support.append(dict(candidate=label,train_cells=len(a),train_sources=int(a.source.nunique()),train_material_units=int(a.library.nunique()),heldout_cells=len(b),heldout_material_units=int(b.library.nunique())))
            record['class_support']=support
            if any(s['train_cells']<20 or s['train_sources']<2 or s['train_material_units']<3 for s in support):
                record['reason']='FULL_CLASS_TRAIN_SUPPORT_GATE';results.append(record);continue
            assert not set(frame.loc[train,'library'])&set(frame.loc[test,'library'])
            w=balanced_weights(yt,frame.loc[train,'source'].to_numpy(),frame.loc[train,'library'].to_numpy())
            v=balanced_weights(yv,frame.loc[test,'source'].to_numpy(),frame.loc[test,'library'].to_numpy())
            qs=scale(QC[train_ix],QC[test_ix],w)
            if qs is None: record['reason']='CONSTANT_REQUIRED_QC';results.append(record);continue
            coords,reason=RNA_coordinates(logX[train_ix],logX[test_ix],w,cfg,seed(str(cfg['seed'])+'|'+task+'|'+heldout))
            if coords is None: record['reason']=reason;results.append(record);continue
            pc=scale(coords['train'],coords['test'],w)
            if pc is None: record['reason']='CONSTANT_REQUIRED_TRAIN_COMPONENT';results.append(record);continue
            predictions={};metric={};losses={};failed=''
            for model,(a,b) in [('QC_only',qs),('RNA_PCs_plus_QC',(np.c_[pc[0],qs[0]],np.c_[pc[1],qs[1]]))]:
                fit,reason=fit_model(a,b,yt,w,classes)
                if fit is None: failed=reason;break
                metric[model],losses[model],raw=metrics(yv,fit['probabilities'],v,classes)
                predictions[model]=fit['probabilities']
                numeric.append(dict(task=task,heldout_source=heldout,model=model,training_scaled_KKT_error=fit['training_scaled_KKT_error'],n_iter=fit['n_iter'],feature_indices_SHA256=coords['feature_indices_SHA256'],heldout_source_excluded_from_training=True))
                for ai,label in enumerate(classes):
                    take=yv==label
                    class_scores.append(dict(task=task,heldout_source=heldout,model=model,candidate=label,heldout_cells=int(take.sum()),heldout_material_units=int(frame.loc[test].loc[take,'library'].nunique()),recall=float(raw[ai,ai]/raw[ai].sum()) if raw[ai].sum() else None,balanced_material_logloss=float(np.average(losses[model][take],weights=v[take])) if take.any() else None))
                    for bi,pred in enumerate(classes): confusions.append(dict(task=task,heldout_source=heldout,model=model,observed_candidate=label,predicted_candidate=pred,cells=int(raw[ai,bi])))
            if failed: record['reason']=failed;results.append(record);continue
            full_support=all(s['heldout_cells']>=20 and s['heldout_material_units']>=2 for s in support)
            if not full_support:
                for m in metric.values(): m['balanced_accuracy']=None;m['macroF1']=None
            delta=losses['QC_only']-losses['RNA_PCs_plus_QC']
            resampling=bootstrap_delta(yv,frame.loc[test,'library'].to_numpy(),delta,classes,seed('bootstrap|'+task+'|'+heldout))
            record.update(status='COMPLETE_FIXED_SOURCE_HELDOUT_DIAGNOSTIC',all_heldout_classes_support_full_macro=full_support,metrics=metric,paired_class_material_balanced_logloss_improvement=float(np.average(delta,weights=v)),conditional_library_bootstrap=resampling,training_feature_indices_SHA256=coords['feature_indices_SHA256'])
            results.append(record)
            (out/'private_progress.json').write_text(json.dumps(dict(phase='FITTING_SOURCE_HELDOUT_REFERENCE',completed_contexts=len(results),registered_contexts=8))+'\n')
    assert len(results)==8
    def save(name,value): (out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    save('private_all8_reference_contexts.json',results);save('private_all_source_class_prediction_descriptions.json',class_scores);save('private_all_source_candidate_confusion_counts.json',confusions);save('private_all_reference_numeric_proofs.json',numeric)
    files={p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(p['bytes'] for p in files.values())<8*1024**2
    save('private_terminal.json',dict(status='COMPLETE_ALL8_MARKER_DISJOINT_REFERENCE_DIAGNOSTICS_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),files=files,all_real_RNA_learned_features_weights_predictions_and_bootstrap_RAM_only=True,original_project_read_only=True,new_clinical_P_FDR_or_formal_identity=False))
    print(json.dumps(dict(status='COMPLETE_REFERENCE_REVIEW_PENDING',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__': main()
