"""Held-source RNA-mixture recovery; profiles/weights/mixtures stay in RAM."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.linalg import solve_triangular
from scipy.optimize import nnls
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def seed(s): return int(hashlib.sha256(s.encode()).hexdigest()[:8], 16)


def geometry(reference, weights):
    a = reference * weights[:, None]
    singular = np.linalg.svd(a, compute_uv=False)
    if singular[-1] <= singular[0] * 1e-10: return None, 'FULL_BASIS_RANK_DEFICIENT'
    condition = float(singular[0] / singular[-1])
    if condition > 1e4: return None, 'FULL_BASIS_NUMERICAL_CONDITION_GATE'
    gram = a.T @ a
    chol = np.linalg.cholesky(gram)
    np.testing.assert_allclose(chol @ chol.T, gram, rtol=1e-10, atol=1e-12)
    return dict(A=a, G=gram, L=chol, condition=condition), ''


def recover(reference, teacher, truth, weights, problem, direct_cases):
    cross = reference.T @ ((weights**2)[:, None] * teacher)
    teacher_gram = teacher.T @ ((weights**2)[:, None] * teacher)
    rhs = cross @ truth
    target_norm2 = np.sum(truth * (teacher_gram @ truth), axis=0)
    output = np.zeros((reference.shape[1], truth.shape[1]))
    residuals = []; errors = []; direct_errors = []; amplitude_sums = []
    for i in range(truth.shape[1]):
        c = rhs[:, i]
        b = solve_triangular(problem['L'], c, lower=True)
        x, _ = nnls(problem['L'].T, b, maxiter=1000)
        gradient = problem['G'] @ x - c
        scale = max(float(np.max(np.abs(c))), 1e-15)
        kkt = max(float(np.maximum(-gradient, 0).max()), float(np.abs(x * gradient).max())) / scale
        assert kkt < 1e-7 and (x >= 0).all()
        objective = float(x @ problem['G'] @ x - 2 * x @ c + target_norm2[i])
        assert objective >= -1e-8 * max(float(target_norm2[i]), 1e-15)
        total = float(x.sum()); assert total > 1e-12
        output[:, i] = x / total
        residuals.append(float(np.sqrt(max(objective, 0) / max(float(target_norm2[i]), 1e-30))))
        errors.append(kkt); amplitude_sums.append(total)
        if i in direct_cases:
            full_target = (teacher @ truth[:, i]) * weights
            direct, _ = nnls(problem['A'], full_target, maxiter=1000)
            np.testing.assert_allclose(x, direct, atol=1e-6, rtol=1e-6)
            np.testing.assert_allclose(np.sum((problem['A'] @ x - full_target)**2), max(objective, 0), atol=1e-10, rtol=1e-7)
            direct_errors.append(float(np.max(np.abs(x - direct))))
    np.testing.assert_allclose(output.sum(axis=0), 1, atol=1e-12)
    return output, dict(fits=truth.shape[1], full_objective_KKT_max=float(max(errors)), direct_full_feature_checks=len(direct_errors), direct_coefficient_max_error=float(max(direct_errors, default=0)), relative_residual_mean=float(np.mean(residuals)), relative_residual_max=float(max(residuals)), amplitude_sum_mean=float(np.mean(amplitude_sums)), amplitude_sum_min=float(min(amplitude_sums)), amplitude_sum_max=float(max(amplitude_sums)))


def target_values(fractions, cfg):
    names = list(cfg['reference_classes']); values = [fractions[i] for i in range(len(names))]
    for name, members in cfg['marginal_targets'].items():
        names.append(name); values.append(fractions[[cfg['reference_classes'].index(x) for x in members]].sum(axis=0))
    lookup = dict(zip(names, values))
    names.append('paired_R_only_E_minimum')
    values.append(np.minimum(lookup['T_anchored_reg_only'], lookup['T_effector_without_anchored_reg']))
    assert len(names) == cfg['total_declared_targets']
    return names, np.asarray(values)


def summarize(truth, estimate):
    error = estimate - truth
    centered = truth - truth.mean(); predicted_centered = estimate - estimate.mean()
    variance = float(centered @ centered)
    if variance > 1e-14 and float(predicted_centered @ predicted_centered) > 1e-14:
        pearson = float(np.corrcoef(truth, estimate)[0, 1])
        alternate = float((centered @ predicted_centered) / np.sqrt((centered @ centered) * (predicted_centered @ predicted_centered)))
        assert abs(pearson - alternate) < 1e-12
        rt, re = rankdata(truth), rankdata(estimate)
        spearman = float(np.corrcoef(rt, re)[0, 1])
    else: pearson = None; spearman = None
    positive = truth >= .01; zero = truth <= 1e-12
    return dict(mixtures=len(truth), true_mean=float(truth.mean()), estimate_mean=float(estimate.mean()), MAE=float(np.abs(error).mean()), RMSE=float(np.sqrt(np.mean(error**2))), bias=float(error.mean()), Pearson=pearson, Spearman=spearman, calibration_slope=float(centered @ predicted_centered / variance) if variance > 1e-14 else None, positive_RNA_ge001_mixtures=int(positive.sum()), positive_RNA_ge001_MAE=float(np.abs(error[positive]).mean()) if positive.any() else None, designed_zero_mixtures=int(zero.sum()), designed_zero_leakage_gt005=float((estimate[zero] > .05).mean()) if zero.any() else None, designed_zero_max_contribution=float(estimate[zero].max()) if zero.any() else None)


def fixture():
    rng = np.random.default_rng(2026102008)
    reference = rng.uniform(.01, 1, (80, 4)); reference /= reference.sum(axis=0)
    weights = 1 / np.sqrt(reference.mean(axis=1) + 1e-6)
    problem, reason = geometry(reference, weights); assert not reason
    truth = rng.dirichlet(np.ones(4), 20).T
    estimated, proof = recover(reference, reference, truth, weights, problem, set(range(20)))
    np.testing.assert_allclose(estimated, truth, atol=1e-8, rtol=1e-8)
    bad = reference.copy(); bad[:, 1] = bad[:, 0]
    rejected, reason = geometry(bad, weights); assert rejected is None and reason == 'FULL_BASIS_RANK_DEFICIENT'
    perfect = summarize(truth[0], truth[0]); assert perfect['MAE'] == 0 and abs(perfect['Pearson'] - 1) < 1e-12
    constant = summarize(np.zeros(20), np.zeros(20)); assert constant['Pearson'] is None
    return dict(status='SYNTHETIC_EXACT_MIXTURE_FULL_VS_REDUCED_NNLS_DUPLICATE_BASIS_AND_ZERO_TARGET_ACCEPTED', real_inputs_used=False)


def main():
    cfg = json.load(sys.stdin); prior = cfg['prior_reference_plan']
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    old_control = Path(os.environ['PANCANCER_PRIOR_REFERENCE_CONTROL'])
    old_output = Path(os.environ['PANCANCER_PRIOR_REFERENCE_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    assert sha(old_output/'private_terminal.json') == cfg['prior_reference_terminal_SHA256']
    parent = json.loads((old_output/'private_terminal.json').read_text())
    assert sha(old_control/'run_CM_marker_disjoint_reference_diagnostic_v1.py') == parent['source_SHA256']
    parent_receipt = json.loads((old_control/'private_launch_receipt.json').read_text())
    assert sha(old_control/'qualify_CM_lineage_state_reference_inputs_v1.py') == parent_receipt['helper_SHA256']
    p = old_output/'private_plan_before_RNA_and_model_fitting.json'
    assert sha(p) == parent['files'][p.name]['SHA256']
    prior_plan = json.loads(p.read_text()); assert prior_plan['plan'] == prior
    sys.path.insert(0, str(old_control))
    from run_CM_marker_disjoint_reference_diagnostic_v1 import select_rows, digest_rows
    from qualify_CM_lineage_state_reference_inputs_v1 import label_T
    atlas = root.parents[1]/'runs/20261013_immunotherapy_atlas_v1'
    assert sha(atlas/'manifests/immune_authority_release_public.json') == prior['release_SHA256']
    z = zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    md = pd.DataFrame({k:np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset','library_id','major']})
    md['row'] = np.arange(len(md)); md['source'] = md.dataset; md['library'] = md.dataset+'::'+md.library_id
    cp = atlas/'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    assert sha(cp) == prior['CM_material_SHA256']
    cm = pd.read_csv(cp,sep='\t',usecols=['sampleID'],keep_default_na=False)
    material = md.library.isin(cm.sampleID).to_numpy()
    old = root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    assert sha(old/'private_fixed34_gene_T_counts_and_depth.npz') == prior['T_cache_SHA256']
    cache = dict(np.load(old/'private_fixed34_gene_T_counts_and_depth.npz',allow_pickle=False))
    tr = cache['global_original_T_row_index']; np.testing.assert_array_equal(tr,np.flatnonzero(md.major.eq('T')))
    fixed = root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed) == prior['fixed_panel_contract_SHA256']
    panels = json.loads(fixed.read_text())['panels']
    md['coarse_lineage'] = 'Other_assigned'
    md.loc[md.major.isin(['T','NK']),'coarse_lineage'] = md.loc[md.major.isin(['T','NK']),'major']
    md.loc[md.major.isin(prior['boundary_majors']),'coarse_lineage'] = 'Boundary_unresolved'
    md['T_capture_state'] = ''
    md.loc[tr,'T_capture_state'] = label_T(cache['raw_gene_counts'],cache['genes'].tolist(),panels,prior['anchor_gene'],prior['effector_genes'],1)
    selections = {}
    for task in prior['tasks']:
        use = material & (md.major.eq('T').to_numpy() if task=='T_capture_state' else True)
        selections[task] = select_rows(md.loc[use],task,prior)
        assert digest_rows(selections[task]) == prior_plan['selection'][task]['native_rows_SHA256']
    coarse = selections['coarse_lineage']; coarse = coarse[md.major.iloc[coarse].ne('T').to_numpy()]
    rows = np.unique(np.r_[selections['T_capture_state'],coarse])
    frame = md.iloc[rows].copy(); frame['candidate'] = frame.coarse_lineage
    frame.loc[frame.major.eq('T'),'candidate'] = frame.loc[frame.major.eq('T'),'T_capture_state']
    assert set(frame.candidate) == set(cfg['reference_classes'])
    groups = frame.groupby(['source','library','candidate'],sort=True).size().reset_index(name='cells')
    keys = ['source','library','candidate']; gi = pd.MultiIndex.from_frame(groups[keys]).get_indexer(pd.MultiIndex.from_frame(frame[keys]))
    source_order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured = np.asarray(read_elem(z['varm/measured_in_source']),bool)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    status = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    maps_path = old/'private_source_gene_measurement_contract.json'
    manifest = json.loads((old/'private_terminal.json').read_text()); assert sha(maps_path) == manifest['files'][maps_path.name]['sha256']
    maps = json.loads(maps_path.read_text()); excluded = {int(v) for x in maps for v in x['physical_columns'] if v>=0}
    good = measured[:,[source_order.index(s) for s in prior['source_scope']]].all(axis=1)
    good &= ~np.char.startswith(canonical,'AMBIGUOUS::') & (status!='ambiguous_source_feature')
    good &= ~pd.Index(canonical).duplicated(keep=False); good[list(excluded)] = False
    common = np.flatnonzero(good); assert len(common)==prior['expected_source_common_disjoint_genes']
    assert digest_rows(common) == prior_plan['common_physical_gene_indices_SHA256']
    out.mkdir(mode=0o700)
    (out/'private_plan_before_expression_and_mixture_effects.json').write_text(json.dumps(dict(plan=cfg,created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),selected_rows_SHA256=digest_rows(rows),selected_cells=len(rows),material_class_groups=len(groups)),indent=2)+'\n')
    aggregate = np.zeros((len(groups),len(common)),np.float64); total = np.zeros(len(groups)); num = np.zeros(len(groups),np.int64)
    independent_gene_total = np.zeros(len(common)); counts=z['layers/counts']; nc,nf=map(int,counts.attrs['shape']); done=0
    for start in range(0,nc,cfg['chunk_rows']):
        stop=min(nc,start+cfg['chunk_rows']);lo=int(np.searchsorted(rows,start));hi=int(np.searchsorted(rows,stop))
        if hi==lo:continue
        ptr=np.asarray(counts['indptr'][start:stop+1],np.int64);aa,bb=int(ptr[0]),int(ptr[-1])
        x=csr_matrix((np.asarray(counts['data'][aa:bb]),np.asarray(counts['indices'][aa:bb]),ptr-ptr[0]),shape=(stop-start,nf))[rows[lo:hi]-start]
        x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
        assert np.isfinite(x.data).all() and (x.data>0).all() and (x.data==np.floor(x.data)).all()
        src=frame.source.to_numpy()[lo:hi]
        for s in set(src):assert measured[x[src==s].indices,source_order.index(s)].all()
        native=np.asarray(x.sum(axis=1)).ravel().astype(np.int64);tt=frame.major.to_numpy()[lo:hi]=='T'
        if tt.any():
            pos=np.searchsorted(tr,rows[lo:hi][tt]);np.testing.assert_array_equal(tr[pos],rows[lo:hi][tt]);np.testing.assert_array_equal(native[tt],cache['full_native_total_counts'][pos]);np.testing.assert_array_equal(np.diff(x.indptr)[tt],cache['full_native_detected_genes'][pos])
            for s in set(src[tt]):
                take=tt&(src==s);columns=next(m for m in maps if m['dataset']==s)['physical_columns'];assert min(columns)>=0
                np.testing.assert_array_equal(x[take][:,columns].toarray(),cache['raw_gene_counts'][np.searchsorted(tr,rows[lo:hi][take])])
        sub=x[:,common];local_groups,inverse=np.unique(gi[lo:hi],return_inverse=True)
        incidence=csr_matrix((np.ones(hi-lo),(inverse,np.arange(hi-lo))),shape=(len(local_groups),hi-lo))
        aggregate[local_groups]+=np.asarray((incidence@sub).toarray())
        total[local_groups]+=np.bincount(inverse,weights=native,minlength=len(local_groups));num[local_groups]+=np.bincount(inverse,minlength=len(local_groups))
        independent_gene_total+=np.asarray(sub.sum(axis=0)).ravel();done+=hi-lo
        (out/'private_progress.json').write_text(json.dumps(dict(phase='AGGREGATING_EXISTING_RNA_IN_RAM',selected_cells_read=done,total_selected_cells=len(rows)))+'\n')
    np.testing.assert_array_equal(num,groups.cells);np.testing.assert_array_equal(aggregate.sum(axis=0),independent_gene_total)
    assert done==len(rows) and (total>0).all()
    classes=cfg['reference_classes'];k=len(classes);sources=prior['source_scope'];profiles={};coverage=[]
    assert k==10 and classes[:7]==prior['T_classes']
    for source in sources:
        p=np.zeros((len(common),k))
        for ci,candidate in enumerate(classes):
            take=groups.source.eq(source)&groups.candidate.eq(candidate); n=int(num[take].sum());libs=int(take.sum());coverage.append(dict(source=source,candidate=candidate,selected_cells=n,material_units=libs,source_profile_eligible=n>=20 and libs>=2))
            if n and total[take].sum()>0:p[:,ci]=aggregate[take].sum(axis=0)/total[take].sum()
        profiles[source]=p
    contexts=[];metrics=[];calibration=[];applicability=[];proofs=[]
    for heldout in sources:
        train=[s for s in sources if s!=heldout]
        bad=[]
        training_sources={}
        for candidate in classes:
            eligible_train=[s for s in train if next(v for v in coverage if v['source']==s and v['candidate']==candidate)['source_profile_eligible']]
            training_sources[candidate]=eligible_train
            held=next(v for v in coverage if v['source']==heldout and v['candidate']==candidate)
            if len(eligible_train)<2 or not held['source_profile_eligible']:bad.append(candidate)
        if bad:
            for scaling in cfg['scalings']:contexts.append(dict(source=heldout,scaling=scaling,status='UNESTIMABLE_FULL_BASIS',reason='REQUIRED_SOURCE_CLASS_SUPPORT',unsupported_source_class_positions=len(bad)))
            continue
        reference=np.column_stack([np.mean([profiles[s][:,ci] for s in training_sources[candidate]],axis=0) for ci,candidate in enumerate(classes)]); teacher=profiles[heldout]
        assert np.isfinite(reference).all() and (reference>=0).all() and np.all(reference.sum(axis=0)<=1+1e-12)
        rng=np.random.default_rng(seed(str(cfg['seed'])+'|'+heldout))
        truth_regimes={'dense':rng.dirichlet(np.ones(k),200).T,'sparse':rng.dirichlet(np.full(k,.2),200).T}
        dilution=np.zeros((k,200));mass=rng.uniform(.02,.4,200);dilution[:7]=rng.dirichlet(np.ones(7),200).T*mass;dilution[7:]=rng.dirichlet(np.ones(3),200).T*(1-mass)
        truth_regimes['T_dilution']=dilution;truth_regimes['pure_profile_anchors']=np.eye(k)
        for scaling in cfg['scalings']:
            weights=1/np.sqrt(reference.mean(axis=1)+1e-6) if scaling=='train_mean_weighted' else np.ones(len(common))
            problem,reason=geometry(reference,weights)
            if problem is None:contexts.append(dict(source=heldout,scaling=scaling,status='UNESTIMABLE_FULL_BASIS',reason=reason));continue
            self_truth=np.c_[np.eye(k),rng.dirichlet(np.ones(k),20).T]
            self_est,self_proof=recover(reference,reference,self_truth,weights,problem,{0,k})
            np.testing.assert_allclose(self_est,self_truth,atol=1e-6,rtol=1e-6)
            estimates={}; truth_values={}; numerical=[self_proof]
            for regime,truth in truth_regimes.items():
                checks=set(range(k)) if regime=='pure_profile_anchors' else {0}
                estimate,proof=recover(reference,teacher,truth,weights,problem,checks);numerical.append(proof)
                names,actual=target_values(truth,cfg);_,predicted=target_values(estimate,cfg);estimates[regime]=predicted;truth_values[regime]=actual
            truth_values['combined_mixed']=np.concatenate([truth_values[n] for n in ['dense','sparse','T_dilution']],axis=1)
            estimates['combined_mixed']=np.concatenate([estimates[n] for n in ['dense','sparse','T_dilution']],axis=1)
            lib_true=[];lib_pred=[];positive_libs=np.zeros(len(names),np.int64);eligible_libraries=0
            for library,parts in groups[groups.source.eq(heldout)].groupby('library',sort=True):
                usable=parts[parts.cells>=5];available=[classes.index(c) for c in usable.candidate]
                if len(available)<2:continue
                ids=usable.index.to_numpy();b=(aggregate[ids]/total[ids,None]).T
                lrng=np.random.default_rng(seed(str(cfg['seed'])+'|'+heldout+'|'+library))
                local_truth=lrng.dirichlet(np.ones(len(available)),20).T
                estimate,proof=recover(reference,b,local_truth,weights,problem,{0} if eligible_libraries==0 else set());numerical.append(proof)
                full=np.zeros((k,20));full[available]=local_truth
                _,actual=target_values(full,cfg);_,predicted=target_values(estimate,cfg);lib_true.append(actual);lib_pred.append(predicted);positive_libs+=(actual>=.01).any(axis=1);eligible_libraries+=1
            if lib_true:
                truth_values['library_stress']=np.concatenate(lib_true,axis=1);estimates['library_stress']=np.concatenate(lib_pred,axis=1)
            for regime,actual in truth_values.items():
                for ti,target in enumerate(names):
                    item=dict(source=heldout,scaling=scaling,regime=regime,target=target,**summarize(actual[ti],estimates[regime][ti]))
                    if regime=='library_stress':item.update(material_units=eligible_libraries,positive_RNA_ge001_material_units=int(positive_libs[ti]))
                    metrics.append(item)
            for ti,target in enumerate(names):
                combined=summarize(truth_values['combined_mixed'][ti],estimates['combined_mixed'][ti]);anchor=summarize(truth_values['pure_profile_anchors'][ti],estimates['pure_profile_anchors'][ti])
                eligible=(combined['positive_RNA_ge001_mixtures']>=30 and combined['MAE']<=.05 and combined['Pearson'] is not None and combined['Pearson']>=.7 and anchor['designed_zero_max_contribution'] is not None and anchor['designed_zero_max_contribution']<=.1)
                applicability.append(dict(source=heldout,scaling=scaling,target=target,status='EXPLORATORY_MEASUREMENT_CANDIDATE' if eligible else 'RECOVERY_LIMITATION',primary_scaling=scaling==cfg['primary_scaling'],MAE=combined['MAE'],Pearson=combined['Pearson'],pure_anchor_zero_max=anchor['designed_zero_max_contribution'],not_clinical_or_identity_validation=True))
                actual=truth_values['combined_mixed'][ti];predicted=estimates['combined_mixed'][ti];bin_id=np.searchsorted(cfg['calibration_bins'],actual,side='right')-1
                for bi in range(len(cfg['calibration_bins'])-1):
                    take=bin_id==bi;calibration.append(dict(source=heldout,scaling=scaling,target=target,bin=bi,mixtures=int(take.sum()),true_mean=float(actual[take].mean()) if take.any() else None,estimate_mean=float(predicted[take].mean()) if take.any() else None))
            contexts.append(dict(source=heldout,scaling=scaling,status='COMPLETE_KNOWN_RNA_MIXTURE_CALIBRATION',full_basis_classes=k,genes=len(common),basis_condition=problem['condition'],heldout_pooled_mixtures=610,library_stress_material_units=eligible_libraries,library_stress_mixtures=eligible_libraries*20))
            proofs.append(dict(source=heldout,scaling=scaling,total_fits=sum(v['fits'] for v in numerical),full_objective_KKT_max=max(v['full_objective_KKT_max'] for v in numerical),direct_full_feature_checks=sum(v['direct_full_feature_checks'] for v in numerical),direct_coefficient_max_error=max(v['direct_coefficient_max_error'] for v in numerical),mean_relative_residual=float(np.average([v['relative_residual_mean'] for v in numerical],weights=[v['fits'] for v in numerical])),min_amplitude_sum=min(v['amplitude_sum_min'] for v in numerical),max_amplitude_sum=max(v['amplitude_sum_max'] for v in numerical),self_reference_exact_mixture_recovery=True))
            (out/'private_progress.json').write_text(json.dumps(dict(phase='FITTING_KNOWN_RNA_MIXTURES',contexts_completed=len(contexts),registered_contexts=8))+'\n')
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    assert len(contexts)==8
    save('private_all8_mixture_contexts.json',contexts);save('private_all_source_class_profile_coverage.json',coverage);save('private_all_source_target_recovery_metrics.json',metrics);save('private_all_source_target_applicability.json',applicability);save('private_all_source_target_binned_calibration.json',calibration);save('private_full_NNLS_numeric_proofs.json',proofs)
    files={p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values())<cfg['output_cap_bytes']
    save('private_terminal.json',dict(status='COMPLETE_ALL_SOURCE_KNOWN_RNA_MIXTURE_CALIBRATION_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),files=files,all_real_RNA_profiles_weights_truth_predictions_random_arrays_RAM_only=True,original_project_read_only=True,no_new_clinical_P_FDR_or_formal_identity=True))
    print(json.dumps(dict(status='COMPLETE_MIXTURE_REVIEW_PENDING',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
