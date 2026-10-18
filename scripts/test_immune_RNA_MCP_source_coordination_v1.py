"""All six frozen heldout-source tests, retaining missing and failed designs."""
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def fit_rank_partial(x, y, controls):
    n = len(x)
    if n < 20:
        return None, 'FEWER_THAN20_COMPLETE_UNITS'
    if not np.isfinite(np.c_[x, y, controls]).all():
        return None, 'NONFINITE_COMPLETE_DESIGN'
    if min(np.ptp(x), np.ptp(y)) == 0:
        return None, 'CONSTANT_SCORE'
    ranked = rankdata(controls, axis=0)
    sd = ranked.std(axis=0, ddof=1)
    if np.any(sd == 0):
        return None, 'CONSTANT_REQUIRED_COVARIATE_NO_AUTOMATIC_DELETION'
    D = np.c_[np.ones(n), (ranked - ranked.mean(axis=0)) / sd]
    if np.linalg.matrix_rank(D) != D.shape[1]:
        return None, 'COLLINEAR_REQUIRED_DESIGN_NO_AUTOMATIC_DELETION'
    if not np.linalg.cond(D) < 1e8:
        return None, 'ILL_CONDITIONED_REQUIRED_DESIGN'
    if n - D.shape[1] - 1 < 10:
        return None, 'FEWER_THAN10_FULL_MODEL_RESIDUAL_DF'
    Q = np.linalg.qr(D, mode='reduced')[0]
    rx0 = rankdata(x); ry0 = rankdata(y)
    rx = rx0 - Q @ (Q.T @ rx0); ry = ry0 - Q @ (Q.T @ ry0)
    norm = np.linalg.norm(rx) * np.linalg.norm(ry)
    if norm < 1e-12:
        return None, 'SCORE_EXPLAINED_BY_REQUIRED_COVARIATES'
    return dict(effect=float(rx @ ry / norm), design=D, Q=Q, rx=rx, ry=ry,
                condition=float(np.linalg.cond(D)), residual_df=n-D.shape[1]-1), 'ESTIMABLE'


def adjust_family(p, by=False):
    p = np.asarray(p, float); order = np.argsort(p, kind='stable'); n = len(p)
    scale = np.sum(1 / np.arange(1, n + 1)) if by else 1.
    ranked = p[order] * n * scale / np.arange(1, n + 1)
    values = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    result = np.empty(n); result[order] = values
    return result


def resample(x, y, C, fitted, slot, seed):
    n = len(x)
    rng = np.random.default_rng(np.random.SeedSequence([seed, slot, 0]))
    permutation = np.stack([rng.permutation(n) for _ in range(19999)]).astype(np.int16)
    Q = fitted['Q']; rx = fitted['rx']; ry = fitted['ry']
    permuted_rho = np.empty(19999)
    for start in range(0, 19999, 1024):
        values = ry[permutation[start:start+1024]].T
        values -= Q @ (Q.T @ values)
        denominator = np.linalg.norm(rx) * np.linalg.norm(values, axis=0)
        assert np.all(denominator > 1e-12)
        permuted_rho[start:start+len(denominator)] = rx @ values / denominator
    assert np.isfinite(permuted_rho).all()
    p = (1 + np.count_nonzero(np.abs(permuted_rho) >= abs(fitted['effect']) - 1e-12)) / 20000
    rng = np.random.default_rng(np.random.SeedSequence([seed, slot, 1]))
    bootstrap = rng.integers(0, n, size=(2000, n), dtype=np.int16)
    rho = np.full(2000, np.nan); reasons = []
    for i, take in enumerate(bootstrap):
        fit, reason = fit_rank_partial(x[take], y[take], C[take])
        reasons.append(reason)
        if fit is not None:
            rho[i] = fit['effect']
    valid = np.isfinite(rho)
    ci = np.quantile(rho[valid], [.025, .975]) if valid.sum() >= 1000 else [np.nan, np.nan]
    return dict(permutation_p=float(p), lower=float(ci[0]), upper=float(ci[1]),
        valid_bootstrap=int(valid.sum()), bootstrap_reasons=dict(Counter(reasons))), dict(
        permutation_indices=permutation, permutation_effects=permuted_rho,
        bootstrap_indices=bootstrap, bootstrap_effects=rho, bootstrap_reasons=np.asarray(reasons))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); budget = resource(root, 64 * 1024**2)
    contract_path = root / 'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json'
    assert sha(contract_path) == '1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09'
    contract = json.loads(contract_path.read_text())
    audit_path = root / 'immune_RNA_MCP_frozen_projection_audit_20261018_v1/private_terminal.json'
    audit = json.loads(audit_path.read_text())
    assert audit['status'] == 'ALL_PLANNED_SLOTS_AND_ALL_SUPPORTED_PROJECTION_VALUES_ACCEPTED'
    projection = root / 'immune_RNA_MCP_frozen_projection_20261018_v1'
    assert sha(projection / 'private_terminal.json') == audit['producer_terminal_sha256']
    scores_path = projection / 'private_R_projection/private_all_case_mean_scores.tsv'
    columns = ['fit_id', 'view', 'program', 'partition', 'sample', 'dataset', 'score']
    scores = pd.read_csv(scores_path, sep='\t', dtype={'sample': str}) if scores_path.exists() else pd.DataFrame(columns=columns)
    scores['score'] = scores.score.astype(float)
    scores = scores[scores.fit_id.eq(contract['anchor'])].copy()
    original = root / 'immune_RNA_MCP_discovery_scope_20261018_v1/private_metadata_only_case_source_partition.tsv'
    scope = pd.read_csv(original, sep='\t', keep_default_na=False)
    scope = scope[~scope.partition.eq('small_source_projection_only')].copy()
    assert scope.fit_sample.is_unique
    reviewed = scope.set_index('fit_sample').reviewed_patient
    assert set(reviewed.unique()) <= {True, False}
    metadata = []
    for view in contract['view_to_official_name']:
        for part in ['discovery', 'TNBC_holdout']:
            filename = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs' / (view+'_'+part+'_metadata.tsv')
            md = pd.read_csv(filename, sep='\t', dtype={'sample': str})
            for sample, g in md.groupby('sample', sort=True):
                assert g.dataset.nunique() == 1 and g.fraction.nunique() == 1 and g.logN.nunique() == 1
                metadata.append(dict(sample=sample,dataset=g.dataset.iloc[0],view=view,partition=part,
                    fraction=float(g.fraction.iloc[0]),logN=float(g.logN.iloc[0]),
                    cellQ=float(g.cellQ.mean()),geneQ=float(g.geneQ.mean()),selected_cells=len(g),
                    reviewed_patient=bool(reviewed.loc[sample])))
    metadata = pd.DataFrame(metadata)
    assert len(metadata) == 198 and not metadata.duplicated(['sample','view']).any()
    output = root / 'immune_RNA_MCP_source_coordination_20261018_v1'
    assert not output.exists(); output.mkdir(mode=0o700)
    (output / 'plan_before_all_frozen_source_tests.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),budget=budget,
        contract_sha256=sha(contract_path),projection_acceptance_sha256=sha(audit_path),
        main_slots=6,reviewed_sensitivity_slots=6,discovery_description_slots=12,
        fixed_anchor=contract['anchor'],residual_df='n - nuisance columns including intercept - one score predictor',
        independent_patient_claim=False,independent_original_CM_validation=False),indent=2)+'\n')
    metadata.to_csv(output / 'private_all_case_view_controls.tsv',sep='\t',index=False)
    rows=[]; sens=[]; discovery=[]; covariates=['fraction','logN','cellQ','geneQ']
    for index, slot in enumerate(contract['primary_source_association']['slots']):
        a,b,program=slot['view_A'],slot['view_B'],slot['program']
        sa=scores[scores.view.eq(a)&scores.program.eq(program)]
        sb=scores[scores.view.eq(b)&scores.program.eq(program)]
        ma=metadata[metadata.view.eq(a)];mb=metadata[metadata.view.eq(b)]
        joined=ma.merge(mb,on=['sample','dataset','partition','reviewed_patient'],suffixes=('_A','_B'),validate='one_to_one')
        joined=joined.merge(sa[['sample','score']].rename(columns={'score':'score_A'}),on='sample',how='left',validate='one_to_one')
        joined=joined.merge(sb[['sample','score']].rename(columns={'score':'score_B'}),on='sample',how='left',validate='one_to_one')
        control_names=[c+'_'+side for side in ['A','B'] for c in covariates]
        hold=joined[joined.partition.eq('TNBC_holdout')].copy()
        assert len(hold)==28 and hold.dataset.eq('GSE305078').all()
        for kind, table in [('main_sample_cases',hold),('reviewed_patient_only',hold[hold.reviewed_patient])]:
            complete=table[['score_A','score_B']+control_names].apply(np.isfinite).all(axis=1)
            data=table[complete];n=len(data)
            x=data.score_A.to_numpy();y=data.score_B.to_numpy();C=data[control_names].to_numpy(float)
            fit,reason=fit_rank_partial(x,y,C)
            row=dict(slot=index,program=program,view_A=a,view_B=b,dataset='GSE305078',partition=kind,
                eligible_units=len(table),complete_units=n,missing_score_or_covariate_units=int((~complete).sum()),
                reviewed_patient_units=int(data.reviewed_patient.sum()),status=reason,
                unadjusted_rho=float(spearmanr(x,y).statistic) if n>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan,
                adjusted_rho=np.nan,lower=np.nan,upper=np.nan,permutation_p=np.nan,
                valid_bootstrap=0,design_condition=np.nan,residual_df=np.nan,
                biological_or_independent_patient_validation=False)
            if fit is not None:
                row.update(adjusted_rho=fit['effect'],design_condition=fit['condition'],residual_df=fit['residual_df'])
                if kind=='main_sample_cases':
                    result,draws=resample(x,y,C,fit,index,contract['primary_source_association']['seed'])
                    row.update(result)
                    np.savez_compressed(output/('private_slot_%02d_resampling.npz'%index),
                        samples=data['sample'].to_numpy(str),x=x,y=y,controls=C,control_names=np.asarray(control_names),**draws)
            (rows if kind=='main_sample_cases' else sens).append(row)
        for source,g in joined[joined.partition.eq('discovery')].groupby('dataset',sort=True):
            complete=g[['score_A','score_B']+control_names].apply(np.isfinite).all(axis=1);d=g[complete]
            x=d.score_A.to_numpy();y=d.score_B.to_numpy();C=d[control_names].to_numpy(float)
            fit,reason=fit_rank_partial(x,y,C)
            rho=float(spearmanr(x,y).statistic) if len(d)>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan
            loo=[]
            for i in range(len(d)):
                keep=np.arange(len(d))!=i
                if keep.sum()>=3 and min(np.ptp(x[keep]),np.ptp(y[keep]))>0:
                    loo.append(float(spearmanr(x[keep],y[keep]).statistic))
            discovery.append(dict(slot=index,program=program,view_A=a,view_B=b,dataset=source,
                eligible_units=len(g),complete_units=len(d),unadjusted_rho=rho,
                adjusted_rho=fit['effect'] if fit is not None else np.nan,adjusted_status=reason,
                leave_one_unadjusted_min=min(loo) if loo else np.nan,leave_one_unadjusted_max=max(loo) if loo else np.nan,
                fitted_discovery_data_not_validation=True))
    main=pd.DataFrame(rows);assert len(main)==6 and len(sens)==6 and len(discovery)==12
    main['family_adjustment_p']=main.permutation_p.fillna(1)
    main['BH6']=adjust_family(main.family_adjustment_p);main['BY6']=adjust_family(main.family_adjustment_p,by=True)
    main.loc[main.permutation_p.isna(),['BH6','BY6']]=np.nan
    main.to_csv(output/'private_all6_primary_source_tests.tsv',sep='\t',index=False)
    pd.DataFrame(sens).to_csv(output/'private_all6_reviewed_patient_sensitivity.tsv',sep='\t',index=False)
    pd.DataFrame(discovery).to_csv(output/'private_all12_fitted_source_descriptions.tsv',sep='\t',index=False)
    files={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(
        status='ALL_FROZEN_SOURCE_SLOTS_COMPLETE_INDEPENDENT_RANDOM_STREAM_ACCEPTANCE_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        primary_slots=6,reviewed_sensitivity_slots=6,discovery_descriptive_slots=12,
        actual_estimable_main_slots=int(main.status.eq('ESTIMABLE').sum()),
        no_clinical_or_direct_communication_validation=True),indent=2)+'\n')
    print('ALL_SOURCE_SLOTS_DONE_INDEPENDENT_ACCEPTANCE_PENDING')


if __name__=='__main__':
    main()
