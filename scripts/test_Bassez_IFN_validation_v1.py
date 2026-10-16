"""Three prespecified independent-cohort questions and full sensitivity families."""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

ALPHA = 'HALLMARK_INTERFERON_ALPHA_RESPONSE'
GAMMA = 'HALLMARK_INTERFERON_GAMMA_RESPONSE'
CLINICAL = ['My_logRNA', 'My_loggenes', 'My_logcells']
COUPLING = ['T_logRNA', 'T_loggenes', 'My_logRNA', 'My_loggenes', 'T_logcells', 'My_logcells', 'My_IFNG']


def seed(key):
    return 20261016 + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def rankz(x):
    r = rankdata(np.asarray(x, float))
    if r.std() < 1e-12:
        raise ValueError('constant continuous variable')
    return (r-r.mean())/r.std()


def nuisance(frame, numeric, levels, adjusted=True):
    cols = [np.ones(len(frame))]
    if adjusted:
        for field, categories in levels.items():
            cols.extend(frame[field].eq(category).to_numpy(float) for category in categories[1:])
        cols.extend(rankz(frame[field].to_numpy()) for field in numeric)
    return np.column_stack(cols)


def fitted(frame, question, levels, adjusted):
    y = rankz(frame[question['y']])
    x = frame.expanded.to_numpy(float) if question['binary'] else rankz(frame[question['x']])
    N = nuisance(frame, question['numeric'], levels, adjusted)
    D = np.column_stack((N, x))
    df = len(frame)-D.shape[1]
    if len(frame) < 20 or df < 10 or np.linalg.matrix_rank(D) < D.shape[1]:
        raise ValueError('sample coverage, residual df or full design rank')
    solver = np.linalg.solve(D.T@D, D.T)
    beta = solver@y
    residual = y-D@beta
    variance = np.linalg.inv(D.T@D)[-1, -1]
    se = np.sqrt(max(float(residual@residual)/df*variance, 0))
    if se < 1e-12:
        raise ValueError('zero residual variance')
    return float(beta[-1]), float(beta[-1]/se), N, D, y, solver, variance


def model_test(frame, question, levels, adjusted, key):
    b, t, N, D, y, solver, variance = fitted(frame, question, levels, adjusted)
    df = len(frame)-D.shape[1]
    reduced = N@np.linalg.lstsq(N, y, rcond=None)[0]
    residual = y-reduced
    rng = np.random.default_rng(seed(key))
    extreme = 0
    for start in range(0, 20000, 250):
        perm = np.column_stack([reduced+residual[rng.permutation(len(frame))] for _ in range(min(250, 20000-start))])
        pb = solver@perm
        errors = perm-D@pb
        se = np.sqrt(np.maximum(np.sum(errors*errors, axis=0)/df*variance, 0))
        pt = np.divide(pb[-1], se, out=np.full(len(se), np.inf), where=se > 1e-12)
        extreme += np.count_nonzero(np.abs(pt) >= abs(t)-1e-12)
    bootstrap = []
    group_indices = [np.flatnonzero(frame.expanded.to_numpy() == value) for value in [0, 1]]
    for _ in range(2000):
        ix = np.concatenate([rng.choice(group, len(group), replace=True) for group in group_indices])
        try:
            bootstrap.append(fitted(frame.iloc[ix], question, levels, adjusted)[0])
        except (ValueError, np.linalg.LinAlgError):
            pass
    lower, upper = np.quantile(bootstrap, [.025, .975]) if len(bootstrap) >= 1800 else (np.nan, np.nan)
    return dict(status='estimated_exploratory_partial_function', beta=b, t=t, p=(extreme+1)/20001,
        lower=float(lower), upper=float(upper), valid_bootstrap=len(bootstrap), residual_df=df)


def bh(p):
    p = np.asarray(p, float); order = np.argsort(p); q = np.empty(len(p))
    q[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1])
    return q


def stratified_rb(values, labels, strata):
    total = 0.; weight = 0
    for group in sorted(set(strata)):
        take = strata == group; x = values[take]; y = labels[take]
        n1 = int(y.sum()); n0 = len(y)-n1
        if min(n1, n0) == 0:
            raise ValueError('treatment stratum lacks either expansion category')
        U = rankdata(x)[y == 1].sum()-n1*(n1+1)/2
        total += len(y)*(2*U/(n1*n0)-1); weight += len(y)
    return total/weight


def wilcoxon(frame, column, key):
    x = frame[column].to_numpy(float); y = frame.expanded.to_numpy(int); strata = frame.treatment.to_numpy(str)
    if min(int(y.sum()), len(y)-int(y.sum())) < 5:
        raise ValueError('expansion category coverage insufficient')
    effect = stratified_rb(x, y, strata)
    rng = np.random.default_rng(seed(key)); extreme = 0
    groups = [np.flatnonzero(strata == group) for group in sorted(set(strata))]
    for _ in range(20000):
        label = y.copy()
        for group in groups:
            label[group] = rng.permutation(y[group])
        extreme += abs(stratified_rb(x, label, strata)) >= abs(effect)-1e-12
    boot = []
    # Treatment x expansion groups retain the observed comparison in each stratum.
    blocks = [np.flatnonzero((strata == group) & (y == value)) for group in sorted(set(strata)) for value in [0, 1]]
    for _ in range(2000):
        ix = np.concatenate([rng.choice(group, len(group), replace=True) for group in blocks])
        boot.append(stratified_rb(x[ix], y[ix], strata[ix]))
    lo, hi = np.quantile(boot, [.025, .975])
    return dict(status='estimated_treatment_stratified_Wilcoxon_sensitivity', beta=effect, p=(1+extreme)/20001,
        lower=float(lo), upper=float(hi), valid_bootstrap=len(boot), residual_df=np.nan)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs = root/'Bassez_IFN_validation_preparation_20261016_v1'
    out = root/'Bassez_IFN_validation_results_20261016_v1'
    assert not out.exists()
    terminal = json.loads((inputs/'private_terminal.json').read_text())
    assert terminal['status'] == 'BASSEZ_FIXED_IFN_SCORES_PREPARED_NO_EFFECTS'
    frame = pd.read_csv(inputs/'private_patient_scores_and_covariates.tsv', sep='\t')
    assert hashlib.sha256((inputs/'private_patient_scores_and_covariates.tsv').read_bytes()).hexdigest() == terminal['patient_table_sha256']
    assert frame.patient_token.is_unique and len(frame) == terminal['patient_blocks']
    questions = [
        dict(name='myeloid_IFNa_vs_T_clonal_expansion', y='My_'+ALPHA, x='expanded', binary=True, numeric=CLINICAL),
        dict(name='myeloid_IFNg_vs_T_clonal_expansion', y='My_'+GAMMA, x='expanded', binary=True, numeric=CLINICAL),
        dict(name='T_IFNG_vs_myeloid_IFNg', y='My_'+GAMMA, x='T_IFNG', binary=False, numeric=COUPLING),
    ]
    levels = {field: sorted(frame[field].unique().tolist()) for field in ['treatment', 'cancer_type']}
    out.mkdir()
    execution = dict(created_UTC=datetime.now(timezone.utc).isoformat(), questions=questions, categorical_levels=levels,
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        preparation_plan_sha256=hashlib.sha256((inputs/'plan_before_scores_or_effects.json').read_bytes()).hexdigest(),
        main_family=3, unadjusted_family=3, Wilcoxon_family=2,
        Wilcoxon_precision='Treatment-stratum weighted rank-biserial statistic centered at zero under within-stratum label permutation; patient bootstrap within treatment x expansion.',
        permutation_seed='20261016 + SHA256(model|question) first8 hex', significance_tests_before_this_execution=False)
    (out/'execution_contract_before_effects.json').write_text(json.dumps(execution, indent=2)+'\n')
    rows = []
    interface_ok = all(terminal['interface_source_positive_confirmed'].values())
    for model in ['primary_adjusted', 'sensitivity_unadjusted']:
        for question in questions:
            row = dict(model=model, question=question['name'], n_patients=len(frame),
                effect_definition='rank-SD difference E-minus-NE' if question['binary'] else 'standardized-rank regression slope',
                status='missing_measurement', beta=np.nan, p=np.nan, lower=np.nan, upper=np.nan, valid_bootstrap=0)
            if question['y'] in frame and question['x'] in frame and (question['binary'] or interface_ok):
                try:
                    row.update(model_test(frame, question, levels, model == 'primary_adjusted', model+'|'+question['name']))
                except (ValueError, np.linalg.LinAlgError) as exc:
                    row['status'] = 'unestimable_'+str(exc)
            rows.append(row)
    for question in questions[:2]:
        row = dict(model='sensitivity_Wilcoxon', question=question['name'], n_patients=len(frame),
            effect_definition='treatment-stratum size weighted rank-biserial E-minus-NE',
            status='missing_measurement', beta=np.nan, p=np.nan, lower=np.nan, upper=np.nan, valid_bootstrap=0)
        if question['y'] in frame:
            try:
                row.update(wilcoxon(frame, question['y'], 'sensitivity_Wilcoxon|'+question['name']))
            except ValueError as exc:
                row['status'] = 'unestimable_'+str(exc)
        rows.append(row)
    result = pd.DataFrame(rows)
    result['BH'] = np.nan
    for _, group in result.groupby('model'):
        q = bh(group.p.fillna(1).to_numpy())
        result.loc[group.index, 'BH'] = np.where(group.p.notna(), q, np.nan)
    assert len(result) == 8 and not result.duplicated(['model', 'question']).any()
    result.to_csv(out/'all8_private_IFN_validation_statistics.tsv', sep='\t', index=False, na_rep='NA')
    descriptive = dict(interface_joint_receptor_range=terminal['myeloid_joint_receptor_fraction_range'],
        two_myeloid_IFN_readout_rho=float(spearmanr(frame['My_'+ALPHA], frame['My_'+GAMMA]).statistic) if all('My_'+r in frame for r in [ALPHA, GAMMA]) else None,
        note='Descriptive expression compatibility and correlated readouts, not direct binding, mediation or clinical prediction.')
    (out/'private_descriptive_interface_and_readout.json').write_text(json.dumps(descriptive, indent=2)+'\n')
    report = dict(status='BASSEZ_THREE_FIXED_IFN_QUESTIONS_AND_WILCOXON_COMPLETE_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(), records=len(result),
        result_sha256=hashlib.sha256((out/'all8_private_IFN_validation_statistics.tsv').read_bytes()).hexdigest(),
        direct_interaction_verified=False, pathological_response_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report)); print(result.to_json(orient='records')); print(json.dumps(descriptive))


if __name__ == '__main__':
    main()
