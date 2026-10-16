"""Independent manual-rank/QR and Mann-Whitney replay of all fixed tests."""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu


def ranks(x):
    x = np.asarray(x, float)
    return np.array([np.count_nonzero(x < v)+(np.count_nonzero(x == v)+1)/2 for v in x], float)


def rz(x):
    r = ranks(x)
    if r.std() < 1e-12:
        raise ValueError('constant')
    return (r-r.mean())/r.std()


def design(f, question, levels, adjusted):
    cols = [np.ones(len(f))]
    if adjusted:
        for name, values in levels.items():
            cols.extend(f[name].eq(value).to_numpy(float) for value in values[1:])
        cols.extend(rz(f[c]) for c in question['numeric'])
    return np.column_stack(cols)


def fit(f, question, levels, adjusted):
    N = design(f, question, levels, adjusted)
    x = f.expanded.to_numpy(float) if question['binary'] else rz(f[question['x']])
    y = rz(f[question['y']]); D = np.column_stack((N, x)); df = len(f)-D.shape[1]
    if len(f) < 20 or df < 10 or np.linalg.matrix_rank(D) != D.shape[1]:
        raise ValueError('rank/coverage')
    Q, R = np.linalg.qr(D, mode='reduced')
    coef = np.linalg.solve(R, Q.T@y); error = y-D@coef
    variance = np.linalg.inv(D.T@D)[-1, -1]
    se = np.sqrt(float(error@error)/df*variance)
    if se < 1e-12:
        raise ValueError('variance')
    return coef[-1], coef[-1]/se, N, D, y, Q, R, variance


def close(a, b):
    assert np.isclose(a, b, rtol=2e-7, atol=1e-9), (a, b)


def rb(x, y, strata):
    result = 0.
    for label in sorted(set(strata)):
        k = strata == label; a = x[k & (y == 1)]; b = x[k & (y == 0)]
        if min(len(a), len(b)) == 0:
            raise ValueError('missing comparison')
        # Pairwise signs, separately checked against SciPy U.
        effect = float(np.sign(a[:, None]-b[None, :]).mean())
        close(effect, 2*mannwhitneyu(a, b, alternative='two-sided').statistic/(len(a)*len(b))-1)
        result += np.count_nonzero(k)*effect
    return result/len(x)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inp = root/'Bassez_IFN_validation_preparation_20261016_v1'
    out = root/'Bassez_IFN_validation_results_20261016_v1'
    prep = json.loads((inp/'private_terminal.json').read_text())
    frame = pd.read_csv(inp/'private_patient_scores_and_covariates.tsv', sep='\t')
    assert hashlib.sha256((inp/'private_patient_scores_and_covariates.tsv').read_bytes()).hexdigest() == prep['patient_table_sha256']
    plan = json.loads((out/'execution_contract_before_effects.json').read_text())
    results = pd.read_csv(out/'all8_private_IFN_validation_statistics.tsv', sep='\t')
    assert len(results) == 8
    questions = {q['name']: q for q in plan['questions']}
    for row in results.itertuples():
        q = questions[row.question]
        rng = np.random.default_rng(20261016+int(hashlib.sha256((row.model+'|'+row.question).encode()).hexdigest()[:8], 16))
        if row.model == 'sensitivity_Wilcoxon':
            x = frame[q['y']].to_numpy(float); y = frame.expanded.to_numpy(int); s = frame.treatment.to_numpy(str)
            effect = rb(x, y, s); close(effect, row.beta)
            groups = [np.flatnonzero(s == value) for value in sorted(set(s))]
            extreme = 0
            for _ in range(20000):
                labels = y.copy()
                for group in groups:
                    labels[group] = rng.permutation(y[group])
                extreme += abs(rb(x, labels, s)) >= abs(effect)-1e-12
            close((1+extreme)/20001, row.p)
            groups = [np.flatnonzero((s == value) & (y == cat)) for value in sorted(set(s)) for cat in [0, 1]]
            samples = []
            for _ in range(2000):
                ix = np.concatenate([rng.choice(g, len(g), replace=True) for g in groups])
                samples.append(rb(x[ix], y[ix], s[ix]))
            lo, hi = np.quantile(samples, [.025, .975]); close(lo, row.lower); close(hi, row.upper)
            assert row.valid_bootstrap == 2000
            continue
        adjusted = row.model == 'primary_adjusted'
        try:
            b, t, N, D, y, Q, R, variance = fit(frame, q, plan['categorical_levels'], adjusted)
        except ValueError:
            assert pd.isna(row.p) and row.status.startswith('unestimable_')
            continue
        close(b, row.beta); close(t, row.t)
        df = len(frame)-D.shape[1]; assert df == row.residual_df
        reduced = N@np.linalg.lstsq(N, y, rcond=None)[0]; residual = y-reduced
        permutations = np.column_stack([reduced+residual[rng.permutation(len(frame))] for _ in range(20000)])
        coef = np.linalg.solve(R, Q.T@permutations); error = permutations-D@coef
        se = np.sqrt(np.sum(error*error, axis=0)/df*variance)
        perm_t = np.divide(coef[-1], se, out=np.full(20000, np.inf), where=se > 1e-12)
        close((1+np.count_nonzero(abs(perm_t) >= abs(t)-1e-12))/20001, row.p)
        blocks = [np.flatnonzero(frame.expanded.to_numpy() == value) for value in [0, 1]]
        samples = []
        for _ in range(2000):
            ix = np.concatenate([rng.choice(g, len(g), replace=True) for g in blocks])
            try:
                samples.append(fit(frame.iloc[ix], q, plan['categorical_levels'], adjusted)[0])
            except (ValueError, np.linalg.LinAlgError):
                pass
        assert len(samples) == row.valid_bootstrap
        if len(samples) >= 1800:
            lo, hi = np.quantile(samples, [.025, .975]); close(lo, row.lower); close(hi, row.upper)
        else:
            assert pd.isna(row.lower) and pd.isna(row.upper)
    for model, group in results.groupby('model'):
        assert len(group) == (2 if model == 'sensitivity_Wilcoxon' else 3)
        p = group.p.fillna(1).to_numpy()
        expected = [min(1, min(value*len(p)/(p <= value).sum() for value in p if value >= own)) for own in p]
        for observed, calculated, raw in zip(group.BH, expected, group.p):
            if pd.isna(raw):
                assert pd.isna(observed)
            else:
                close(observed, calculated)
    report = dict(status='PASS_ALL8_BASSEZ_IFN_EFFECTS_PERMUTATIONS_BOOTSTRAPS_AND_BH',
        finished_UTC=datetime.now(timezone.utc).isoformat(),
        result_sha256=hashlib.sha256((out/'all8_private_IFN_validation_statistics.tsv').read_bytes()).hexdigest(),
        point_engine='independent manual ranks / full QR; Wilcoxon pairwise signs and SciPy U',
        old_MCP_weights_changed=False, direct_interaction_verified=False)
    (out/'private_independent_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
