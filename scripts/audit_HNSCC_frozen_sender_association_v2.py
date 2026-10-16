"""Independent full-QR replay of every primary randomization and descriptor."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def ranks(matrix):
    return (matrix[:, None, :] > matrix[None, :, :]).sum(axis=1)+((matrix[:, None, :] == matrix[None, :, :]).sum(axis=1)+1)/2


def fit(data, require_t=True):
    n = len(data)
    if n < 20 or n-14 < 10 or not np.isfinite(data).all():
        raise ValueError('n_df_or_missing')
    rr = ranks(np.column_stack((data[:, :2], data[:, 4:])))
    sd = np.std(rr, axis=0)
    if np.min(sd) < 1e-12:
        raise ValueError('constant')
    rr = (rr-rr.mean(axis=0))/sd
    D = np.column_stack((np.ones(n), data[:, 2:4], rr[:, 2:], rr[:, 0]))
    singular = np.linalg.svd(D, compute_uv=False)
    if np.count_nonzero(singular > singular[0]*max(D.shape)*np.finfo(float).eps) != 14 or singular[0]/singular[-1] > 1e8:
        raise ValueError('rank_condition')
    Q, R = np.linalg.qr(D, mode='reduced')
    if float(R[-1, -1]**2) < 1e-12:
        raise ValueError('zero_residual_exposure')
    solver = np.linalg.solve(R, Q.T)
    y = rr[:, 1]; b = np.linalg.solve(R, Q.T@y)
    ri = np.linalg.solve(R, np.eye(14))
    variance = float(ri[-1]@ri[-1]); df = n-14
    if not require_t:
        return float(b[-1]), np.nan, D, solver, variance, y, df
    residual = y-Q@(Q.T@y)
    se = np.sqrt(float(residual@residual)/df*variance)
    if se < 1e-12:
        raise ValueError('residual_variance')
    return float(b[-1]), float(b[-1]/se), D, solver, variance, y, df


def same(a, b):
    if not np.isfinite(a):
        assert pd.isna(b), (a, b)
    else:
        assert np.isclose(a, b, atol=2e-8, rtol=2e-6), (a, b)


def descriptor(data, expected):
    try:
        effect = fit(data, require_t=False)[0]
    except (ValueError, np.linalg.LinAlgError):
        assert pd.isna(expected.beta) and expected.LOO_valid == 0
        return
    same(effect, expected.beta)
    values = []
    for i in range(len(data)):
        try:
            values.append(fit(np.delete(data, i, axis=0), require_t=False)[0])
        except (ValueError, np.linalg.LinAlgError):
            pass
    assert len(values) == expected.LOO_valid
    if values:
        same(min(values), expected.LOO_min); same(max(values), expected.LOO_max)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    out = root/'HNSCC_frozen_sender_association_results_20261016_v2'
    path = out/'private_independent_audit.json'; assert not path.exists()
    deadline = time.monotonic()+3600
    while not (out/'private_terminal.json').exists():
        receipt = json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/revise_HNSCC_sender_coefficient_bootstrap_v2/launch_receipt.json').read_text())
        process = Path('/proc')/str(receipt['pid'])/'stat'
        assert process.exists(), 'Original coefficient revision handle missing before terminal'
        state = process.read_text().rsplit(')', 1)[1].split()
        assert state[19] == receipt['proc_start_ticks'] and state[0] != 'Z', 'Revision no longer live before terminal'
        assert time.monotonic() < deadline, 'Revision observation window expired; preserve original job'
        time.sleep(10)
    terminal = json.loads((out/'private_terminal.json').read_text())
    result_path = out/'all400_private_primary_association_statistics.tsv'
    descriptor_path = out/'all800_private_control_and_interface_descriptors.tsv'
    assert hashlib.sha256(result_path.read_bytes()).hexdigest() == terminal['result_sha256']
    assert hashlib.sha256(descriptor_path.read_bytes()).hexdigest() == terminal['descriptor_sha256']
    prep = root/'HNSCC_frozen_sender_association_preflight_20261016_v1'
    config = json.loads((prep/'contract_before_any_association_effects.json').read_text())
    frames = pd.read_csv(prep/'private_all_planned_patient_designs.tsv', sep='\t')
    fields = ['exposure', 'response']+config['categorical_nuisance']+config['numeric_nuisance']
    stats = pd.read_csv(result_path, sep='\t')
    old_stats = pd.read_csv(root/'HNSCC_frozen_sender_association_results_20261016_v1/all400_private_primary_association_statistics.tsv', sep='\t')
    assert stats[['ligand','sender_group']].equals(old_stats[['ligand','sender_group']])
    assert np.allclose(stats[['beta','t','p','BH']], old_stats[['beta','t','p','BH']], equal_nan=True, atol=1e-13, rtol=1e-13)

    descriptors = pd.read_csv(descriptor_path, sep='\t').set_index(['ligand', 'sender_group', 'variant'])
    preflight = pd.read_csv(prep/'all400_private_design_identifiability.tsv', sep='\t').set_index(['ligand', 'sender_group'])
    assert len(stats) == 400 and len(descriptors) == 800
    # Independently rebuild each frozen original receiver gene contribution.
    g, source_rows, cells, var, eligible, measured = hnscc_inputs(root, project)
    selected = pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_selected_cells.tsv', sep='\t')
    selected = selected[selected.subtype_candidate.eq('Macro_C1QC')]
    pos = cells.index.get_indexer(selected.cell); assert (pos >= 0).all()
    raw = ad.io.sparse_dataset(g['layers/counts'])[source_rows[pos], :].tocsr(); raw.eliminate_zeros()
    total = np.asarray(raw.sum(axis=1)).ravel()
    patients = sorted(frames.patient_id.unique()); ids = cells.iloc[pos].patient_id.to_numpy()
    codes = np.array([patients.index(p) for p in ids]); sizes = np.bincount(codes, minlength=len(patients))
    w = pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv', sep='\t')
    w = w[w.native_column.ge(0)]
    contributions = {}
    for row in w.itertuples():
        col = int(row.native_column)
        v = raw[:, col].toarray().ravel()
        value = float(row.weight)*(np.log1p(10000*v/total)-float(row.mean))/float(row.sd)
        contributions[col] = np.bincount(codes, weights=value, minlength=len(patients))/sizes
    original = pd.read_csv(out/'private_original_receiver_gene_contributions.tsv', sep='\t', index_col=0).loc[patients]
    assert set(original.columns.astype(int)) == set(contributions)
    for col, value in contributions.items():
        assert np.allclose(original[str(col)], value, atol=1e-10, rtol=1e-10)
    removed = pd.read_csv(out/'private_all_interface_contribution_removals.tsv', sep='\t').set_index(['ligand', 'patient_id'])
    expression = pd.read_csv(root/'HNSCC_all_sender_ligand_measurements_20261016_v1/private_all_patient_group_gene_expression.tsv', sep='\t')
    columns = expression[['gene', 'native_column']].drop_duplicates().set_index('gene').native_column
    edges = pd.read_csv(root/'HNSCC_prior_receptor_qualification_20261016_v1/private_all_shortlist_receptor_edges.tsv', sep='\t')
    for ligand in sorted(stats.ligand.unique()):
        genes = {ligand}
        for receptor in edges[edges.ligand.eq(ligand) & edges.prior_reference_and_RNA_qualified].receptor:
            genes.update(receptor.split('_'))
        overlap = sorted({int(columns[g]) for g in genes if pd.notna(columns[g])} & set(contributions))
        values = sum((contributions[c] for c in overlap), np.zeros(len(patients)))
        assert np.allclose(removed.loc[[(ligand, p) for p in patients], 'removed_original_contribution'], values, atol=1e-10, rtol=1e-10)
        assert (removed.loc[[(ligand, p) for p in patients], 'removed_gene_count'] == len(overlap)).all()
    checked = 0
    with threadpool_limits(limits=1):
        for j, row in enumerate(stats.itertuples()):
            plan = preflight.loc[(row.ligand, row.sender_group)]
            sub = frames[frames.ligand.eq(row.ligand) & frames.sender_group.eq(row.sender_group) & frames.included_in_full_design].sort_values('patient_id')
            assert len(sub) == row.n_patients
            if plan.status != 'full_design_estimable_effects_pending':
                assert pd.isna(row.beta) and pd.isna(row.p) and pd.isna(row.BH)
                for variant in ['primary_full_control', 'interface_disjoint_full_control']:
                    d = descriptors.loc[(row.ligand, row.sender_group, variant)]
                    assert pd.isna(d.beta) and d.LOO_valid == 0 and pd.isna(d.LOO_min) and pd.isna(d.LOO_max)
                continue
            data = sub[fields].to_numpy(float)
            b, t, D, solver, variance, y, df = fit(data)
            same(b, row.beta); same(t, row.t); assert row.residual_df == df
            N = D[:, :-1]; nuisance_Q, _ = np.linalg.qr(N, mode='reduced')
            reduced = nuisance_Q@(nuisance_Q.T@y); residual = y-reduced
            full_Q, _ = np.linalg.qr(D, mode="reduced")
            strata = data[:, 2]+2*data[:, 3]
            groups = [np.flatnonzero(strata == k) for k in sorted(set(strata))]
            key = row.ligand+'|'+row.sender_group
            rng = np.random.default_rng(20261016+int(hashlib.sha256(key.encode()).hexdigest()[:8], 16))
            extreme = 0
            for _ in range(80):
                order = np.tile(np.arange(len(data))[:, None], (1, 250))
                for group in groups:
                    for col in range(250):
                        order[group, col] = rng.permutation(group)
                perm = reduced[:, None]+residual[order]
                coefficients = solver@perm; error = perm-full_Q@(full_Q.T@perm)
                se = np.sqrt(np.sum(error*error, axis=0)/df*variance)
                pt = np.divide(coefficients[-1], se, out=np.full(250, np.inf), where=se > 1e-12)
                extreme += np.count_nonzero(np.abs(pt) >= abs(t)-1e-12)
            assert (int(extreme)+1)/20001 == row.p or abs((int(extreme)+1)/20001-row.p) < 1e-15
            draws = []
            for _ in range(2000):
                ix = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
                try:
                    draws.append(fit(data[ix], require_t=False)[0])
                except (ValueError, np.linalg.LinAlgError):
                    pass
            assert len(draws) == row.valid_bootstrap
            lo, hi = np.quantile(draws, [.025, .975]) if len(draws) >= 1800 else (np.nan, np.nan)
            same(float(lo), row.lower); same(float(hi), row.upper)
            base = descriptors.loc[(row.ligand, row.sender_group, 'primary_full_control')]
            descriptor(data, base)
            disjoint = data.copy()
            disjoint[:, 1] -= removed.loc[[(row.ligand, p) for p in sub.patient_id], 'removed_original_contribution'].to_numpy()
            descriptor(disjoint, descriptors.loc[(row.ligand, row.sender_group, 'interface_disjoint_full_control')])
            r = ranks(data[:, :2]); a = r[:, 0]-r[:, 0].mean(); c = r[:, 1]-r[:, 1].mean()
            same(float(a@c/np.sqrt((a@a)*(c@c))), base.unadjusted_rho)
            checked += 1
            (out/'private_independent_audit_progress.json').write_text(json.dumps(dict(checked_primary_models=checked, observed_UTC=datetime.now(timezone.utc).isoformat())))
    assert checked == terminal['estimated']
    p = stats.p.fillna(1).to_numpy(); expected_q = np.ones(400)
    # Definition by all thresholds, rather than producer's sorted reverse cummin.
    for i, value in enumerate(p):
        eligible = np.flatnonzero(p >= value)
        expected_q[i] = min(1, min(400*p[j]/np.count_nonzero(p <= p[j]) for j in eligible))
    assert np.allclose(stats.BH[stats.p.notna()], expected_q[stats.p.notna()], atol=1e-12)
    report = dict(status='PASS_ALL400_V2_FAMILY_FULL_QR_PRIMARY_PERMUTATIONS_COEFFICIENT_BOOTSTRAP_NATIVE_INTERFACE_ALL_LOO',
        finished_UTC=datetime.now(timezone.utc).isoformat(), checked_primary_models=checked,
        result_sha256=terminal['result_sha256'], descriptor_sha256=terminal['descriptor_sha256'],
        direct_interaction_verified=False, no_new_selection_or_tests=True)
    path.write_text(json.dumps(report, indent=2)+'\n'); print(json.dumps(report))


if __name__ == '__main__':
    main()
