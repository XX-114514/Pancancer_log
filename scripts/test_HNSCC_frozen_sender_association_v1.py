"""One frozen full-control ligand/programme family; private native evidence."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def seed(key):
    return 20261016+int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)


def fitted(data):
    n = len(data)
    if n < 20 or n-14 < 10 or not np.isfinite(data).all():
        raise ValueError('sample_or_residual_df_or_missingness')
    c = np.column_stack((data[:, :2], data[:, 4:]))
    r = rankdata(c, axis=0, method='average')
    sd = r.std(axis=0)
    if (sd < 1e-12).any():
        raise ValueError('constant_required_continuous_field')
    r = (r-r.mean(axis=0))/sd
    x, y = r[:, 0], r[:, 1]
    N = np.column_stack((np.ones(n), data[:, 2:4], r[:, 2:]))
    D = np.column_stack((N, x))
    if np.linalg.matrix_rank(D) != 14 or np.linalg.cond(D) > 1e8:
        raise ValueError('full_design_rank_or_condition')
    Q, _ = np.linalg.qr(N, mode='reduced')
    rx = x-Q@(Q.T@x); ry = y-Q@(Q.T@y)
    xx = float(rx@rx)
    if xx < 1e-12:
        raise ValueError('zero_residual_exposure')
    beta = float(rx@ry/xx)
    error = ry-beta*rx; df = n-14
    se = np.sqrt(float(error@error)/df/xx)
    if se < 1e-12:
        raise ValueError('zero_residual_outcome_variance')
    return beta, beta/se, Q, rx, xx, y, df


def main_test(data, key):
    beta, t, Q, rx, xx, y, df = fitted(data)
    reduced = Q@(Q.T@y); residual = y-reduced
    strata = data[:, 2]+2*data[:, 3]
    groups = [np.flatnonzero(strata == value) for value in sorted(set(strata))]
    rng = np.random.default_rng(seed(key)); extreme = 0
    for _ in range(80):
        order = np.tile(np.arange(len(data))[:, None], (1, 250))
        for group in groups:
            for j in range(250):
                order[group, j] = rng.permutation(group)
        perm = reduced[:, None]+residual[order]
        yr = perm-Q@(Q.T@perm)
        b = rx@yr/xx
        ss = np.maximum(np.sum((yr-rx[:, None]*b)**2, axis=0), 0)
        se = np.sqrt(ss/df/xx)
        pt = np.divide(b, se, out=np.full(len(b), np.inf), where=se > 1e-12)
        extreme += int(np.count_nonzero(np.abs(pt) >= abs(t)-1e-12))
    draws = []
    for _ in range(2000):
        ix = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        try:
            draws.append(fitted(data[ix])[0])
        except (ValueError, np.linalg.LinAlgError):
            pass
    low, high = np.quantile(draws, [.025, .975]) if len(draws) >= 1800 else (np.nan, np.nan)
    return dict(beta=beta, t=float(t), p=(extreme+1)/20001, lower=float(low), upper=float(high),
        valid_bootstrap=len(draws), residual_df=df, status='estimated_exploratory_full_control')


def descriptors(data):
    try:
        effect = fitted(data)[0]
    except (ValueError, np.linalg.LinAlgError):
        return dict(beta=np.nan, LOO_valid=0, LOO_min=np.nan, LOO_max=np.nan, status='unestimable_full_control_descriptor')
    values = []
    for i in range(len(data)):
        try:
            values.append(fitted(np.delete(data, i, axis=0))[0])
        except (ValueError, np.linalg.LinAlgError):
            pass
    return dict(beta=effect, LOO_valid=len(values), LOO_min=min(values) if values else np.nan,
        LOO_max=max(values) if values else np.nan, status='estimated_descriptive_no_new_P')


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    prep = root/'HNSCC_frozen_sender_association_preflight_20261016_v1'
    audit = json.loads((prep/'private_independent_audit.json').read_text())
    assert audit['status'].startswith('PASS_ALL400_KEYS_NATIVE_COVARIATES')
    config = json.loads((prep/'contract_before_any_association_effects.json').read_text())
    path = prep/'private_all_planned_patient_designs.tsv'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == audit['design_sha256']
    frames = pd.read_csv(path, sep='\t')
    slots = pd.read_csv(prep/'all400_private_design_identifiability.tsv', sep='\t')
    out = root/'HNSCC_frozen_sender_association_results_20261016_v1'; assert not out.exists(); out.mkdir()
    contract = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        preflight_contract_sha256=hashlib.sha256((prep/'contract_before_any_association_effects.json').read_bytes()).hexdigest(),
        design_sha256=audit['design_sha256'], primary_family=400, nuisance=config['numeric_nuisance']+config['categorical_nuisance'],
        resampling='20000 Freedman-Lane residual permutations within treatment x HPV, batch250: group then column RNG order; then2000 patient bootstrap draws within the same sorted strata. Seed20261016+SHA256(ligand|sender) first8 hex.',
        interval='Rank transformations and all required controls refit in each bootstrap; >=1800 valid draws required for interval.',
        descriptor_scope='Every slot full-control and interface-disjoint full-control point and full LOO, plus unadjusted rank correlation; no new descriptor P/BH.',
        interface_disjoint='Subtract original frozen centred weighted contribution of ligand and all RNA-qualified receptor-chain physical columns from original programme, no refit. Same original receiver cells and full native RNA denominator. Full original programme remains primary.',
        interpretation='Prior and clinical lead selected receiver; exploratory, no mediation, direct binding or independent clinical prediction.',
        no_threshold_reselection_or_covariate_deletion=True)
    (out/'contract_before_effect_tests.json').write_text(json.dumps(contract, indent=2)+'\n')
    # Compute exact frozen gene contributions on the original receiver cell draws.
    source = root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1'
    selected = pd.read_csv(source/'private_selected_cells.tsv', sep='\t')
    selected = selected[selected.subtype_candidate.eq('Macro_C1QC')]
    g, rows, cells, var, eligible, measured = hnscc_inputs(root, project)
    pos = cells.index.get_indexer(selected.cell); assert (pos >= 0).all()
    raw = ad.io.sparse_dataset(g['layers/counts'])[rows[pos], :].tocsr().astype(float)
    raw.sum_duplicates(); raw.eliminate_zeros()
    total = np.asarray(raw.sum(axis=1)).ravel(); assert (total > 0).all()
    weight_path = root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv'
    prior_terminal = json.loads((root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1/private_terminal.json').read_text())
    assert hashlib.sha256(weight_path.read_bytes()).hexdigest() == prior_terminal['programme_weights_sha256']
    w = pd.read_csv(weight_path, sep='\t'); w = w[w.native_column.ge(0)]
    columns = w.native_column.to_numpy(int); assert measured[columns].all() and len(set(columns)) == len(columns)
    value = np.log1p(raw[:, columns].toarray()*10000/total[:, None])
    contribution = ((value-w['mean'].to_numpy())/w.sd.to_numpy())*w.weight.to_numpy()
    means = pd.DataFrame(contribution, index=cells.iloc[pos].patient_id.to_numpy(), columns=columns).groupby(level=0).mean()
    means.to_csv(out/'private_original_receiver_gene_contributions.tsv', sep='\t')
    expression = pd.read_csv(root/'HNSCC_all_sender_ligand_measurements_20261016_v1/private_all_patient_group_gene_expression.tsv', sep='\t')
    column_map = expression[['gene', 'native_column']].drop_duplicates().set_index('gene').native_column
    assert column_map.index.is_unique
    receptors = pd.read_csv(root/'HNSCC_prior_receptor_qualification_20261016_v1/private_all_shortlist_receptor_edges.tsv', sep='\t')
    response_rows = []; masks = {}
    for ligand in sorted(slots.ligand.unique()):
        genes = {ligand}
        for r in receptors[receptors.ligand.eq(ligand) & receptors.prior_reference_and_RNA_qualified].receptor:
            genes.update(r.split('_'))
        overlap = sorted({int(column_map[g]) for g in genes if pd.notna(column_map[g])} & set(columns))
        masks[ligand] = overlap
        for patient in sorted(frames.patient_id.unique()):
            response_rows.append(dict(ligand=ligand, patient_id=patient, removed_physical_columns='|'.join(map(str, overlap)),
                removed_gene_count=len(overlap), removed_original_contribution=float(means.loc[patient, overlap].sum()) if overlap else 0.))
    response = pd.DataFrame(response_rows).set_index(['ligand', 'patient_id'])
    response.to_csv(out/'private_all_interface_contribution_removals.tsv', sep='\t')
    results = []; desc = []
    fields = ['exposure', 'response']+config['categorical_nuisance']+config['numeric_nuisance']
    with threadpool_limits(limits=1):
        for j, slot in enumerate(slots.itertuples()):
            sub = frames[frames.ligand.eq(slot.ligand) & frames.sender_group.eq(slot.sender_group) & frames.included_in_full_design].sort_values('patient_id')
            data = sub[fields].to_numpy(float)
            row = dict(ligand=slot.ligand, sender_group=slot.sender_group, programme_sign=slot.programme_sign,
                n_patients=len(data), beta=np.nan, t=np.nan, p=np.nan, lower=np.nan, upper=np.nan,
                valid_bootstrap=0, residual_df=slot.residual_df, status=slot.status,
                preflight_reasons=slot.reasons, full_family_slots=400)
            if slot.status == 'full_design_estimable_effects_pending':
                row.update(main_test(data, slot.ligand+'|'+slot.sender_group))
                base = descriptors(data)
                altered = data.copy()
                altered[:, 1] -= response.loc[[(slot.ligand, p) for p in sub.patient_id], 'removed_original_contribution'].to_numpy()
                disjoint = descriptors(altered)
                r = rankdata(data[:, :2], axis=0)
                unadjusted = float(np.corrcoef(r.T)[0, 1])
            else:
                base = disjoint = dict(beta=np.nan, LOO_valid=0, LOO_min=np.nan, LOO_max=np.nan, status='primary_unestimable_or_ineligible')
                unadjusted = np.nan
            for variant, values in [('primary_full_control', base), ('interface_disjoint_full_control', disjoint)]:
                desc.append(dict(ligand=slot.ligand, sender_group=slot.sender_group, variant=variant,
                    n_patients=len(data), removed_interface_genes=len(masks[slot.ligand]),
                    unadjusted_rho=unadjusted, **values))
            results.append(row)
            (out/'private_progress.json').write_text(json.dumps(dict(completed_slots=j+1,
                completed_primary_models=sum(r['status'] == 'estimated_exploratory_full_control' for r in results),
                observed_UTC=datetime.now(timezone.utc).isoformat())))
    result = pd.DataFrame(results); p = result.p.fillna(1).to_numpy(); order = np.argsort(p); q = np.empty(len(p))
    q[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1])
    result['BH'] = np.where(result.p.notna(), q, np.nan)
    assert len(result) == 400 and len(desc) == 800
    result.to_csv(out/'all400_private_primary_association_statistics.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(desc).to_csv(out/'all800_private_control_and_interface_descriptors.tsv', sep='\t', index=False, na_rep='NA')
    terminal = dict(status='ALL400_FROZEN_FULL_CONTROL_SENDER_TESTS_AND_ALL800_DESCRIPTORS_COMPLETE_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(), planned_family=400, estimated=int(result.p.notna().sum()),
        nominal_P_below05=int(result.p.lt(.05).sum()), primary_BH_below05=int(result.BH.lt(.05).sum()),
        result_sha256=hashlib.sha256((out/'all400_private_primary_association_statistics.tsv').read_bytes()).hexdigest(),
        descriptor_sha256=hashlib.sha256((out/'all800_private_control_and_interface_descriptors.tsv').read_bytes()).hexdigest(),
        direct_interaction_verified=False, independent_audit_pending=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(terminal))
    print(result[result.p.notna()].sort_values('p').head(10).to_json(orient='records'))


if __name__ == '__main__':
    main()
