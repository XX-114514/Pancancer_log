"""All regulator RNA, target redundancy and sender-shared RNA descriptions.

No regulator expression filter, new P/BH or independent-validation claim.
"""
import ast
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import qr
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def rho(x, y, nuisance):
    xx = rankdata(x).astype(float)
    yy = rankdata(y).astype(float)
    if np.std(xx) <= 1e-12 or np.std(yy) <= 1e-12:
        return np.nan
    full = np.column_stack([nuisance, (xx - np.mean(xx)) / np.std(xx)])
    if np.linalg.matrix_rank(nuisance) != nuisance.shape[1] or np.linalg.matrix_rank(full) != full.shape[1] or np.linalg.cond(full) >= 1e8:
        return np.nan
    q = qr(nuisance, mode='economic')[0]
    xx -= q @ (q.T @ xx)
    yy -= q @ (q.T @ yy)
    den = np.linalg.norm(xx) * np.linalg.norm(yy)
    return float(xx @ yy / den) if den > 1e-12 else np.nan


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    desc = root / 'NSCLC_gene_disjoint_TF_response_descriptions_20261016_v1'
    da = json.loads((desc / 'private_independent_audit.json').read_text())
    assert da['status'].startswith('PASS_ALL_FIXED_TF_POINTS_')
    dt = json.loads((desc / 'private_terminal.json').read_text())
    dp = desc / 'private_all2568_TF_functional_descriptions.tsv'
    assert sha(dp) == dt['files'][dp.name]
    table = pd.read_csv(dp, sep='\t')
    active = table[table.description_status.eq('DESCRIPTIVE_FULL_CONTROL_POINT_ONLY')]
    inputs = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    mp = json.loads((inputs / 'plan_before_any_TF_measurement.json').read_text())
    pp = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1/plan_before_any_TF_activity_or_association.json'
    sp = json.loads(pp.read_text())
    params_path = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    net_path = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/human_collectri_signed_complexes.tsv'
    official = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt'
    assert sha(params_path) == mp['private_frozen_parameters_sha256'] and sha(net_path) == mp['signed_TF_network_sha256'] and sha(official) == sp['official_ULM_source_sha256']
    params = pd.read_csv(params_path, sep='\t')
    net = pd.read_csv(net_path, sep='\t')
    nodes = ast.parse(official.read_text())
    functions = [n for n in nodes.body if isinstance(n, ast.FunctionDef) and n.name in {'mat_cov', 'mat_cor', 't_val'}]
    assert len(functions) == 3
    ns = {'np': np}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(official), 'exec'), ns)
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    genes = pd.Index(pd.read_csv(native / 'private_native_gene_axis.tsv', sep='\t').gene)
    units = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t')
    means = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    positive = np.load(native / 'positive_cell_counts.npy', mmap_mode='r')
    parents = pd.read_csv(root / 'NSCLC_automatic_ligand_response_results_20261016_v1/private_all_unique_ligand_response_results.tsv', sep='\t').set_index('record')
    hashes = pd.read_csv(inputs / 'private_ready_design_SHA256.tsv', sep='\t').set_index('TF_record').sha256
    out = root / 'NSCLC_TF_specificity_shared_response_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), functional_audit_sha256=sha(desc / 'private_independent_audit.json'), functional_table_sha256=sha(dp),
        all_fixed_slots=len(table), all_measured_regulator_slots=len(active),
        regulator_RNA='All literal native regulator genes: pooled >=1 cell RNA and case coverage in both original endpoints. Nonliteral complexes unresolved, not absent. No expression-based exclusion or relabel.',
        redundancy='Every within-parent regulator pair: signed native target cosine and receiver residual-rank score correlation after full original Z. First PC variance fraction of standardized residual signatures; orient largest loading positive, no ligand-directed sign selection.',
        shared_response='For EVERY measured regulator compute its SAME disjoint signed target score in original sender cells; add sender same-regulator score to original full Z, retaining all controls. This is competing/shared-programme RNA description, not proof of confounding versus mediation.',
        sender_scoring='Pinned official decoupler1.9.2 score algebra, all native genes outside both original programmes and original interface. All original cases and full-RNA normalization.',
        limits='Regulator target signatures are not regulator protein or DNA-binding activity. Direct-TF versus cofactor membership not established here; mRNA nondetection does not prove protein absence. No new P/BH or independent samples.')
    (out / 'plan_before_specificity_diagnostics.json').write_text(json.dumps(plan, indent=2) + '\n')
    records, pairs, pcs, case_scores = [], [], [], []
    for parent_id, f in active.groupby('record', sort=True):
        parent = parents.loc[parent_id]
        loaded = []
        for r in f.itertuples():
            p = inputs / 'private_frozen_designs' / (r.TF_record + '.npz')
            assert sha(p) == hashes.loc[r.TF_record]
            loaded.append(np.load(p, allow_pickle=False))
        z = loaded[0]
        tf_names = pd.Index(f.TF)
        pp = params[params.context.eq(parent.context) & params.program.eq(parent.program)]
        forbidden = set(pp.gene) | {parent.ligand} | {g for receptor in parent.qualified_receptors.split('|') for g in receptor.split('_')}
        bg = pd.Index([g for g in genes if g not in forbidden])
        weights = np.zeros((len(bg), len(f)))
        for j, tf in enumerate(tf_names):
            group = net[net.source.eq(tf) & net.target.isin(bg)]
            assert len(group) >= 20
            weights[bg.get_indexer(group.target), j] = group.weight.to_numpy(float)
        sender = np.asarray(means[np.ix_(z['sender_units'].astype(int), genes.get_indexer(bg))], float)
        sender_scores = ns['t_val'](ns['mat_cor'](weights, sender.T), len(bg) - 2)
        assert np.isfinite(sender_scores).all()
        receiver_scores = np.column_stack([zz['Y'] for zz in loaded])
        Z = z['Z']
        q = qr(Z, mode='economic')[0]
        ranks = rankdata(receiver_scores, axis=0).astype(float)
        residual = ranks - q @ (q.T @ ranks)
        standard = residual / np.linalg.norm(residual, axis=0)
        score_correlation = standard.T @ standard
        w = csr_matrix(weights)
        norm = np.linalg.norm(weights, axis=0)
        target_cosine = (w.T @ w).toarray() / (norm[:, None] * norm[None, :])
        u, singular, vt = np.linalg.svd(standard, full_matrices=False)
        orientation = 1 if vt[0, np.argmax(np.abs(vt[0]))] >= 0 else -1
        pc = u[:, 0] * singular[0] * orientation
        xr = rankdata(z['X']).astype(float)
        xr -= q @ (q.T @ xr)
        pc_rho = float(xr @ pc / (np.linalg.norm(xr) * np.linalg.norm(pc)))
        pcs.append(dict(record=parent_id, regulators=len(f), cases=len(sender), first_PC_variance_fraction=float(singular[0] ** 2 / np.sum(singular ** 2)), residual_ligand_first_PC_correlation=pc_rho, PC_sign_chosen_without_ligand=True))
        for j, row in enumerate(f.itertuples()):
            same_score = sender_scores[:, j]
            sd = np.std(same_score)
            zz = np.column_stack([Z, (same_score - np.mean(same_score)) / sd]) if sd > 1e-12 else None
            adjusted = rho(z['X'], receiver_scores[:, j], zz) if zz is not None else np.nan
            observed_gene = genes.get_indexer([row.TF])[0]
            record = dict(TF_record=row.TF_record, record=parent_id, TF=row.TF, ligand=parent.ligand, literal_regulator_gene_measured=observed_gene >= 0, regulator_RNA_status='LITERAL_NATIVE_GENE' if observed_gene >= 0 else 'NONLITERAL_OR_COMPLEX_UNRESOLVED_NOT_ABSENT', sender_same_regulator_adjusted_rho=adjusted, sender_regulator_to_ligand_rho=rho(z['X'], same_score, Z), sender_to_receiver_same_regulator_rho=rho(same_score, receiver_scores[:, j], Z), original_partial_rank_rho=row.partial_rank_rho, original_sender_programme_adjusted_rho=row.sender_programme_adjusted_rho,
                maximum_absolute_residual_score_correlation=float(np.max(np.abs(np.delete(score_correlation[j], j)))) if len(f) > 1 else None, maximum_absolute_target_cosine=float(np.max(np.abs(np.delete(target_cosine[j], j)))) if len(f) > 1 else None)
            if observed_gene >= 0:
                for role in ['sender', 'receiver']:
                    indices = z[role + '_units'].astype(int)
                    counts = np.asarray(positive[indices, observed_gene], dtype=np.int64)
                    total = int(units.n_cells.to_numpy(int)[indices].sum())
                    record[role + '_regulator_RNA_positive_cells'] = int(counts.sum())
                    record[role + '_regulator_RNA_positive_fraction'] = int(counts.sum()) / total
                    record[role + '_regulator_RNA_positive_cases'] = int(np.count_nonzero(counts))
            records.append(record)
            case_scores.extend(dict(TF_record=row.TF_record, case_index=k, sender_same_regulator_score=float(same_score[k])) for k in range(len(sender)))
            for k in range(j + 1, len(f)):
                pairs.append(dict(record=parent_id, TF_A=str(tf_names[j]), TF_B=str(tf_names[k]), signed_target_cosine=float(target_cosine[j, k]), receiver_residual_rank_correlation=float(score_correlation[j, k])))
        for item in loaded:
            item.close()
    pd.DataFrame(records).to_csv(out / 'private_all_regulator_specificity_and_shared_response.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(pairs).to_csv(out / 'private_all_regulator_pair_redundancy.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(pcs).to_csv(out / 'private_all_parent_residual_PC_descriptions.tsv', sep='\t', index=False)
    pd.DataFrame(case_scores).to_csv(out / 'private_all_sender_same_regulator_case_scores.tsv.gz', sep='\t', index=False, compression='gzip')
    result = dict(status='COMPLETE_ALL_REGULATOR_RNA_REDUNDANCY_AND_SENDER_SHARED_RESPONSE_DESCRIPTIONS_AUDIT_PENDING', finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), all_fixed_scope=len(table), described_regulators=len(records), all_regulator_pairs=len(pairs), all_parent_PCs=len(pcs), all_sender_score_case_records=len(case_scores), files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_new_P_BH_or_expression_based_selection=True, regulator_protein_activity_or_independent_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
