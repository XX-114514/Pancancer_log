"""V4: official ULM unit-by-regulator orientation and complete raw replay.

Raw replay covers selected units, not all1.7M cells. Every planned source/slot
and all original fine-unit mappings are checked; no effects or revised q.
"""
import ast
import hashlib
import json
import os
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'existing_immune_fine_native_signatures_20261016_v4'
    receipt = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2/measure_existing_immune_fine_native_signatures_v4/launch_receipt.json').read_text())
    deadline = time.monotonic() + 21600
    while not (out / 'private_terminal.json').exists():
        p = Path('/proc') / str(receipt['pid']) / 'stat'
        assert p.exists(), 'Producer missing without terminal; inspect original job, never restart automatically'
        s = p.read_text().rsplit(')', 1)[1].split()
        assert s[19] == str(receipt['proc_start_ticks']) and s[0] != 'Z' and time.monotonic() < deadline
        time.sleep(15)
    terminal = json.loads((out / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    plan = json.loads((out / 'plan_before_native_expression_measurement.json').read_text())
    assert sha(out / 'private_prespecified64_raw_audit_units.tsv') == plan['metadata_only64_raw_audit_unit_sha256']
    run = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    assert sha(run / 'manifests/immune_authority_release_public.json') == plan['atlas_release_sha256']
    z = zarr.open_group(str(run / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    obs = pd.DataFrame({k: np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset', 'library_id', 'node_id', 'major']})
    units = pd.read_csv(out / 'private_all_original_fine_units_and_quality.tsv', sep='\t',
                        dtype={'dataset': str, 'library_id': str, 'node_id': str})
    mapping = np.load(root / 'existing_immune_fine_native_measurement_scope_20261016_v1/private_original_cell_to_fine_node_unit.npy')
    use = obs.major.isin(['T', 'Myeloid'])
    counts_by_key = Counter(tuple(x) for x in obs.loc[use, ['dataset', 'library_id', 'node_id', 'major']].itertuples(index=False, name=None))
    assert len(counts_by_key) == len(units) == 6957
    key_to_unit = {tuple(getattr(r, k) for k in ['dataset', 'library_id', 'node_id', 'major']): r.unit_index for r in units.itertuples()}
    assert set(key_to_unit) == set(counts_by_key)
    expected = np.fromiter((key_to_unit.get(tuple(row), -1) for row in obs.itertuples(index=False, name=None)), dtype=np.int32)
    assert np.array_equal(expected, mapping)
    assert np.array_equal(np.bincount(mapping[mapping >= 0], minlength=len(units)), units.n_cells)
    assert np.array_equal(units.later_n20_eligible, units.n_cells.ge(20))
    assert not units.exact_NSCLC_endpoint_mapping_validated.any()
    fraction = units.n_cells / units.groupby(['dataset', 'library_id', 'major']).n_cells.transform('sum')
    assert np.allclose(fraction, units.within_original_major_fraction, rtol=1e-12, atol=1e-12)
    selected = pd.read_csv(out / 'private_prespecified64_raw_audit_units.tsv', sep='\t')
    audit_ids = selected.unit_index.to_numpy(int)
    assert len(audit_ids) == len(set(audit_ids)) == 64
    assert set(selected.dataset) == set(units.dataset) and set(selected.node_id) == set(units.node_id)
    # Independent selector reconstruction from pre-expression metadata.
    priority = {r.unit_index: hashlib.sha256(('20261016|' + str(r.dataset) + '|' + str(r.library_id) + '|' + str(r.node_id)).encode()).hexdigest() for r in units.itertuples()}
    chosen = set()
    for field in ['dataset', 'node_id']:
        for _, f in units.groupby(field, sort=True):
            chosen.update(sorted(f.unit_index, key=lambda u: priority[u])[:2 if field == 'dataset' else 1])
    chosen.update(sorted(units.loc[units.n_cells.lt(20), 'unit_index'], key=lambda u: priority[u])[:2])
    for u in sorted(priority, key=priority.get):
        if len(chosen) == 64:
            break
        chosen.add(u)
    assert set(audit_ids) == chosen
    source_order = np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
    mask = np.asarray(read_elem(z['varm']['measured_in_source']), bool)
    symbol = np.asarray(read_elem(z['var']['gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var']['canonical_gene_identity'])).astype(str)
    aliases = np.asarray(read_elem(z['var']['source_gene_aliases'])).astype(str)
    status = np.asarray(read_elem(z['var']['feature_identity_status'])).astype(str)
    maps = {}
    for si, source in enumerate(source_order):
        names = defaultdict(set)
        for col in np.flatnonzero(mask[:, si]):
            if canonical[col].startswith('AMBIGUOUS::') or status[col] == 'ambiguous_source_feature':
                continue
            labels = [canonical[col], symbol[col]]
            if status[col] == 'stable_Ensembl_ID':
                labels += aliases[col].split('|')
            for gene in labels:
                if gene not in ['', 'nan', 'None', '__SOURCE_VALUE_MISSING__']:
                    names[gene].add(int(col))
        maps[source] = names
    panel_axis = pd.read_csv(out / 'private_all_source_panel_gene_axis.tsv.gz', sep='\t')
    panel = np.load(out / 'private_limited_native_panel.npz', allow_pickle=False)
    physical = panel['physical_columns']
    assert len(set(physical)) == len(physical)
    pc = {col: j for j, col in enumerate(physical)}
    for row in panel_axis.itertuples():
        found = maps[row.dataset][row.gene]
        column = next(iter(found)) if len(found) == 1 else -1
        assert row.physical_matches == len(found) and row.native_column == column
        assert row.source_native_unambiguous == (column >= 0)
        assert row.panel_column == pc.get(column, -1)
    shape = (64, len(symbol))
    raw_log = np.zeros(shape); raw_sum = np.zeros(shape)
    raw_positive = np.zeros(shape, dtype=np.int32); raw_ge2 = np.zeros(shape, dtype=np.int32)
    raw_quality = np.zeros((64, 4))
    local_unit = {u: j for j, u in enumerate(audit_ids)}
    raw = z['layers']['counts']
    for start in range(0, len(obs), 8192):
        stop = min(len(obs), start + 8192)
        selected_rows = np.isin(mapping[start:stop], audit_ids)
        if not selected_rows.any():
            continue
        ptr = np.asarray(raw['indptr'][start:stop + 1], dtype=np.int64)
        x = csr_matrix((np.asarray(raw['data'][int(ptr[0]):int(ptr[-1])]), np.asarray(raw['indices'][int(ptr[0]):int(ptr[-1])]), ptr - ptr[0]), shape=(stop - start, len(symbol)))[selected_rows].copy()
        x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
        assert (x.data > 0).all() and np.equal(x.data, np.floor(x.data)).all()
        ids = mapping[start:stop][selected_rows]
        total = np.asarray(x.sum(axis=1)).ravel().astype(float)
        detected = np.diff(x.indptr)
        logged = x.astype(float, copy=True)
        logged.data = np.log1p(logged.data * np.repeat(10000 / total, detected))
        for u in np.unique(ids):
            k = local_unit[u]; use_rows = ids == u
            # Per-unit direct row reductions, independent of producer group matrix.
            raw_sum[k] += np.asarray(x[use_rows].sum(axis=0)).ravel()
            raw_log[k] += np.asarray(logged[use_rows].sum(axis=0)).ravel()
            raw_positive[k] += np.asarray((x[use_rows] >= 1).sum(axis=0)).ravel().astype(np.int32)
            raw_ge2[k] += np.asarray((x[use_rows] >= 2).sum(axis=0)).ravel().astype(np.int32)
            raw_quality[k] += [int(use_rows.sum()), float(total[use_rows].sum()), float(np.log1p(total[use_rows]).sum()), float(np.log1p(detected[use_rows]).sum())]
        if start // 8192 % 32 == 0:
            (out / 'private_independent_audit_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), stage='independent64_raw_rows', rows_read=stop, selected_raw_cells=int(raw_quality[:, 0].sum()))) + '\n')
    subset = units.set_index('unit_index').loc[audit_ids]
    assert np.array_equal(raw_quality[:, 0], subset.n_cells)
    raw_log /= raw_quality[:, :1]
    for k, field in [(1, 'mean_full_RNA'), (2, 'mean_logRNA'), (3, 'mean_loggenes')]:
        assert np.allclose(raw_quality[:, k] / raw_quality[:, 0], subset[field], rtol=1e-10, atol=1e-10)
    for key, a in [('mean_logCP10K', raw_log[:, physical]), ('sum_raw_counts', raw_sum[:, physical]),
                   ('positive_cell_counts', raw_positive[:, physical]), ('ge2_cell_counts', raw_ge2[:, physical])]:
        assert np.allclose(panel[key][audit_ids], a, rtol=1e-10, atol=1e-10), key
    assert (panel['positive_cell_counts'] <= units.n_cells.to_numpy()[:, None]).all()
    assert (panel['ge2_cell_counts'] <= panel['positive_cell_counts']).all()
    for row in panel_axis.itertuples():
        if row.native_column >= 0:
            continue
        # Missing gene remains an axis reason; no synthetic negative column.
        assert row.panel_column == -1 and not row.source_native_unambiguous
    params_file = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    parents_file = root / 'NSCLC_automatic_ligand_response_review_20261016_v1/private_all_primary_BH_supported_rows.tsv'
    slots_file = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1/private_all_deduplicated_TF_slots.tsv'
    net_file = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/human_collectri_signed_complexes.tsv'
    axis_file = root / 'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv'
    for p, key in [(params_file, 'frozen_parameters_sha256'), (parents_file, 'frozen_parents_sha256'), (slots_file, 'frozen_TF_scope_sha256'), (net_file, 'signed_network_sha256'), (axis_file, 'discovery_gene_axis_sha256')]:
        assert sha(p) == plan[key]
    params = pd.read_csv(params_file, sep='\t')
    parents = pd.read_csv(parents_file, sep='\t')
    slots = pd.read_csv(slots_file, sep='\t').set_index(['record', 'TF'])
    net = pd.read_csv(net_file, sep='\t')
    targets = {s: f for s, f in net.groupby('source')}
    discovery = pd.read_csv(axis_file, sep='\t').gene.tolist()
    pm = pd.read_csv(out / 'private_all_source_parent_programme_scope.tsv.gz', sep='\t')
    tm = pd.read_csv(out / 'private_all_source_parent_TF_scope.tsv.gz', sep='\t')
    assert len(pm) == 9 * len(parents) * 3 and not pm.duplicated(['dataset', 'record', 'role']).any()
    assert len(tm) == 9 * len(slots) and not tm.duplicated(['dataset', 'record', 'TF']).any()
    assert set(zip(tm.dataset, tm.record, tm.TF)) == {(s, p, t) for s in source_order for p, t in slots.index}
    scores = np.load(out / 'private_frozen_programme_and_TF_scores.npz', allow_pickle=False)
    ps = scores['programme']; ts = scores['TF']
    official = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt'
    original_TF_plan = json.loads((root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1/plan_before_any_TF_activity_or_association.json').read_text())
    assert sha(official) == original_TF_plan['official_ULM_source_sha256']
    tree = ast.parse(official.read_text())
    ns = {'np': np}
    exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ['mat_cov', 'mat_cor', 't_val']], type_ignores=[]), str(official), 'exec'), ns)
    test_W = np.column_stack([np.arange(61) % 3 - 1, np.arange(61) % 5 - 2]).astype(float)
    test_Y = np.random.default_rng(20261016).normal(size=(61, 5))
    test_C = ns['mat_cor'](test_W, test_Y)
    assert test_C.shape == (5, 2)
    for i in range(5):
        for j in range(2):
            assert np.isclose(test_C[i, j], np.corrcoef(test_W[:, j], test_Y[:, i])[0, 1])
    max_programme = 0.; max_TF = 0.; compared_programmes = 0; compared_TFs = 0
    for source in source_order:
        names = maps[source]
        lookup = lambda gene: next(iter(names[gene])) if len(names[gene]) == 1 else -1
        all_units = units.index[units.dataset.eq(source)].to_numpy()
        audited_local = np.flatnonzero(subset.dataset.to_numpy() == source)
        audited_global = audit_ids[audited_local]
        for parent in parents.itertuples():
            pp = params[params.context.eq(parent.context) & params.program.eq(parent.program)]
            sender_type = 'typeB' if parent.receiver_type == 'typeA' else 'typeA'
            sender_genes = set(pp.loc[pp.type.eq(sender_type), 'gene'])
            interfaces = {parent.ligand} | {g for receptor in parent.qualified_receptors.split('|') for g in receptor.split('_')}
            forbidden = set(pp.gene) | interfaces
            for row in pm[pm.dataset.eq(source) & pm.record.eq(parent.record)].itertuples():
                p = pp[pp.type.eq(row.type)]
                total_l1 = float(p.weight.abs().sum())
                if row.role == 'receiver_masked':
                    p = p[~p.gene.isin(sender_genes | interfaces)]
                cols = [lookup(g) for g in p.gene]
                cc = Counter(cols)
                kept = np.array([c >= 0 and cc[c] == 1 for c in cols], dtype=bool)
                fraction = float(p.weight.abs().to_numpy()[kept].sum() / p.weight.abs().sum())
                retained = float(p.weight.abs().sum() / total_l1)
                valid = fraction >= .95 and (row.role != 'receiver_masked' or (retained >= .8 and kept.sum() >= 20))
                assert np.isclose(row.source_measured_L1_fraction, fraction) and np.isclose(row.original_L1_retained_after_mask, retained)
                assert row.measured_genes == kept.sum() and row.status.startswith('MEASURED_') == valid
                if valid:
                    assert np.isfinite(ps[all_units, row.score_column]).all()
                    p = p.loc[kept]; c = np.asarray(cols)[kept]
                    values = np.sum(((raw_log[np.ix_(audited_local, c)] - p.gene_mean.to_numpy()) / p.gene_sd.to_numpy()) * p.weight.to_numpy(), axis=1)
                    delta = np.max(np.abs(values - ps[audited_global, row.score_column]))
                    assert np.allclose(values, ps[audited_global, row.score_column], rtol=1e-9, atol=1e-8)
                    max_programme = max(max_programme, float(delta)); compared_programmes += len(values)
                else:
                    assert np.isnan(ps[all_units, row.score_column]).all()
            bg_genes = [g for g in discovery if g not in forbidden]
            bg_cols = [lookup(g) for g in bg_genes]; cc = Counter(bg_cols)
            remove = {c for gene in forbidden for c in names[gene]}
            bg = sorted({c for c in bg_cols if c >= 0 and cc[c] == 1 and c not in remove})
            bg_index = {c: i for i, c in enumerate(bg)}
            active_rows = []; weights = []
            for row in tm[tm.dataset.eq(source) & tm.record.eq(parent.record)].itertuples():
                f = targets[row.TF]; f = f[f.target.isin(bg_genes)]
                cols = [lookup(g) for g in f.target]; cc = Counter(cols)
                keep = np.array([c in bg_index and cc[c] == 1 for c in cols], dtype=bool)
                original_valid = slots.loc[(row.record, row.TF), 'status'] == 'READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED'
                valid = original_valid and keep.sum() >= 20 and len(bg) >= 500
                assert row.original_disjoint_targets == len(f) and row.source_disjoint_targets == keep.sum()
                assert row.source_disjoint_background == len(bg) and row.original_background == len(bg_genes)
                assert row.status.startswith('MEASURED_') == valid
                if valid:
                    assert np.isfinite(ts[all_units, row.score_column]).all()
                    w = np.zeros(len(bg))
                    for c, sign in zip(np.asarray(cols)[keep], f.weight.to_numpy()[keep]):
                        w[bg_index[c]] = sign
                    weights.append(w); active_rows.append(row)
                else:
                    assert np.isnan(ts[all_units, row.score_column]).all()
            if active_rows:
                values = ns['t_val'](ns['mat_cor'](np.column_stack(weights), raw_log[np.ix_(audited_local, bg)].T), len(bg) - 2)
                assert values.shape == (len(audited_local), len(active_rows))
                for j, row in enumerate(active_rows):
                    old = ts[audited_global, row.score_column]
                    assert np.allclose(values[:, j], old, rtol=1e-8, atol=1e-8), (source, row.record, row.TF)
                    max_TF = max(max_TF, float(np.max(np.abs(values[:, j] - old)))); compared_TFs += len(old)
    proof = dict(status='PASS_ALL_FINE_UNIT_KEYS_ALL_SOURCE_SLOT_COVERAGE_AND_PRESPECIFIED64_RAW_UNIT_PROGRAMMES_OFFICIAL_TF_AND_PANEL_REPLAY',
                 finished_UTC=datetime.now(timezone.utc).isoformat(), audit_source_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'),
                 all_units=len(units), selected_raw_units=64, selected_raw_cells=int(raw_quality[:, 0].sum()),
                 programme_values_compared=compared_programmes, TF_values_compared=compared_TFs,
                 maximum_programme_difference=max_programme, maximum_official_TF_difference=max_TF,
                 full_raw_all_cell_independent_replay=False, official_ULM_source_sha256=sha(official), official_orientation_synthetic_Pearson_contract_pass=True, no_effects_P_BH_or_endpoint_promotion=True)
    dest = out / 'private_independent_audit.json'
    assert not dest.exists()
    dest.write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
