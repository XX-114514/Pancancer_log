"""Full selected native RNA, exact frequency, PCA and design readback before fit."""
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.io import mmread
from scipy.sparse import coo_matrix, csr_matrix, vstack

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget = resource(root, 8 * 1024**2)
    run = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    terminal = json.loads((run / 'private_terminal.json').read_text())
    assert sha(run / 'private_terminal.json') == '2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    for n, v in terminal['files'].items(): assert sha(run / n) == v['sha256']
    scope = root / 'immune_RNA_MCP_discovery_scope_20261018_v1'
    chosen = pd.read_csv(scope / 'private_metadata_only_case_source_partition.tsv', sep='\t', keep_default_na=False)
    chosen = chosen[~chosen.partition.eq('small_source_projection_only')].reset_index(drop=True)
    common = pd.read_csv(scope / 'private_marker_excluded_native_gene_map.tsv', sep='\t')
    universe = common.gene.to_numpy(str);ng = len(universe)
    md = pd.read_csv(run / 'private_selected_native_cell_depth.tsv', sep='\t', keep_default_na=False)
    freeze = pd.read_csv(run / 'private_frozen_original_cell_selection.tsv', sep='\t', keep_default_na=False)
    assert md[freeze.columns].equals(freeze)
    rows = md.global_row.to_numpy(np.int64);assert len(rows) == len(set(rows)) == terminal['selected_cells']
    out = root / 'immune_RNA_MCP_native_expression_input_audit_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out / 'plan_before_complete_native_input_acceptance.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        producer_terminal_sha256=sha(run / 'private_terminal.json'),
        scope='Every selected cell and common-gene raw opportunity; original obs material/major and deterministic cap; source-native gene physical mapping and all grouping aliases excluded; every written expression value and PC; exact integer1% eligibility/ranking; training-only PCA centering and quality/covariate design. No model or clinical outcome fit.'), indent=2) + '\n')
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    obs = {k: np.asarray(read_elem(z['obs'][k])) for k in ['dataset', 'library_id', 'major']}
    material = pd.Series(obs['dataset']).astype(str) + '::' + pd.Series(obs['library_id']).astype(str)
    material_counts = material.value_counts()
    for field in ['dataset', 'partition', 'fit_sample']:
        assert np.array_equal(md[field], chosen.iloc[md.material_index.to_numpy()][field])
    assert np.array_equal(material.iloc[rows], chosen.iloc[md.material_index.to_numpy()].sampleID)
    assert np.array_equal(md.annotated_scope_cells, material_counts.loc[material.iloc[rows]].to_numpy())
    assert (obs['major'][rows[md.view.eq('Myeloid').to_numpy()]] == 'Myeloid').all()
    assert (obs['major'][rows[~md.view.eq('Myeloid').to_numpy()]] == 'T').all()
    source = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    cache = np.load(source / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False)
    raw_marker = cache['raw_gene_counts'];marker_genes = cache['genes'].tolist();tr = cache['global_original_T_row_index']
    gene_contract = json.loads((source / 'private_source_gene_measurement_contract.json').read_text())
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed = json.loads(fixed_path.read_text());group_genes = set(marker_genes)
    core = np.count_nonzero(raw_marker[:, [marker_genes.index(g) for g in fixed['panels']['pan_T']]] >= 1, axis=1) >= 2
    reg = (np.count_nonzero(raw_marker[:, [marker_genes.index(g) for g in fixed['panels']['regulatory_associated']]] >= 1, axis=1) >= 2) & (raw_marker[:, marker_genes.index('FOXP3')] >= 1)
    eff = np.count_nonzero(raw_marker[:, [marker_genes.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] >= 1, axis=1) >= 2
    t_case = pd.Index(chosen.sampleID).get_indexer(material.iloc[tr])
    my = np.flatnonzero(obs['major'] == 'Myeloid');my_case = pd.Index(chosen.sampleID).get_indexer(material.iloc[my])
    ready = pd.read_csv(root / 'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1/private_all_original_libraries_multicell_view_counts.tsv', sep='\t')
    ready = ready[ready.threshold.eq(1)].set_index('sampleID')
    for view, candidate_rows, candidate_case in [
        ('R_FOXP3_capture', tr[core & reg & ~eff], t_case[core & reg & ~eff]),
        ('E_effector_capture', tr[core & eff & ~reg], t_case[core & eff & ~reg]),
        ('Myeloid', my, my_case)]:
        for case, selected in chosen.iterrows():
            candidates = candidate_rows[candidate_case == case]
            field = {'R_FOXP3_capture': 'FOXP3_supported_R_only', 'E_effector_capture': 'effector_only', 'Myeloid': 'original_myeloid_cells'}[view]
            assert len(candidates) == int(ready.loc[selected.sampleID, field])
            expected = sorted(candidates, key=lambda row: hashlib.sha256(('2026101802|' + view + '|' + str(int(row))).encode()).hexdigest())[:100]
            observed = md[md.view.eq(view) & md.material_index.eq(case)]
            assert sorted(observed.global_row.tolist()) == sorted(expected)
            assert (observed.full_view_cells == len(candidates)).all() and len(observed) >= 20
    source_order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    symbol = np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identities = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    aliases = np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    labels = [{s, c} | (set(a.split('|')) if st == 'stable_Ensembl_ID' else set()) for s, c, a, st in zip(symbol, canonical, aliases, identities)]
    for dataset in chosen.dataset.unique():
        columns = common[dataset].to_numpy(np.int64);si = source_order.index(dataset)
        assert len(set(columns)) == len(columns) and (columns >= 0).all() and measured[columns, si].all()
        for gene, col in zip(universe, columns):
            assert gene in labels[col] and not (labels[col] & group_genes)
            assert not canonical[col].startswith('AMBIGUOUS::') and identities[col] != 'ambiguous_source_feature'
    c = z['layers/counts'];nc, nf = map(int, c.attrs['shape'])
    totals = np.zeros(len(md), np.int64);detected = np.zeros(len(md), np.int32);blocks = []
    for start in range(0, nc, 16384):
        stop = min(nc, start + 16384);lo, hi = np.searchsorted(rows, [start, stop])
        if hi > lo:
            ptr = np.asarray(c['indptr'][start:stop + 1], np.int64)
            x = csr_matrix((np.asarray(c['data'][int(ptr[0]):int(ptr[-1])]), np.asarray(c['indices'][int(ptr[0]):int(ptr[-1])]), ptr - ptr[0]), shape=(stop - start, nf))[rows[lo:hi] - start]
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert np.isfinite(x.data).all() and (x.data > 0).all() and (x.data == np.floor(x.data)).all()
            totals[lo:hi] = np.add.reduceat(x.data.astype(np.int64), x.indptr[:-1]);detected[lo:hi] = np.diff(x.indptr)
            rr, cc, vv = [], [], []
            for dataset in md.iloc[lo:hi].dataset.unique():
                ix = np.flatnonzero(md.iloc[lo:hi].dataset.to_numpy() == dataset)
                assert measured[x[ix].indices, source_order.index(dataset)].all()
                data = x[ix][:, common[dataset].to_numpy(np.int64)].tocoo()
                rr.append(ix[data.row]);cc.append(data.col);vv.append(data.data)
            blocks.append(coo_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(hi - lo, ng)).tocsr())
        if start // 16384 % 16 == 0 or stop == nc:
            (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), original_rows_read=stop, selected_cells_replayed=int(hi))) + '\n')
    X = vstack(blocks, format='csr')
    assert X.shape == (len(md), ng) and np.array_equal(totals, md.fullRNA) and np.array_equal(detected, md.full_detected_genes)
    log = X.astype(float);log.data = np.log1p(log.data * np.repeat(10000 / totals, np.diff(log.indptr)))
    summary = []
    for view in ['R_FOXP3_capture', 'E_effector_capture', 'Myeloid']:
        take = md.view.eq(view).to_numpy();vmd = md[take].reset_index(drop=True)
        train = vmd.partition.eq('reviewed_case_discovery').to_numpy();xr = X[take];xl = log[take]
        assert vmd.loc[train, 'fit_sample'].nunique() == 38 and vmd.loc[~train, 'fit_sample'].nunique() == 28
        means = np.asarray(xl[train].sum(axis=0)).ravel() / train.sum()
        variances = np.asarray(xl[train].power(2).sum(axis=0)).ravel() / train.sum() - means**2
        positives = np.asarray((xr[train] > 0).sum(axis=0)).ravel().astype(np.int64)
        exact = (100 * positives >= int(train.sum())) & (variances > 1e-12)
        raw_quality = np.column_stack([np.log1p(totals[take]), np.log1p(detected[take]), vmd.full_view_cells.to_numpy() / vmd.annotated_scope_cells.to_numpy(), np.log1p(vmd.full_view_cells)])
        mu = raw_quality[train].mean(axis=0);sd = raw_quality[train].std(axis=0);q = (raw_quality - mu) / sd
        projection = np.load(run / 'private_inputs' / (view + '_frozen_projection.npz'), allow_pickle=False)
        assert np.allclose(projection['quality_mean'], mu) and np.allclose(projection['quality_sd'], sd)
        candidates = projection['candidate_genes'].tolist();cols = pd.Index(universe).get_indexer(candidates)
        assert (cols >= 0).all() and len(set(candidates)) == 5000 and exact[cols].all()
        excluded = exact.copy();excluded[cols] = False
        assert not (variances[excluded] > variances[cols].min() + 1e-10).any(), 'Exact integer1% top-variance candidate set differs'
        hvg = projection['hvg_genes'].tolist();hcols = pd.Index(universe).get_indexer(hvg)
        assert len(set(hvg)) == 2000 and set(hvg) <= set(candidates)
        excluded = exact.copy();excluded[hcols] = False
        assert not (variances[excluded] > variances[hcols].min() + 1e-10).any(), 'Exact integer1% HVG set differs'
        assert np.allclose(projection['training_candidate_log_mean'], means[cols], rtol=1e-10, atol=1e-12)
        assert np.allclose(projection['training_candidate_log_variance'], variances[cols], rtol=1e-9, atol=1e-11)
        assert np.allclose(projection['training_candidate_detection_fraction'], positives[cols] / train.sum(), rtol=1e-12, atol=1e-14)
        assert np.allclose(projection['pca_mean'], means[hcols], rtol=1e-10, atol=1e-12)
        comp = projection['pca_components'];assert comp.shape == (20, 2000) and np.allclose(comp @ comp.T, np.eye(20), atol=1e-8)
        pc = (xl[:, hcols].toarray() - projection['pca_mean']) @ comp.T
        for label, mask in [('discovery', train), ('TNBC_holdout', ~train)]:
            prefix = run / 'private_inputs' / (view + '_' + label)
            metadata = pd.read_csv(str(prefix) + '_metadata.tsv', sep='\t')
            stored = mmread(gzip.open(str(prefix) + '_tpm.mtx.gz', 'rb')).tocsr().T.tocsr()
            expected = xl[mask][:, cols];difference = stored - expected
            assert difference.nnz == 0 or np.max(np.abs(difference.data)) < 1e-12
            assert Path(str(prefix) + '_genes.tsv').read_text().splitlines() == candidates
            assert np.array_equal(metadata.cell, vmd.loc[mask, 'cell']) and np.array_equal(metadata['sample'], vmd.loc[mask, 'fit_sample'])
            for j, field in enumerate(['cellQ', 'geneQ', 'fraction', 'logN']): assert np.allclose(metadata[field], q[mask, j], rtol=1e-12, atol=1e-12)
            if label == 'discovery': assert np.array_equal(metadata.sourceCC, metadata.dataset.eq('GSE294109').to_numpy(float))
            else: assert metadata.sourceCC.isna().all()
            pcs = pd.read_csv(str(prefix) + '_PC.tsv', sep='\t', index_col=0)
            assert list(pcs.index) == metadata.cell.tolist() and np.allclose(pcs, pc[mask], rtol=1e-10, atol=1e-10)
        design = np.column_stack([np.ones(train.sum()), q[train, :2], vmd.loc[train, 'dataset'].eq('GSE294109').to_numpy(float), q[train, 2:]])
        assert np.linalg.matrix_rank(design) == 6 and np.linalg.cond(design) < 1e8
        summary.append(dict(view=view, cells=len(vmd), common_raw_gene_opportunities=int(xr.shape[0] * xr.shape[1]), exact_integer_frequency_and_ranking_pass=True,
            all_saved_gene_values_and_PCs_replayed=True, all_quality_controls_full_original_not_capped=True, candidate_genes=5000, HVG=2000, PCs=20, rank=6))
    (out / 'private_complete_input_acceptance.json').write_text(json.dumps(dict(status='ACCEPTED_ALL_SELECTED_NATIVE_RNA_FEATURE_PCA_AND_DESIGN_INPUTS',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(run / 'private_terminal.json'),
        selected_cells=len(md), views=summary, discovery_case_selection_not_new_identity_validation=True, no_MCP_or_clinical_fit=True), indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    t = dict(status='ACCEPTED_ALL_THREE_NATIVE_VIEWS_38_DISCOVERY_28_SOURCE_HOLDOUT_INPUTS',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(run / 'private_terminal.json'),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}, biological_or_new_patient_independence_claim=False)
    (out / 'private_terminal.json').write_text(json.dumps(t, indent=2) + '\n')
    print(json.dumps(dict(status=t['status'], completed_UTC=t['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
