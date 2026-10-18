"""Outcome-free, source-native expression inputs for three frozen RNA views."""
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.io import mmwrite
from scipy.sparse import coo_matrix, csr_matrix, vstack
from sklearn.decomposition import PCA

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget = resource(root, 256 * 1024**2)
    scope = root / 'immune_RNA_MCP_discovery_scope_20261018_v1'
    terminal = json.loads((scope / 'private_terminal.json').read_text())
    assert sha(scope / 'private_terminal.json') == '0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    for n, v in terminal['files'].items(): assert sha(scope / n) == v['sha256']
    all_chosen = pd.read_csv(scope / 'private_metadata_only_case_source_partition.tsv', sep='\t', keep_default_na=False)
    chosen = all_chosen[~all_chosen.partition.eq('small_source_projection_only')].copy().reset_index(drop=True)
    gene_map = pd.read_csv(scope / 'private_marker_excluded_native_gene_map.tsv', sep='\t')
    genes = gene_map.gene.to_numpy(str);ng = len(genes)
    source = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    raw_audit = root / 'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1/private_terminal.json'
    audit = json.loads(raw_audit.read_text());assert audit['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_')
    assert sha(source / 'private_terminal.json') == audit['producer_terminal_sha256']
    for n, v in json.loads((source / 'private_terminal.json').read_text())['files'].items(): assert sha(source / n) == v['sha256']
    cache = np.load(source / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False)
    raw = cache['raw_gene_counts'];marker_names = cache['genes'].tolist();tu = cache['original_T_unit_index']
    T = pd.read_csv(source / 'private_original_T_unit_scope.tsv', sep='\t', keep_default_na=False)
    ready = root / 'immune_CM_RNA_supported_multicell_views_preflight_20261018_v1/private_all_original_libraries_multicell_view_counts.tsv'
    full = pd.read_csv(ready, sep='\t', keep_default_na=False);full = full[full.threshold.eq(1)].set_index('sampleID')
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed = json.loads(fixed_path.read_text())
    assert not (set(genes) & set(marker_names))
    core = (raw[:, [marker_names.index(g) for g in fixed['panels']['pan_T']]] > 0).sum(axis=1) >= 2
    reg = ((raw[:, [marker_names.index(g) for g in fixed['panels']['regulatory_associated']]] > 0).sum(axis=1) >= 2) & (raw[:, marker_names.index('FOXP3')] > 0)
    eff = (raw[:, [marker_names.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] > 0).sum(axis=1) >= 2
    sample_lookup = pd.Index(chosen.sampleID)
    t_material = sample_lookup.get_indexer(T.dataset.astype(str) + '::' + T.library_id.astype(str))[tu]
    old_scope = root / 'existing_immune_fine_native_measurement_scope_20261016_v1'
    old_terminal = json.loads((old_scope / 'private_terminal.json').read_text())
    for n, digest in old_terminal['files'].items(): assert sha(old_scope / n) == digest
    old_units = pd.read_csv(old_scope / 'private_all_frozen_T_myeloid_node_units.tsv', sep='\t', keep_default_na=False)
    old_mapping = np.load(old_scope / 'private_original_cell_to_fine_node_unit.npy', mmap_mode='r')
    my_material = sample_lookup.get_indexer(old_units.dataset.astype(str) + '::' + old_units.library_id.astype(str))
    my_material[~old_units.major.eq('Myeloid').to_numpy()] = -1
    my_global = np.flatnonzero(old_mapping >= 0)
    my_case = my_material[old_mapping[my_global]]
    all_rows = []
    for view, positions, cases in [
        ('R_FOXP3_capture', cache['global_original_T_row_index'][core & reg & ~eff & (t_material >= 0)], t_material[core & reg & ~eff & (t_material >= 0)]),
        ('E_effector_capture', cache['global_original_T_row_index'][core & eff & ~reg & (t_material >= 0)], t_material[core & eff & ~reg & (t_material >= 0)]),
        ('Myeloid', my_global[my_case >= 0], my_case[my_case >= 0])]:
        for case in range(len(chosen)):
            candidates = positions[cases == case]
            field = {'R_FOXP3_capture': 'FOXP3_supported_R_only', 'E_effector_capture': 'effector_only', 'Myeloid': 'original_myeloid_cells'}[view]
            expected = int(full.loc[chosen.iloc[case].sampleID, field])
            assert len(candidates) == expected and expected >= 20
            hashed = sorted(candidates, key=lambda row: hashlib.sha256(('2026101802|' + view + '|' + str(int(row))).encode()).hexdigest())
            for row in hashed[:100]:
                all_rows.append(dict(global_row=int(row), view=view, material_index=case, full_view_cells=expected))
    md = pd.DataFrame(all_rows).sort_values('global_row').reset_index(drop=True)
    assert md.global_row.is_unique
    for field in ['fit_sample', 'dataset', 'partition', 'annotated_scope_cells']:
        md[field] = chosen.iloc[md.material_index.to_numpy()][field].to_numpy()
    assert (md.full_view_cells <= md.annotated_scope_cells).all()
    md['cell'] = md.global_row.map(lambda i: 'RNA_MCP_CELL_%09d' % i)
    out = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    inputs = out / 'private_inputs';inputs.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        scope_terminal_sha256=sha(scope / 'private_terminal.json'), original_T_RNA_audit_sha256=sha(raw_audit),
        raw_original_cell_scope_sha256=sha(old_scope / 'private_terminal.json'),
        scope='All38 reviewed discovery cases and28 TNBC source-heldout sample cases. Small RMC source3 retained in scope but optional gene-coverage projection deferred, not assigned zero missing expression.',
        selection='Mutually exclusive original RNA capture views, exactly as accepted counts. Metadata-only hash caps100 original cells per material/view, independent of RNA magnitude or outcome. No new annotation.',
        expression='All common32207 source-measured marker-excluded opportunities, then discovery-only >=1% raw positive and positive logCP10K variance; retain up to5000 variance-ranked candidate genes per view, first2000 for PCA20. Full original native RNA denominator.',
        covariates=['cellQ', 'geneQ', 'sourceCC', 'fraction', 'logN'],
        controls='cellQ=log fullRNA; geneQ=log detected genes; sourceCC discovery ccRCC indicator; fraction=full original view cells/all original material cells; logN=log1p full view cells. Continuous standardization fitted on discovery only. No capped-count fraction. Heldout sourceCC is NA and is never passed to discovery fitting.',
        unknown_phase='Discovery includes pre/post/unknown, latest stage overlay kept. No pretreatment-only clinical or treatment-specific claim.',
        output_cap_bytes=256 * 1024**2, input_acceptance_before_official_fit=True,
        no_clinical_outcomes_loaded=True)
    (out / 'plan_before_selected_native_expression_read.json').write_text(json.dumps(plan, indent=2) + '\n')
    md.to_csv(out / 'private_frozen_original_cell_selection.tsv', sep='\t', index=False)
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    counts = z['layers/counts'];nc, nf = map(int, counts.attrs['shape'])
    rows = md.global_row.to_numpy();total = np.zeros(len(md), np.int64);detected = np.zeros(len(md), np.int32)
    block_list = []
    for start in range(0, nc, 8192):
        stop = min(nc, start + 8192);a, b = np.searchsorted(rows, [start, stop])
        if b > a:
            ptr = np.asarray(counts['indptr'][start:stop + 1], np.int64)
            x = csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]), np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]), ptr - ptr[0]), shape=(stop - start, nf))[rows[a:b] - start].copy()
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert (x.data > 0).all() and np.isfinite(x.data).all() and (x.data == np.floor(x.data)).all()
            total[a:b] = np.asarray(x.sum(axis=1)).ravel().astype(np.int64);detected[a:b] = np.diff(x.indptr)
            rr, cc, vv = [], [], []
            for source_name in md.iloc[a:b].dataset.unique():
                ix = np.flatnonzero(md.iloc[a:b].dataset.to_numpy() == source_name)
                col = gene_map[source_name].to_numpy(np.int64)
                assert (col >= 0).all() and len(set(col)) == len(col)
                v = x[ix][:, col].tocoo();rr.append(ix[v.row]);cc.append(v.col);vv.append(v.data)
            block_list.append(coo_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(b - a, ng)).tocsr())
        if start // 8192 % 16 == 0 or stop == nc:
            (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), original_rows_read=stop, selected_cells_read=int(b), selected_cells=len(md))) + '\n')
    X = vstack(block_list, format='csr');assert X.shape == (len(md), ng) and (total > 0).all()
    t = ~md.view.eq('Myeloid').to_numpy();where = np.searchsorted(cache['global_original_T_row_index'], rows[t])
    assert np.array_equal(cache['global_original_T_row_index'][where], rows[t])
    assert np.array_equal(cache['full_native_total_counts'][where], total[t]) and np.array_equal(cache['full_native_detected_genes'][where], detected[t])
    logged = X.astype(np.float64);logged.data = np.log1p(logged.data * np.repeat(10000 / total, np.diff(logged.indptr)))
    md['fullRNA'] = total;md['full_detected_genes'] = detected
    summaries = []
    for view in ['R_FOXP3_capture', 'E_effector_capture', 'Myeloid']:
        take = md.view.eq(view).to_numpy();vmd = md.loc[take].copy().reset_index(drop=True)
        train = vmd.partition.eq('reviewed_case_discovery').to_numpy()
        xp = X[take];lp = logged[take]
        detection_fraction = np.asarray((xp[train] > 0).mean(axis=0)).ravel()
        mean = np.asarray(lp[train].mean(axis=0)).ravel();variance = np.asarray(lp[train].power(2).mean(axis=0)).ravel() - mean**2
        eligible = np.flatnonzero((detection_fraction >= .01) & (variance > 1e-12))
        keep = eligible[np.lexsort((genes[eligible], -variance[eligible]))[:5000]]
        assert len(keep) >= 2000
        hvg = keep[:2000]
        pca = PCA(n_components=20, svd_solver='randomized', random_state=2026101802)
        pca.fit(lp[train][:, hvg].toarray());pc = pca.transform(lp[:, hvg].toarray())
        quality = np.column_stack([np.log1p(vmd.fullRNA), np.log1p(vmd.full_detected_genes),
            vmd.full_view_cells.to_numpy() / vmd.annotated_scope_cells.to_numpy(), np.log1p(vmd.full_view_cells)])
        mu = quality[train].mean(axis=0);sd = quality[train].std(axis=0);assert (sd > 0).all()
        scaled = (quality - mu) / sd
        for j, key in enumerate(['cellQ', 'geneQ', 'fraction', 'logN']):vmd[key] = scaled[:, j]
        vmd['sourceCC'] = np.where(train, vmd.dataset.eq('GSE294109').to_numpy(float), np.nan)
        design = np.column_stack([np.ones(train.sum()), scaled[train, :2], vmd.loc[train, 'sourceCC'], scaled[train, 2:]])
        assert np.linalg.matrix_rank(design) == design.shape[1] and np.linalg.cond(design) < 1e8
        vmd['sample'] = vmd.fit_sample
        for partition in ['reviewed_case_discovery', 'TNBC_source_holdout_sample_identity_pending']:
            part = vmd.partition.eq(partition).to_numpy();label = 'discovery' if partition == 'reviewed_case_discovery' else 'TNBC_holdout'
            prefix = inputs / (view + '_' + label)
            with gzip.open(str(prefix) + '_tpm.mtx.gz', 'wb', compresslevel=3) as f:mmwrite(f, lp[part][:, keep].T)
            vmd.loc[part].to_csv(str(prefix) + '_metadata.tsv', sep='\t', index=False)
            pd.DataFrame(pc[part], index=vmd.loc[part, 'cell'], columns=['PC%d' % (j + 1) for j in range(20)]).to_csv(str(prefix) + '_PC.tsv', sep='\t')
            pd.Series(genes[keep]).to_csv(str(prefix) + '_genes.tsv', sep='\t', index=False, header=False)
            summaries.append(dict(view=view, partition=partition, cells=int(part.sum()), case_or_sample_units=int(vmd.loc[part, 'sample'].nunique()),
                candidate_genes=len(keep), PCA_genes=len(hvg), PCs=20, discovery_covariate_rank=int(np.linalg.matrix_rank(design)), discovery_covariate_condition=float(np.linalg.cond(design))))
        np.savez_compressed(inputs / (view + '_frozen_projection.npz'), candidate_genes=genes[keep], hvg_genes=genes[hvg],
            pca_mean=pca.mean_, pca_components=pca.components_, pca_explained_variance=pca.explained_variance_,
            quality_mean=mu, quality_sd=sd, training_candidate_log_mean=mean[keep], training_candidate_log_variance=variance[keep],
            training_candidate_detection_fraction=detection_fraction[keep])
    md.to_csv(out / 'private_selected_native_cell_depth.tsv', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(out / 'private_input_scope_and_design_summary.tsv', sep='\t', index=False)
    files = [p for p in out.rglob('*') if p.is_file()];assert sum(p.stat().st_size for p in files) < 256 * 1024**2
    t = dict(status='COMPLETE_THREE_NATIVE_RNA_VIEWS_DISCOVERY_AND_TNBC_INPUTS_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), selected_cells=len(md),
        files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
        clinical_outcomes_used=False, new_MCP_trained=False)
    (out / 'private_terminal.json').write_text(json.dumps(t, indent=2) + '\n')
    print(json.dumps(dict(status=t['status'], completed_UTC=t['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
