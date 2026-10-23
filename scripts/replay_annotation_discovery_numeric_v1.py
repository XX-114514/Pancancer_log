"""Independent native/panel and representative neighborhood numerical replay.

Reads CSR arrays directly and reruns the declared PCA independently. No raw
matrix is saved. The smallest frozen source is selected by size, not its result.
"""
import argparse
import datetime
import gzip
import json
from pathlib import Path
import time

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import LINEAGE_PANELS


def read_rows_direct(csr, requested):
    requested = np.asarray(requested, dtype=np.int64)
    unique, inverse = np.unique(requested, return_inverse=True)
    pieces = []
    low = 0
    width = int(csr.attrs['shape'][1])
    while low < len(unique):
        high = np.searchsorted(unique, unique[low] + 8192)
        begin, end = int(unique[low]), int(unique[high - 1]) + 1
        pointers = np.asarray(csr['indptr'][begin:end + 1], dtype=np.int64)
        first, last = int(pointers[0]), int(pointers[-1])
        block = sp.csr_matrix((np.asarray(csr['data'][first:last]),
                               np.asarray(csr['indices'][first:last]), pointers - first),
                              shape=(end - begin, width))
        pieces.append(block[unique[low:high] - begin]);low = high
    return sp.vstack(pieces, format='csr')[inverse]


def main(args):
    started = time.monotonic();run, pilot, out = Path(args.run), Path(args.pilot), Path(args.output)
    cfg = json.loads((run / 'config/model_spec_v1.json').read_text())
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    features = json.loads((pilot / 'config/model_features_and_heldouts_v1.json').read_text())
    selected = json.loads((run / 'outputs/discovery_feature_selection_v1.json').read_text())['selected_feature_names']
    readout = features['heldout_RNA_readouts']
    evidence_genes = sorted({g for panels in readout.values() for panel in panels.values() for g in panel} |
                            {g for panel in LINEAGE_PANELS.values() for g in panel})
    chosen = min(inventory, key=lambda item: item['n_selected_cells'])['dataset']
    reader = FullContextAtlas(args.manifest)
    with np.load(run / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        source, route, modeled = f['source'].copy(), f['inherited_route'].copy(), f['modeled'].copy()
    blocks, totals, source_axis = [], [], []
    evidence, detection, available = None, None, None
    for item in inventory:
        ds = item['dataset'];selector = pilot / item['selector_file']
        if sha256_file(selector) != item['selector_sha256']:raise ValueError('Changed frozen selector')
        positions = np.load(selector, allow_pickle=False)
        lookup = reader.feature_lookup(ds)
        if any(len(lookup.get(g, set())) != 1 for g in selected):raise ValueError('Changed measured PCA feature')
        columns = np.asarray([next(iter(lookup[g])) for g in selected], dtype=np.int64)
        z = reader.store(ds)
        key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        if ds == chosen:
            available = np.array([len(lookup.get(g, set())) == 1 for g in evidence_genes])
            ev_cols = np.array([next(iter(lookup[g])) if available[j] else 0 for j, g in enumerate(evidence_genes)], dtype=np.int64)
        local_x, local_t, local_y, local_d = [], [], [], []
        for low in range(0, len(positions), 2048):
            rows = reader.rows(ds)[positions[low:low + 2048]]
            raw = read_rows_direct(z[key], rows)
            if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                    not np.equal(raw.data, np.floor(raw.data)).all()):raise ValueError('Native count domain changed')
            total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
            scale = np.zeros(len(total), np.float32);scale[total > 0] = 10000 / total[total > 0]
            x = raw[:, columns].astype(np.float32).multiply(scale[:, None]).tocsr()
            x.data = np.log1p(x.data);local_x.append(x);local_t.append(total)
            if ds == chosen:
                y = raw[:, ev_cols].toarray().astype(np.float32);y[:, ~available] = 0
                local_d.append(y > 0);local_y.append(np.log1p(y * scale[:, None]))
        blocks.append(sp.vstack(local_x, format='csr'));totals.append(np.concatenate(local_t));source_axis.extend([ds] * len(positions))
        if ds == chosen:evidence, detection = np.vstack(local_y), np.vstack(local_d)
    x, native_total = sp.vstack(blocks, format='csr'), np.concatenate(totals)
    if not np.array_equal(source, np.asarray(source_axis)) or not np.array_equal(modeled, native_total > 0):
        raise ValueError('Independent count/source/model eligibility replay differs')
    source_rows = np.flatnonzero(source == chosen)
    panel_values_checked, largest_panel_error = 0, 0.
    for seed in cfg['model_seeds']:
        with np.load(run / 'outputs' / ('candidate_partitions_' + str(seed) + '.npz'), allow_pickle=False) as f:
            for method in ['A', 'B', 'C']:
                labels = f[method][source_rows]
                path = run / 'outputs' / ('native_partition_panels_' + str(seed) + '_' + method + '_' + chosen + '.json.gz')
                with gzip.open(path, 'rt') as handle:table = json.load(handle)
                for group in table:
                    member = labels == group['partition_code']
                    if int(member.sum()) != group['n_cells']:raise ValueError('Independent panel group size differs')
                    for name, genes in LINEAGE_PANELS.items():
                        cols = [evidence_genes.index(g) for g in genes if available[evidence_genes.index(g)]]
                        original = group['lineage_panel_evidence'][name]
                        if len(cols) != original['n_measured']:raise ValueError('Independent panel mask differs')
                        if cols:
                            mean = float(evidence[member][:, cols].astype(np.float64).mean())
                            detected = float(detection[member][:, cols].mean())
                            np.testing.assert_allclose(mean, original['mean_native_logCP10K'], rtol=2e-6, atol=2e-6)
                            np.testing.assert_allclose(detected, original['mean_detection_fraction'], rtol=1e-12, atol=1e-12)
                            largest_panel_error = max(largest_panel_error, abs(mean - original['mean_native_logCP10K']))
                            panel_values_checked += 2
    seed = cfg['model_seeds'][0]
    pca = PCA(n_components=cfg['components'], svd_solver='randomized', iterated_power=cfg['pca_iterations'],
              random_state=seed, copy=False)
    global_embedding = pca.fit_transform(x[modeled].toarray().astype(np.float32, copy=False))
    with np.load(run / 'models' / ('global_discovery_PCA_' + str(seed) + '.npz'), allow_pickle=False) as original:
        for a, b in [(pca.mean_, original['mean']), (pca.components_, original['components']),
                     (pca.explained_variance_, original['explained_variance'])]:
            np.testing.assert_allclose(a, b, rtol=2e-5, atol=2e-5)
    if not modeled.all():raise ValueError('Representative replay requires audited positive-RNA pilot; do not invent omitted row positions')
    local_pca = PCA(n_components=cfg['components'], svd_solver='randomized', iterated_power=cfg['pca_iterations'],
                    random_state=seed, copy=False)
    local_embedding = local_pca.fit_transform(x[source_rows].toarray().astype(np.float32, copy=False))
    rng = np.random.default_rng(cfg['diagnostic_sampling_seed']);anchors = None
    for ds in np.unique(source):
        for parent in np.unique(route):
            positions = np.flatnonzero((source == ds) & (route == parent) & modeled)
            if len(positions) >= cfg['minimum_diagnostic_stratum']:
                selected_rows = np.sort(rng.choice(positions, min(len(positions), cfg['diagnostic_cells_per_stratum']), replace=False))
                if ds == chosen and parent == 'T_NK':anchors = selected_rows
    if anchors is None:raise ValueError('Size-selected source lacks the prespecified T/NK replay scope')
    local_positions = np.searchsorted(source_rows, anchors)
    source_values = {}
    for name, panel in readout['state'].items():
        cols = [evidence_genes.index(g) for g in panel if available[evidence_genes.index(g)]]
        if len(cols) >= cfg['minimum_state_genes_measured']:source_values[name] = evidence[:, cols].mean(axis=1)[local_positions]
    original = json.loads((run / 'outputs' / ('representation_diagnostics_' + str(seed) + '.json')).read_text())
    comparisons, largest_correlation_error = [], 0.
    for representation, embedding in [('raw_global_control', global_embedding[anchors]), ('A_local', local_embedding[local_positions])]:
        _, near = NearestNeighbors(n_neighbors=cfg['diagnostic_neighbors'] + 1, n_jobs=4).fit(embedding).kneighbors(embedding)
        near = np.asarray([indices[indices != j][:cfg['diagnostic_neighbors']] for j, indices in enumerate(near)])
        expected = next(r for r in original['conditional_RNA_diagnostics']
                        if r['source'] == chosen and r['inherited_route'] == 'T_NK' and r['representation'] == representation)
        depth = np.log1p(native_total[anchors]);design = np.column_stack([np.ones(len(depth)), depth - depth.mean()])
        for program, y in source_values.items():
            residual = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
            predicted = residual[near].mean(axis=1)
            old = expected['states'][program]
            if residual.std() < 1e-8 or predicted.std() < 1e-8:
                if old['status'] != 'no_variable_signal':raise ValueError('No-variable RNA case differs')
                continue
            if old['status'] != 'estimated':raise ValueError('Estimable RNA case differs')
            actual = {'depth_adjusted_neighbor_Pearson': float(np.corrcoef(residual, predicted)[0, 1]),
                      'depth_adjusted_neighbor_Spearman': float(np.corrcoef(rankdata(residual), rankdata(predicted))[0, 1])}
            for key, value in actual.items():
                np.testing.assert_allclose(value, old[key], rtol=1e-6, atol=1e-6)
                largest_correlation_error = max(largest_correlation_error, abs(value - old[key]))
            comparisons.append({'representation': representation, 'program': program, 'passed': True})
    result = {'status': 'NATIVE_PANEL_AND_REPRESENTATIVE_RNA_NUMERIC_REPLAY_PASS',
              'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_selected_by_smallest_frozen_cell_count_not_result': chosen,
              'all_three_seed_all_method_panels_in_selected_source_replayed': True,
              'n_native_panel_values_checked': panel_values_checked,
              'largest_native_panel_float32_float64_error': largest_panel_error,
              'full_pilot_global_PCA_recomputed_and_parameters_matched': True,
              'source_local_PCA_recomputed': True, 'representative_RNA_cases': comparisons,
              'largest_RNA_correlation_replay_error': largest_correlation_error,
              'B_C_corrected_neighborhood_numeric_replay': 'not_claimed_by_this_representative_replay',
              'new_identity_or_method_winner_accepted': False,
              'wall_seconds': time.monotonic() - started}
    with (out / 'provenance/terminal_v1.json').open('x') as f:json.dump(result, f, indent=2);f.write('\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['run', 'pilot', 'manifest', 'output']:p.add_argument('--' + name, required=True)
    main(p.parse_args())
