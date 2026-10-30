"""Standard count-HVG source-local A, joint B, and native-coarse-routed C pilot.

Native expression stays in RAM. Same measured feature candidates are used for
all routes. Discovery uses all frozen sources; source-held-out reference mapping
is a separate required stage and is never claimed by this program.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata, pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import gather_csr_rows, LINEAGE_PANELS


ROUTES = {
    'T cell': 'T_NK', 'T/NK lymphoid': 'T_NK', 'T': 'T_NK', 'NK': 'T_NK',
    'T_NK_boundary': 'T_NK', 'B cell': 'B', 'B': 'B',
    'Plasma cell': 'Plasma', 'Plasma': 'Plasma',
    'Myeloid': 'Myeloid', 'Monocyte': 'Myeloid', 'Macrophage': 'Myeloid',
    'DC': 'Myeloid', 'Dendritic cell': 'Myeloid', 'Mast cell': 'Mast', 'Mast': 'Mast',
    'Neutrophil': 'Neutrophil', 'Erythroid': 'Erythroid',
    'Epithelial cell': 'Epithelial', 'Epithelial': 'Epithelial',
    'Endothelial cell': 'Endothelial', 'Endothelial': 'Endothelial',
    'Fibroblast': 'Fibroblast', 'Fibroblast cell': 'Fibroblast',
    'Pericyte': 'Mural', 'Smooth muscle cell': 'Mural', 'Mural': 'Mural',
    'Mural cell': 'Mural', 'Pericyte_SMC': 'Mural',
    'Oligodendrocyte': 'Glial',
}


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def normalized_log_counts(counts, native_totals):
    x = counts.tocsr().astype(np.float32)
    totals = np.asarray(native_totals, dtype=np.float64)
    if (x.shape[0] != len(totals) or not np.isfinite(totals).all() or
            (totals < 0).any() or not np.isfinite(x.data).all() or (x.data < 0).any() or
            not np.equal(x.data, np.floor(x.data)).all()):
        raise ValueError('Invalid native count domain or denominator')
    if np.any(np.asarray(x.sum(axis=1, dtype=np.float64)).ravel() > totals + 1e-5):
        raise ValueError('Feature sum exceeds native RNA denominator')
    positive = totals > 0
    scale = np.zeros(len(totals), dtype=np.float32)
    scale[positive] = 10000.0 / totals[positive]
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x, positive


def feature_columns(lookup, model, readout):
    columns = [next(iter(lookup[g])) if len(lookup.get(g, set())) == 1 else None for g in model]
    if any(c is None for c in columns) or len(set(columns)) != len(columns):
        raise ValueError('Model genes are missing, ambiguous or physical aliases')
    held = {next(iter(lookup[g])) for g in readout if len(lookup.get(g, set())) == 1}
    if set(columns) & held:
        raise ValueError('A held-out RNA readout aliases a physical model feature')
    return np.asarray(columns, dtype=np.int64)


from standard_annotation_order_primitives_v1 import (
    select_features, pca_fit, conditional_readout_metrics,
    coarse_routes_from_local, COARSE_ROUTING_GENES)
from TME_all_lineage_cluster_policy_v1 import GATES, GATE_ORDER, coherent_flags
from TME_assay_aware_T_anchor_policy_v1 import ACCESSORY, SIGNAL, CHAIN, assay_aware_flags


def harmonize(embedding, sources, seed, parameters):
    if len(np.unique(sources)) == 1:
        return embedding.copy(), {'status': 'one_source_no_correction'}
    import harmonypy
    fit = harmonypy.run_harmony(embedding, pd.DataFrame({'source': sources}), ['source'],
                               theta=parameters['theta'], lamb=parameters['lambda'],
                               max_iter_harmony=parameters['max_iter_harmony'],
                               ncores=parameters['threads'], random_state=seed, verbose=False)
    corrected = np.asarray(fit.Z_corr)
    if corrected.shape == embedding.T.shape:
        corrected = corrected.T
    if corrected.shape != embedding.shape or not np.isfinite(corrected).all():
        raise ValueError('Unexpected Harmony orientation or non-finite output')
    return corrected.astype(np.float32), {'status': 'discovery_all_sources_jointly_fit',
                                        'implementation_version': harmonypy.__version__}


def cluster(embedding, seed, parameters):
    import anndata as ad
    import scanpy as sc
    import numba
    sc.settings.n_jobs = parameters['threads']
    numba.set_num_threads(parameters['threads'])
    a = ad.AnnData(sp.csr_matrix((len(embedding), 0)))
    a.obsm['X_pilot'] = embedding
    sc.pp.neighbors(a, use_rep='X_pilot', n_neighbors=min(parameters['neighbors'], len(embedding) - 1),
                    random_state=seed, transformer='pynndescent', metric='euclidean')
    sc.tl.leiden(a, resolution=parameters['resolution'], random_state=seed,
                 flavor='igraph', n_iterations=2, directed=False)
    labels = a.obs['leiden'].to_numpy().astype(str)
    unique, codes = np.unique(labels, return_inverse=True)
    return codes.astype(np.int32), len(unique)



def marker_summary(labels, panel_log, detected, genes, measured):
    result = []
    for code in np.unique(labels):
        member = labels == code
        values = {}
        for name, panel in LINEAGE_PANELS.items():
            columns = [genes.index(g) for g in panel if g in genes and measured[genes.index(g)]]
            values[name] = {'n_measured': len(columns), 'n_declared': len(panel),
                            'mean_native_logCP10K': float(panel_log[member][:, columns].mean()) if columns else None,
                            'mean_detection_fraction': float(detected[member][:, columns].mean()) if columns else None}
        result.append({'partition_code': int(code), 'n_cells': int(member.sum()),
                       'lineage_panel_evidence': values, 'biological_identity_assigned': False})
    return result


def main(args):
    start = time.monotonic()
    root = Path(args.output)
    cfg = json.loads((root / 'config/model_spec_v1.json').read_text())
    pilot = Path(args.pilot)
    for role, item in cfg['bindings'].items():
        if sha256_file(item['path']) != item['sha256']:
            raise ValueError('Changed frozen input: ' + role)
    audit = json.loads((pilot / 'provenance/input_auditor_exit_v1.json').read_text())
    if audit['actual_exit_code'] != 0:
        raise ValueError('Pilot input audit not accepted')
    reader = FullContextAtlas(args.manifest)
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    plan = json.loads((pilot / 'config/model_features_and_heldouts_v1.json').read_text())
    from native_archive_collection_v1 import NativeArchiveCollection
    from native_annotation_review_v1 import native_measured_genes
    from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
    native_scope = NativeArchiveCollection(cfg['native_archive_manifest'])
    assert native_scope.archive_manifest['frozen_context']['sha256'] == sha256_file(args.manifest)
    native_metadata = {}
    common = set(plan['common_model_feature_candidates'])
    for entry in inventory:
        ds = entry['dataset']
        names, columns = native_measured_genes(native_scope, ds)
        native_metadata[ds] = dict(zip(names, map(int, columns)))
        common &= set(names)
        native_scope._fields.clear()
    model_genes = [g for g in plan['common_model_feature_candidates'] if g in common]
    heldout = {g for panels in plan['heldout_RNA_readouts'].values() for panel in panels.values() for g in panel}
    evidence_genes = sorted(heldout | {g for p in LINEAGE_PANELS.values() for g in p})
    pieces, raw_pieces, evidence, detected, totals, sources, routes, fine, measured_per_source = [], [], [], [], [], [], [], [], {}
    native_flags, source_assessed, source_modes = [], {}, {}
    wanted_typing = set(ACCESSORY + SIGNAL + CHAIN + ['CD3D','CD3E']) | {g for modules in GATES.values() for clause in modules for g in clause}
    import psutil
    assert psutil.virtual_memory().available >= cfg['minimum_available_RAM_bytes']
    sys.path.insert(0, str(Path(cfg['joined']) / 'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(root, cfg)
    smoke = __import__('subprocess').run([sys.executable, '-B', str(root / 'tests/test_standard_annotation_order_primitives_v1.py')], capture_output=True, text=True, timeout=180)
    (root / 'provenance/standard_API_test_stdout.txt').write_text(smoke.stdout)
    (root / 'provenance/standard_API_test_stderr.txt').write_text(smoke.stderr)
    assert smoke.returncode == 0 and 'skipped' not in smoke.stderr
    save(root / 'provenance/standard_API_test_receipt_v1.json', {'actual_exit_code': smoke.returncode, 'no_skips': True, 'synthetic_not_biological_validation': True})
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'),
             {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'stage': stage, 'wall_seconds': time.monotonic() - start, **extra})
        storage_gate(root, cfg)
    for row in inventory:
        ds = row['dataset']
        positions = np.load(pilot / row['selector_file'], allow_pickle=False)
        if sha256_file(pilot / row['selector_file']) != row['selector_sha256']:
            raise ValueError('Changed selected cells')
        lookup = reader.feature_lookup(ds)
        model_columns = feature_columns(lookup, model_genes, sorted(heldout))
        available = np.asarray([g in native_metadata[ds] and len(lookup.get(g, set())) == 1 for g in evidence_genes])
        columns = np.asarray([next(iter(lookup[g])) if available[i] else 0 for i, g in enumerate(evidence_genes)])
        measured_per_source[ds] = available
        local, local_counts, ev, det, native = [], [], [], [], []
        typing_genes = [g for g in native_metadata[ds] if g in wanted_typing]
        typing_columns = [native_metadata[ds][g] for g in typing_genes]
        flags_for_source = []
        z = reader.store(ds)
        key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        assert str(z.store.path) == str(native_scope.store(ds).store.path)
        delta = z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None
        for low in range(0, len(positions), cfg['read_block_rows']):
            chosen = positions[low:low + cfg['read_block_rows']]
            raw = gather_csr_rows(z[key], reader.rows(ds)[chosen])
            if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                    not np.equal(raw.data, np.floor(raw.data)).all()):
                raise ValueError('Native raw counts fail integer assumption')
            assert raw.data.max(initial=0) < 2**24
            total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
            local_counts.append(raw[:, model_columns].astype(np.float32).tocsr())
            x, positive = normalized_log_counts(raw[:, model_columns], total)
            analysis = raw.astype(np.float64).tocsr()
            if delta is not None:
                analysis = (analysis - gather_nonnegative_float_delta(delta, reader.rows(ds)[chosen])).tocsr()
            analysis.eliminate_zeros()
            assert np.isfinite(analysis.data).all() and analysis.data.min(initial=0) >= -1e-8
            analysis.data[analysis.data < 0] = 0
            analysis.eliminate_zeros()
            analysis_total = np.asarray(analysis.sum(axis=1)).ravel()
            panel = analysis[:, typing_columns].tocsr()
            flags, assessed = coherent_flags(panel, typing_genes)
            _, _, mode = assay_aware_flags(panel, typing_genes)
            flags_for_source.append(flags)
            source_assessed[ds] = assessed
            source_modes[ds] = mode
            y = analysis[:, columns].toarray().astype(np.float32)
            y[:, ~available] = 0
            fraction = np.zeros(len(total), dtype=np.float32)
            ok_evaluation = analysis_total > 0
            fraction[ok_evaluation] = 10000 / analysis_total[ok_evaluation]
            evaluated = np.log1p(y * fraction[:, None])
            evaluated[~ok_evaluation] = np.nan
            ev.append(evaluated)
            det.append(y > 0)
            native.append(total)
            local.append(x)
        labels = reader.frozen_labels(ds)
        major = labels['major'][positions].astype(str)
        pieces.append(sp.vstack(local, format='csr'))
        raw_pieces.append(sp.vstack(local_counts, format='csr'))
        native_flags.append(np.vstack(flags_for_source))
        evidence.append(np.vstack(ev)); detected.append(np.vstack(det)); totals.append(np.concatenate(native))
        sources.extend([ds] * len(positions))
        routes.extend([ROUTES.get(v, 'unresolved_context') for v in major])
        fine.extend(labels['fine'][positions].astype(str))
        checkpoint('native_pilot_loaded', dataset=ds, selected_cells=len(positions))
    x = sp.vstack(pieces, format='csr'); del pieces
    raw_x = sp.vstack(raw_pieces, format='csr'); del raw_pieces
    all_flags = np.vstack(native_flags); del native_flags
    ev = np.vstack(evidence); det = np.vstack(detected); del evidence, detected
    total = np.concatenate(totals)
    source = np.asarray(sources); route = np.asarray(routes); oldfine = np.asarray(fine)
    fit_rows = total > 0
    if fit_rows.sum() < len(total):
        checkpoint('zero_RNA_rows_preserved_but_not_modeled', n_excluded=int((~fit_rows).sum()))
    # These small position/label arrays are private, not exported to Git/chat.
    with (root / 'outputs/derived_axis_and_model_mask.npz').open('xb') as f:
        np.savez_compressed(f, modeled=fit_rows, source=source, inherited_route=route)
    hvg, hvg_table = select_features(raw_x[fit_rows], source[fit_rows], cfg['n_features'], model_genes, cfg['candidate_HVGs'])
    hvg_table.index.name = 'gene'
    hvg_table.to_csv(root / 'outputs/standard_batch_HVG_v1.tsv', sep='\t')
    del raw_x
    save(root / 'outputs/discovery_feature_selection_v1.json',
         {'all_source_discovery_not_holdout': True, 'method': 'official_Scanpy_seurat_v3_on_raw_integer_counts_source_batch', 'candidate_gene_universe': model_genes, 'reserved_RNA_not_in_model': True,
          'selected_feature_indices': hvg.tolist(), 'selected_feature_names': [model_genes[i] for i in hvg]})
    x = x[:, hvg]
    states, state_contract = {}, {}
    for name, genes in plan['heldout_RNA_readouts']['state'].items():
        used = [g for g in genes if g not in COARSE_ROUTING_GENES and all(
            measured_per_source[ds][evidence_genes.index(g)] for ds in measured_per_source)]
        values = np.full(len(total), np.nan, dtype=np.float32)
        if len(used) >= cfg['minimum_state_genes_measured']:
            values = ev[:, [evidence_genes.index(g) for g in used]].mean(axis=1)
        states[name] = values
        state_contract[name] = dict(used_gene_names=used, removed_coarse_routing_gene_names=sorted(set(genes)&COARSE_ROUTING_GENES),
            common_measured_same_signature_in_all_sources=True, model_and_coarse_routing_unused_RNA_not_orthogonal_assay=True)
    candidate = np.full(len(source), np.nan)
    for ds in source_assessed:
        if source_assessed[ds]['gamma_delta_T_RNA_candidate']:
            member = source == ds
            candidate[member] = all_flags[member, GATE_ORDER.index('gamma_delta_T_RNA_candidate')]
    save(root / 'outputs/evaluation_RNA_contract_v1.json', dict(state=state_contract,
        coarse_routing_gene_names=sorted(COARSE_ROUTING_GENES),source_assessed=source_assessed,source_T_modes=source_modes,
        gene_absence_not_negative=True,model_expression_layer='raw_integer_counts_full_RNA_CP10K_log1p',
        typing_and_evaluation_layer='inherited_adopted_analysis_counts_full_RNA_CP10K_log1p',
        gamma_RNA_overlap_with_coarse_routing_not_independent_identity_validation=True))
    with (root / 'outputs/private_shared_evaluation_vectors_v1.npz').open('xb') as f:
        np.savez_compressed(f, native_total=total, gamma_RNA_candidate=candidate, native_coherent_flags=all_flags,
                            **{'state_'+k:v for k,v in states.items()})
    anchors = {}
    rng = np.random.default_rng(cfg['diagnostic_sampling_seed'])
    for ds in np.unique(source):
        for r in np.unique(route):
            positions = np.flatnonzero((source == ds) & (route == r) & fit_rows)
            if len(positions) >= cfg['minimum_diagnostic_stratum']:
                anchors[(ds, r)] = np.sort(rng.choice(positions, min(len(positions), cfg['diagnostic_cells_per_stratum']), replace=False))
    results = []
    for seed in cfg['model_seeds']:
        partitions = {}
        global_raw, model = pca_fit(x[fit_rows], seed, cfg['components'], cfg['pca_iterations'])
        checkpoint('global_PCA_fit', seed=seed)
        global_corrected, harmony_info = harmonize(global_raw, source[fit_rows], seed, cfg)
        # Aligned A/B use precisely this same joint geometry; A annotations come
        # from local graphs and B from the joint graph. Geometry scores are one
        # observation, never two independent pieces of evidence.
        embeddings = {'raw_global_control': np.full((len(total), global_raw.shape[1]), np.nan, np.float32),
                      'A_B_shared_aligned': np.full((len(total), global_raw.shape[1]), np.nan, np.float32),
                      'A_local': np.full((len(total), cfg['components']), np.nan, np.float32),
                      'C_within_prior_route': np.full((len(total), cfg['components']), np.nan, np.float32)}
        embeddings['raw_global_control'][fit_rows] = global_raw
        embeddings['A_B_shared_aligned'][fit_rows] = global_corrected
        b, nb = cluster(global_corrected, seed, cfg)
        full = np.full(len(total), -1, np.int32); full[fit_rows] = b; partitions['B'] = full
        checkpoint('B_joint_graph_partitioned', seed=seed, n_partitions=nb)
        model_dir = root / 'models'; model_dir.mkdir(exist_ok=True)
        with (model_dir / ('global_discovery_PCA_' + str(seed) + '.npz')).open('xb') as f:
            np.savez_compressed(f, components=model.components_, mean=model.mean_, explained_variance=model.explained_variance_)
        for method, grouping, correct in [('A', source, False), ('C', route, True)]:
            labels = np.full(len(total), -1, np.int32)
            offset = 0
            if method == 'C':
                grouping = source_local_coarse_route
            for group in np.unique(grouping):
                member = np.flatnonzero(fit_rows & (grouping == group))
                if len(member) < cfg['minimum_fit_group']:
                    checkpoint('small_group_not_force_clustered', method=method, group=str(group), n_cells=len(member), seed=seed)
                    continue
                local, fitted = pca_fit(x[member], seed, cfg['components'], cfg['pca_iterations'])
                if correct:
                    local, _ = harmonize(local, source[member], seed, cfg)
                target = 'A_local' if method == 'A' else 'C_within_prior_route'
                embeddings[target][member, :local.shape[1]] = local
                # Small groups can have fewer PCs: pad missing coordinates with
                # zero, not NaN, without inventing extra fitted components.
                embeddings[target][member, local.shape[1]:] = 0
                lab, n = cluster(local, seed, cfg)
                labels[member] = lab + offset; offset += n
                checkpoint('local_graph_partitioned', method=method, group=str(group), n_cells=len(member), n_partitions=n, seed=seed)
            partitions[method] = labels
            if method == 'A':
                source_local_coarse_route, coarse_records = coarse_routes_from_local(source, labels, all_flags, source_assessed)
                save(root / 'outputs' / ('source_local_coarse_RNA_routes_'+str(seed)+'.json'), coarse_records)
        with (root / 'outputs' / ('private_embeddings_and_C_spaces_'+str(seed)+'.npz')).open('xb') as f:
            np.savez_compressed(f, C_space_id=source_local_coarse_route, **embeddings)
        summary = {'seed': seed, 'harmony': harmony_info, 'conditional_RNA_diagnostics': [],
                   'inherited_label_overlap_is_diagnostic_not_accuracy': True,
                   'C_routing_from_source_local_coherent_RNA_not_independent_truth': True, 'PCA_implementation': 'official_Scanpy_centered_ARPACK_no_gene_variance_scaling'}
        for (ds, r), fixed_positions in anchors.items():
            # C coordinates live in separate fitted coarse spaces. Every method
            # receives the identical intersected frozen anchors for comparison.
            for cr in np.unique(source_local_coarse_route[fixed_positions]):
                positions = fixed_positions[source_local_coarse_route[fixed_positions] == cr]
                for name, embedding in embeddings.items():
                    if not np.isfinite(embedding[positions]).all():
                        metrics = {'status': 'representation_unavailable_small_group'}
                    else:
                        spaces = source_local_coarse_route if name == 'C_within_prior_route' else source if name == 'A_local' else None
                        metrics = conditional_readout_metrics(embedding, states, total, candidate, positions, cfg['diagnostic_neighbors'], space_ids=spaces)
                    summary['conditional_RNA_diagnostics'].append(dict(source=str(ds), inherited_route=str(r),
                        source_local_C_RNA_route=str(cr), n_frozen_base_anchors=len(fixed_positions),
                        every_method_uses_same_anchor_subset=True, representation=name, **metrics))
        with (root / 'outputs' / ('candidate_partitions_' + str(seed) + '.npz')).open('xb') as f:
            np.savez_compressed(f, **partitions)
        tables = {}
        for method, labels in partitions.items():
            rows = []
            for code in np.unique(labels):
                member = labels == code
                rows.append({'partition_code': int(code), 'n_cells': int(member.sum()),
                             'inherited_route_counts': dict(collections.Counter(route[member])),
                             'source_counts': dict(collections.Counter(source[member])),
                             'inherited_fine_counts_not_truth': dict(collections.Counter(oldfine[member]))})
            tables[method] = rows
            # Source-specific measurement masks remain explicit; aggregate only
            # within a source. Do not treat source-unmeasured evidence as zero.
            for ds in np.unique(source):
                mask = (source == ds) & (labels >= 0) & fit_rows
                table = marker_summary(labels[mask], ev[mask], det[mask], evidence_genes, measured_per_source[ds])
                path = root / 'outputs' / ('native_partition_panels_' + str(seed) + '_' + method + '_' + ds + '.json.gz')
                with gzip.open(path, 'xt') as f:json.dump(table, f, allow_nan=False)
        save(root / 'outputs' / ('partition_crosswalk_diagnostic_' + str(seed) + '.json'), tables)
        save(root / 'outputs' / ('representation_diagnostics_' + str(seed) + '.json'), summary)
        results.append({'seed': seed, 'partitions': {m: len(set(v[v >= 0])) for m, v in partitions.items()}})
        checkpoint('seed_discovery_complete_not_annotation_acceptance', seed=seed)
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
            raise ValueError('Bounded pilot output cap exceeded')
    save(root / 'provenance/terminal_v1.json',
         {'status': 'STANDARD_NATIVE_COARSE_ROUTED_ABC_PILOT_COMPLETE_REQUIRES_INDEPENDENT_READBACK_BIOLOGICAL_REVIEW_AND_STUDY_HOLDOUT',
          'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'results': results,
          'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
          'native_counts_copied_to_disk': False, 'source_holdout_complete': False,
          'method_winner_selected': False, 'full_ABC_identity_validation_complete': False,
          'existing_authoritative_labels_or_CM_changed': False})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['output', 'pilot', 'manifest']:p.add_argument('--' + name, required=True)
    a = p.parse_args()
    try:main(a)
    except Exception as error:
        save(Path(a.output) / 'provenance/failure_v1.json',
             {'status': 'FAILED_NOT_ACCEPTED', 'error_class': type(error).__name__, 'error': str(error)})
        raise
