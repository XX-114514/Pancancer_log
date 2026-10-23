"""Fit source-local graphs on every qualified row, preserving archive exclusions.

The representative shortcut did not qualify for replacing local fine graphs.
This stage fits source-by-provisional-family anchors using native measured
features; common cross-source identities and native group annotation follow.
"""
import argparse
import datetime
import gc
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import psutil
import scipy.sparse as sp

from annotation_representation_pilot_v1 import normalized_log_counts, select_features, pca_fit, save
from compare_annotation_granularity_v1 import graph_grid
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_panel_aggregate_v1 import gather_csr_rows
from summarize_full_native_annotation_evidence_v1 import read_source


def save_arrays(path, **arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as f: np.savez_compressed(f, **arrays)
    with np.load(path, allow_pickle=False) as f:
        if set(arrays) != set(f.files): raise ValueError('Written array keys differ')
        for key, value in arrays.items(): np.testing.assert_array_equal(value, f[key])


def main(root):
    began = time.monotonic()
    cfg = json.loads((root / 'config/full_local_graph_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']: raise ValueError('Frozen local graph input changed')
    def checkpoint(stage, **values):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'), {
            'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage': stage, 'wall_seconds': time.monotonic() - began, **values})
    for filename, count in [('test_native_local_fit_pool_v1.py', 4), ('test_source_lineage_discovery_v1.py', 4)]:
        test = subprocess.run([sys.executable, '-B', str(root / 'tests' / filename)], capture_output=True, text=True, timeout=180)
        (root / 'provenance' / (filename + '.stderr.txt')).write_text(test.stderr)
        if test.returncode or ('Ran %d tests' % count) not in test.stderr or 'skipped' in test.stderr.lower():
            raise ValueError('Native branch counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 8, 'no_skips': True})
    inputs, reference = Path(cfg['local_inputs']), Path(cfg['reference'])
    for stage in [inputs, reference]:
        terminal = json.loads((stage / 'provenance/supervisor_terminal_v1.json').read_text())
        if terminal['actual_child_exit_code'] != 0 or terminal['wall_timeout']: raise ValueError('Input stage incomplete')
    pools = json.loads((inputs / 'outputs/full_native_local_discovery_input_index_v1.json').read_text())
    refs = {r['source']: r for r in json.loads((reference / 'outputs/full_archive_reference_index_v1.json').read_text())}
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    if set(refs) != {r['source'] for r in pools} or set(refs) != set(reader.archive_partitions): raise ValueError('Archive source scope differs')
    params = json.loads(Path(cfg['original_model_spec']).read_text())
    # Small sources first is an execution order, never biological selection.
    pools = sorted(pools, key=lambda item: (item['n_archive_rows'], item['source']))
    source_index = []
    total_modeled = total_fine = total_native_nnz = 0
    for item in pools:
        written = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        df = subprocess.run(['df', '-P', str(root)], capture_output=True, text=True, timeout=15)
        quota = subprocess.run(['quota', '-g', '-v', '-w', cfg['quota_group']], capture_output=True, text=True, timeout=15)
        if df.returncode or quota.returncode: raise RuntimeError('Fresh storage gate unavailable')
        filesystem = df.stdout.splitlines()[-1].split()[0]
        quota_rows = [line.split() for line in quota.stdout.splitlines() if line.split() and line.split()[0] == filesystem]
        if len(quota_rows) != 1: raise ValueError('Output filesystem quota is ambiguous')
        used, hard = int(quota_rows[0][1].rstrip('*')), int(quota_rows[0][3])
        fs = os.statvfs(root)
        available_bytes = min((hard - used) * 1024, fs.f_bavail * fs.f_frsize)
        required_bytes = max(0, cfg['output_cap_bytes'] - written) + cfg['other_workspace_reserve_bytes']
        if hard <= used or available_bytes < required_bytes: raise OSError('Fresh bounded output storage gate failed')
        source_started = time.monotonic()
        ds = item['source']; info, pool = read_source(inputs, item); rm, evidence = read_source(reference, refs[ds])
        rows = reader.rows(ds); n = len(rows)
        np.testing.assert_array_equal(rows, pool['native_row']); np.testing.assert_array_equal(rows, evidence['native_row'])
        identity_sha = text_axis_sha256(reader.native_ids(ds))
        if info['native_ids_sha256'] != identity_sha or rm['native_ids_sha256'] != identity_sha: raise ValueError('Native identity mismatch')
        feature_path = inputs / 'outputs' / info['feature_candidates']
        if sha256_file(feature_path) != info['feature_candidates_sha256']: raise ValueError('Native feature candidates changed')
        features = json.loads(feature_path.read_text())
        genes, columns = features['genes'], np.asarray(features['native_columns'], np.int64)
        if len(genes) != len(set(genes)) or len(columns) != len(set(columns)) or len(genes) != len(columns):
            raise ValueError('Native feature axes ambiguous')
        lookup = reader.feature_lookup(ds)
        if any(lookup.get(gene, set()) != {int(column)} for gene, column in zip(genes, columns)):
            raise ValueError('Source measured feature axis changed')
        if set(columns) & set(features['heldout_physical_alias_columns_excluded']): raise ValueError('RNA readout leaked into native model features')
        native = reader.store(ds)
        count_key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        groups = np.full(n, -1, np.int32)
        fine = np.zeros(n, bool)
        branches = []
        group_offset = 0
        source_nnz = 0
        # A bounded sparse source cache avoids rereading interleaved native CSR
        # spans for every lineage. It is RAM-only and does not cache all sources.
        fit_rows = np.flatnonzero(pool['fit_pool_candidate'])
        source_x = None
        cache_bound = (int(rm['native_nnz_read']) * 48 +
            len(fit_rows) * min(cfg['n_features'], len(genes)) * 12 + cfg['RAM_workspace_reserve_bytes'])
        if len(fit_rows) and cache_bound < .7 * psutil.virtual_memory().available:
            source_chunks = []
            for low in range(0, len(fit_rows), cfg['read_block_size']):
                selected_rows = fit_rows[low:low + cfg['read_block_size']]
                raw = gather_csr_rows(native[count_key], rows[selected_rows])
                if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                        not np.equal(raw.data, np.floor(raw.data)).all()):
                    raise ValueError('Native integer count domain changed')
                total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
                np.testing.assert_array_equal(total, evidence['native_total_RNA'][selected_rows])
                log, positive = normalized_log_counts(raw[:, columns], total)
                if not positive.all(): raise ValueError('Zero RNA entered source fit pool')
                source_chunks.append(log); source_nnz += raw.nnz
            if source_nnz > rm['native_nnz_read']: raise ValueError('Fitted subset exceeds frozen source native nnz')
            source_x = sp.vstack(source_chunks, format='csr'); del source_chunks, raw, log
            source_position = np.full(n, -1, np.int64); source_position[fit_rows] = np.arange(len(fit_rows))
            checkpoint('source_sparse_RAM_input_loaded', source=ds, n_fit_rows=len(fit_rows),
                       native_nnz_read=source_nnz, estimated_source_cache_peak_bytes=cache_bound)
            gc.collect()
        for family in np.unique(pool['provisional_family_code']):
            branch_mask = pool['provisional_family_code'] == family
            member = np.flatnonzero(branch_mask & pool['fit_pool_candidate'])
            branch_id = 'branch_%02d' % len(branches)
            record = {'branch_id': branch_id, 'provisional_family_code': int(family),
                'n_archive_rows': int(branch_mask.sum()), 'n_fit_pool_rows': len(member),
                'n_unfitted_archive_rows': int(branch_mask.sum()) - len(member)}
            if not len(member):
                record['status'] = 'no_qualified_fit_rows_archive_kept'
            elif len(member) < params['minimum_fit_group']:
                groups[member] = group_offset; group_offset += 1
                record.update(status='small_branch_coarse_only', n_groups=1)
            else:
                # Conservative dense-upper-bound CSR workspace + selected PCA.
                estimated_peak = (len(member) * len(genes) * 16 +
                    len(member) * min(cfg['n_features'], len(genes)) * 12 + cfg['RAM_workspace_reserve_bytes'])
                available = int(psutil.virtual_memory().available)
                if available < estimated_peak: raise MemoryError('Source-local branch failed fresh conservative RAM gate')
                checkpoint('native_branch_input_started', source=ds, branch_id=branch_id,
                           n_fit_rows=len(member), estimated_peak_RAM_bytes=estimated_peak, available_RAM_bytes=available)
                if source_x is not None:
                    if (source_position[member] < 0).any(): raise ValueError('Branch outside sparse source cache')
                    x = source_x[source_position[member]]
                else:
                    chunks = []
                    for low in range(0, len(member), cfg['read_block_size']):
                        selected_rows = member[low:low + cfg['read_block_size']]
                        raw = gather_csr_rows(native[count_key], rows[selected_rows])
                        if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                                not np.equal(raw.data, np.floor(raw.data)).all()):
                            raise ValueError('Native integer count domain changed')
                        total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
                        np.testing.assert_array_equal(total, evidence['native_total_RNA'][selected_rows])
                        log, positive = normalized_log_counts(raw[:, columns], total)
                        if not positive.all(): raise ValueError('Zero RNA entered fitted branch')
                        source_nnz += raw.nnz
                        chunks.append(log)
                    x = sp.vstack(chunks, format='csr'); del chunks, raw, log
                gc.collect()
                mean = np.asarray(x.mean(axis=0)).ravel()
                variance = np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel() - mean * mean, 0)
                variable = np.flatnonzero((variance > 0) & (np.asarray((x > 0).sum(axis=0)).ravel() >= 5))
                record['n_native_variable_candidates'] = len(variable)
                if len(variable) < 3:
                    groups[member] = group_offset; group_offset += 1
                    record.update(status='too_few_variable_features_coarse_only', n_groups=1)
                else:
                    candidate_x = x[:, variable]; del x
                    rank = select_features(candidate_x, np.zeros(len(member), np.int8), cfg['n_features'])
                    chosen = variable[rank]
                    model_x = candidate_x[:, rank]; del candidate_x; gc.collect()
                    embedding, model = pca_fit(model_x, cfg['seed'], params['components'], params['pca_iterations'])
                    del model_x; gc.collect()
                    local = graph_grid(embedding, cfg['seed'], params, [cfg['resolution']])[str(cfg['resolution'])]
                    groups[member] = local + group_offset; fine[member] = True
                    model_path = root / 'models' / ds / (branch_id + '.npz')
                    save_arrays(model_path, mean=model.mean_, components=model.components_, explained_variance=model.explained_variance_,
                        native_columns=columns[chosen], source_row_positions=member, native_rows=rows[member],
                        fit_embedding=embedding, source_local_group_code=groups[member])
                    model_meta_path = model_path.with_suffix('.json')
                    save(model_meta_path, {'source': ds, 'branch_id': branch_id, 'native_feature_names': [genes[j] for j in chosen],
                        'fit_axis_sha256': row_axis_sha256(rows[member]), 'source_local_group_offset': group_offset,
                        'feature_selection': 'within_source_branch_log_dispersion_rank_on_all_qualified_rows',
                        'source_local_transductive_discovery': True, 'clinical_or_CM_selection': False,
                        'reference_only_mapping_not_used_for_local_group_assignment': True})
                    occupied = len(np.unique(local)); group_offset += occupied
                    record.update(status='full_qualified_source_branch_graph', n_groups=occupied, n_selected_features=len(chosen),
                        model=str(model_path.relative_to(root)), model_sha256=sha256_file(model_path),
                        model_metadata=str(model_meta_path.relative_to(root)), model_metadata_sha256=sha256_file(model_meta_path),
                        estimated_peak_RAM_bytes=estimated_peak)
                    del embedding, model, local
                if 'x' in locals(): del x
                gc.collect()
            branches.append(record)
            checkpoint('full_source_local_branch_complete', source=ds, **record)
            if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
                raise ValueError('Full-source local model output cap exceeded')
        np.testing.assert_array_equal(groups >= 0, pool['fit_pool_candidate'])
        for group in np.unique(groups[groups >= 0]):
            if len(np.unique(pool['provisional_family_code'][groups == group])) != 1: raise ValueError('Local group crossed provisional routes')
        output = root / 'outputs' / (ds + '_native_local_groups_v1.npz')
        save_arrays(output, native_row=rows, source_local_group_code=groups, fine_graph_fitted=fine,
                    fit_pool_candidate=pool['fit_pool_candidate'], provisional_family_code=pool['provisional_family_code'])
        meta = {'source': ds, 'n_archive_rows': n, 'native_rows_sha256': row_axis_sha256(rows), 'native_ids_sha256': identity_sha,
            'n_fit_pool_rows': int(pool['fit_pool_candidate'].sum()), 'n_fine_graph_rows': int(fine.sum()),
            'n_source_local_groups': group_offset, 'branches': branches, 'native_nnz_read': source_nnz,
            'wall_seconds': time.monotonic() - source_started, 'local_group_codes_require_source_and_version_namespace': True,
            'local_identity_and_own_group_native_marker_review_pending': True, 'cross_source_common_identity_not_assigned': True,
            'original_labels_and_archive_exclusions_preserved': True, 'full_source_native_candidate_features_used': True,
            'used_bounded_source_sparse_RAM_cache': source_x is not None}
        meta_path = output.with_suffix('.json'); save(meta_path, meta)
        source_index.append({'source': ds, 'n_archive_rows': n, 'arrays': output.name, 'arrays_sha256': sha256_file(output),
                             'metadata': meta_path.name, 'metadata_sha256': sha256_file(meta_path)})
        total_modeled += int(pool['fit_pool_candidate'].sum()); total_fine += int(fine.sum()); total_native_nnz += source_nnz
        checkpoint('full_source_local_discovery_complete', source=ds, completed_sources=len(source_index))
        reader._fields.clear(); del evidence, pool, groups, fine, source_x; gc.collect()
    path = root / 'outputs/full_native_source_local_graph_index_v1.json'; save(path, source_index)
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_NATIVE_SOURCE_LOCAL_GRAPHS_COMPLETE_GROUP_ANNOTATION_AND_AUDIT_PENDING',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(source_index),
        'n_archive_rows': sum(r['n_archive_rows'] for r in source_index), 'n_modeled_rows': total_modeled,
        'n_fine_graph_rows': total_fine, 'native_nnz_read': total_native_nnz, 'index_sha256': sha256_file(path),
        'wall_seconds': time.monotonic() - began, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'native_counts_copied': False, 'new_final_identity': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); root = Path(parser.parse_args().output)
    try: main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
