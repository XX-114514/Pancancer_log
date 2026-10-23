"""Native own-group expression and unit support for every source-local group.

Consumes immutable, checkpointed source outputs while graph discovery runs.
Full completion requires the successful full producer and all expected sources;
final identity release additionally requires independent graph/evidence audits.
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
from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_annotation_review_v1 import native_measured_genes, native_group_evidence, descriptive_markers
from native_marker_unit_support_v1 import summarize_unit_direction
from native_panel_aggregate_v1 import gather_csr_rows
from source_local_graph_audit_rules_v1 import validate_partition
from summarize_full_native_annotation_evidence_v1 import read_source
from audit_full_source_local_graphs_v1 import index_at


def save_arrays(path, **arrays):
    with path.open('xb') as f: np.savez_compressed(f, **arrays)
    with np.load(path, allow_pickle=False) as f:
        if set(f.files) != set(arrays): raise ValueError('Native evidence serialization keys differ')
        for name, value in arrays.items(): np.testing.assert_array_equal(value, f[name])


def exact_process_live(root):
    path = root / 'provenance/child_launch_v1.json'
    if not path.exists(): return False
    value = json.loads(path.read_text()); stat = Path('/proc') / str(value['PID']) / 'stat'
    return stat.exists() and stat.read_text().rsplit(')', 1)[1].split()[19] == value['startticks']


def storage_gate(root, cfg):
    df = subprocess.run(['df', '-P', str(root)], capture_output=True, text=True, timeout=15)
    quota = subprocess.run(['quota', '-g', '-v', '-w', cfg['quota_group']], capture_output=True, text=True, timeout=15)
    if df.returncode or quota.returncode: raise OSError('Marker evidence storage quota unavailable')
    fsname = df.stdout.splitlines()[-1].split()[0]
    q = [line.split() for line in quota.stdout.splitlines() if line.split() and line.split()[0] == fsname]
    if len(q) != 1: raise ValueError('Ambiguous marker evidence quota')
    used, hard = int(q[0][1].rstrip('*')), int(q[0][3]); fs = os.statvfs(root)
    written = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    need = max(0, cfg['output_cap_bytes'] - written) + cfg['workspace_reserve_bytes']
    if hard <= used or min((hard - used) * 1024, fs.f_bavail * fs.f_frsize) < need:
        raise OSError('Marker evidence bounded storage gate failed')
    if written > cfg['output_cap_bytes']: raise ValueError('Marker evidence output cap exceeded')


def source_evidence(root, cfg, reader, producer, record, reference, ref_record):
    started = time.monotonic(); ds = record['source']; meta, groups = read_source(producer, record)
    rm, evidence = read_source(reference, ref_record); rows = reader.rows(ds); n = len(rows)
    for side in [groups, evidence]: np.testing.assert_array_equal(rows, side['native_row'])
    identity = text_axis_sha256(reader.native_ids(ds))
    if identity != meta['native_ids_sha256'] or identity != rm['native_ids_sha256']: raise ValueError('Marker evidence cell identity mismatch')
    codes = groups['source_local_group_code']; fitted = groups['fit_pool_candidate']; ng = meta['n_source_local_groups']
    validate_partition(rows, codes, groups['fine_graph_fitted'], fitted, groups['provisional_family_code'], ng)
    genes, columns = native_measured_genes(reader, ds); p = len(genes)
    if ng < 1: raise ValueError('No fitted local group for marker review')
    cache_bound = int(rm['native_nnz_read']) * 48 + ng * p * 48 + cfg['RAM_reserve_bytes']
    if ng * p * 64 + cfg['RAM_reserve_bytes'] > psutil.virtual_memory().available:
        raise MemoryError('Native group statistics failed fresh RAM gate')
    cache_ok = cache_bound < .6 * psutil.virtual_memory().available
    totals = evidence['native_total_RNA']; selected = np.flatnonzero(fitted); chunks = []
    stat = {'n': np.zeros(ng, np.int64), **{k: np.zeros((ng, p), np.float64) for k in ['cp_sum', 'log_sum', 'detected']}}
    native = reader.store(ds); key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
    nnz_read = 0
    for low in range(0, len(selected), cfg['read_block_size']):
        pos = selected[low:low + cfg['read_block_size']]; raw = gather_csr_rows(native[key], rows[pos])
        if not np.isfinite(raw.data).all() or (raw.data < 0).any() or not np.equal(raw.data, np.floor(raw.data)).all():
            raise ValueError('Invalid native marker count domain')
        np.testing.assert_array_equal(np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel(), totals[pos]); nnz_read += raw.nnz
        measured = raw[:, columns].tocsr(); aggregate = native_group_evidence(measured, genes, totals[pos], codes[pos], ng)
        for name in stat: stat[name] += aggregate[name]
        if cache_ok: chunks.append(measured)
    np.testing.assert_array_equal(stat['n'], np.bincount(codes[fitted], minlength=ng))
    if nnz_read > rm['native_nnz_read']: raise ValueError('Marker fit subset exceeds native source nnz')
    source_cache = sp.vstack(chunks, format='csr') if cache_ok else None
    del chunks, raw, measured, aggregate; gc.collect()
    family = np.full(ng, -999, np.int16)
    for group in range(ng): family[group] = groups['provisional_family_code'][np.flatnonzero(codes == group)[0]]
    group_markers = []
    for group in range(ng):
        same = (family == family[group]); same[group] = False
        count = int(stat['n'][same].sum()); scope = 'other_source_local_groups_in_same_provisional_family'
        if count == 0:
            same = np.arange(ng) != group; count = int(stat['n'][same].sum())
            scope = 'other_source_local_groups_all_families_no_within_family_comparator'
        if count:
            markers = descriptive_markers(genes, stat['cp_sum'][group] / stat['n'][group], stat['detected'][group] / stat['n'][group],
                                          stat['cp_sum'][same].sum(axis=0) / count, stat['detected'][same].sum(axis=0) / count,
                                          minimum_detection=cfg['minimum_marker_detection'], limit=cfg['marker_limit'])
        else: markers = []
        group_markers.append({'source_local_group_code': group, 'provisional_family_code': int(family[group]),
                              'n_fit_cells': int(stat['n'][group]), 'comparator_scope': scope,
                              'n_comparator_cells': count, 'own_group_markers': markers,
                              'no_own_marker_is_not_absence_of_identity': True})
    # A second RAM pass (or a bounded native pass when RAM caching is unsafe)
    # computes group-unit means and an explicitly same-family comparator.
    unit_values = reader.frozen_expression_units(ds).astype(str)
    units, unit = np.unique(unit_values[selected], return_inverse=True); gcodes = codes[selected]
    pairkeys, paircode, pair_n = np.unique(unit.astype(np.int64) * ng + gcodes, return_inverse=True, return_counts=True)
    pair_group, pair_unit = pairkeys % ng, pairkeys // ng
    k = cfg['marker_limit']; marker_index = np.full((ng, k), -1, np.int64)
    gene_index = {g: i for i, g in enumerate(genes)}
    for group, info in enumerate(group_markers):
        for j, m in enumerate(info['own_group_markers']): marker_index[group, j] = gene_index[m['gene']]
    own_cp = np.zeros((len(pairkeys), k), np.float64); own_det = np.zeros_like(own_cp)
    family_specs = {}
    for f in np.unique(family):
        members = np.flatnonzero(family == f); marker_genes = np.unique(marker_index[members][marker_index[members] >= 0])
        n_by_unit = np.bincount(unit[np.isin(gcodes, members)], minlength=len(units))
        family_specs[int(f)] = (marker_genes, n_by_unit)
    support_bytes = sum(len(units) * len(v[0]) * 16 for v in family_specs.values()) + own_cp.nbytes * 2
    if support_bytes + cfg['RAM_reserve_bytes'] > psutil.virtual_memory().available: raise MemoryError('Marker unit support exceeds RAM gate')
    family_stats = {f: {'genes': g, 'n': ns, 'cp': np.zeros((len(units), len(g)), np.float64),
                        'detected': np.zeros((len(units), len(g)), np.float64)} for f, (g, ns) in family_specs.items()}
    for low in range(0, len(selected), cfg['read_block_size']):
        pos = selected[low:low + cfg['read_block_size']]; blockcodes = codes[pos]; ub = unit[low:low + len(pos)]; pb = paircode[low:low + len(pos)]
        if source_cache is None:
            raw = gather_csr_rows(native[key], rows[pos]); nnz_read += raw.nnz
            np.testing.assert_array_equal(np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel(), totals[pos])
            block = raw[:, columns].tocsr()
        else: block = source_cache[low:low + len(pos)]
        scale = 10000. / totals[pos]
        for f, context in family_stats.items():
            member = family[blockcodes] == f
            if not member.any() or not len(context['genes']): continue
            x = block[member][:, context['genes']]
            cp = x.astype(np.float64).multiply(scale[member, None]).toarray()
            detected = (x > 0).toarray().astype(np.float64)
            np.add.at(context['cp'], ub[member], cp); np.add.at(context['detected'], ub[member], detected)
        for group in np.unique(blockcodes):
            member = blockcodes == group; cols = marker_index[group]; keep = cols >= 0
            if not keep.any(): continue
            x = block[member][:, cols[keep]]
            cp = x.astype(np.float64).multiply(scale[member, None]).toarray(); det = (x > 0).toarray().astype(np.float64)
            np.add.at(own_cp, (pb[member, None], np.flatnonzero(keep)[None, :]), cp)
            np.add.at(own_det, (pb[member, None], np.flatnonzero(keep)[None, :]), det)
    support = []
    for group, info in enumerate(group_markers):
        pairs = np.flatnonzero(pair_group == group); u = pair_unit[pairs]; context = family_stats[int(family[group])]
        parent_n = context['n'][u]; other_n = parent_n - pair_n[pairs]
        if (other_n < 0).any(): raise ValueError('Negative within-family comparator size')
        for j, marker in enumerate(info['own_group_markers']):
            col = int(np.searchsorted(context['genes'], marker_index[group, j]))
            other_cp = context['cp'][u, col] - own_cp[pairs, j]
            if (other_cp < -1e-5).any(): raise ValueError('Negative within-family comparator expression')
            support.append({'source_local_group_code': group, 'gene': marker['gene'],
                            'comparator_scope': 'same_provisional_family_other_local_groups_within_same_expression_unit',
                            **summarize_unit_direction(pair_n[pairs], own_cp[pairs, j], other_n, np.maximum(other_cp, 0), cfg['minimum_unit_cells'])})
    arrays = root / 'outputs' / (ds + '_local_native_expression_v1.npz')
    save_arrays(arrays, n=stat['n'], cp_sum=stat['cp_sum'], log_sum=stat['log_sum'], detected=stat['detected'],
                native_columns=columns, provisional_family_code=family)
    units_path = root / 'outputs' / (ds + '_local_marker_unit_evidence_v1.npz')
    save_arrays(units_path, pair_group_code=pair_group, expression_unit_code=pair_unit, n_cells=pair_n,
                marker_gene_index=marker_index, cp_sum=own_cp, detected=own_det)
    output = root / 'outputs' / (ds + '_local_native_marker_evidence_v1.json')
    save(output, {'source': ds, 'native_rows_sha256': row_axis_sha256(rows), 'native_ids_sha256': identity,
                  'n_archive_rows': n, 'n_fitted_rows': len(selected), 'n_unfitted_archive_rows_retained': n - len(selected),
                  'n_local_groups': ng, 'native_measured_gene_names': genes, 'all_native_measured_genes_included_in_evidence': True,
                  'training_excluded_RNA_genes_included_in_annotation_not_independent_validation': True,
                  'group_native_statistics': arrays.name, 'group_native_statistics_sha256': sha256_file(arrays),
                  'marker_unit_arrays': units_path.name, 'marker_unit_arrays_sha256': sha256_file(units_path),
                  'expression_unit_dictionary': units.tolist(), 'expression_unit_dictionary_sha256': text_axis_sha256(units),
                  'groups': group_markers, 'marker_unit_direction_summaries': support,
                  'source_graph_record': record, 'used_bounded_source_sparse_RAM_cache': cache_ok, 'native_nnz_read': int(nnz_read),
                  'wall_seconds': time.monotonic() - started, 'new_identity_assigned': False, 'inferential_DEG_computed': False,
                  'parent_group_markers_borrowed': False, 'patient_independence_claimed': False,
                  'graph_and_native_evidence_independent_acceptance_required_before_identity_release': True})
    return {'source': ds, 'n_archive_rows': n, 'n_fitted_rows': len(selected), 'n_local_groups': ng,
            'metadata': output.name, 'metadata_sha256': sha256_file(output), 'statistics': arrays.name,
            'statistics_sha256': sha256_file(arrays), 'marker_unit_arrays': units_path.name,
            'marker_unit_arrays_sha256': sha256_file(units_path), 'source_graph_record': record}


def main(root):
    began = time.monotonic(); cfg = json.loads((root / 'config/native_marker_evidence_spec_v1.json').read_text())
    producer, reference = Path(cfg['producer']), Path(cfg['reference'])
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']: raise ValueError('Native marker input binding changed')
    test_specs = [(['test_native_marker_unit_support_v1.py'], 3),
                  (['test_native_annotation_review_v1.py',
                    'NativeReviewTests.test_native_total_zero_RNA_and_group_denominators',
                    'NativeReviewTests.test_marker_fallback_is_own_group_expression_not_DEG'], 2)]
    for arguments, count in test_specs:
        test = subprocess.run([sys.executable, '-B', str(root / 'tests' / arguments[0])] + arguments[1:],
                              capture_output=True, text=True, timeout=120)
        (root / 'provenance' / (arguments[0] + '.stderr.txt')).write_text(test.stderr)
        if test.returncode or ('Ran %d tests' % count) not in test.stderr or 'skipped' in test.stderr.lower():
            raise ValueError('Native marker numerical counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 5, 'no_skips': True})
    refs = index_at(reference, 'full_archive_reference_index_v1.json')
    reader = NativeArchiveCollection(cfg['archive_manifest']); expected = set(reader.archive_partitions)
    if set(refs) != expected: raise ValueError('Native marker full scope differs')
    results = {}
    while set(results) != expected:
        complete = set()
        for checkpoint in (producer / 'checkpoints').glob('*.json'):
            try: value = json.loads(checkpoint.read_text())
            except json.JSONDecodeError:
                if time.time() - checkpoint.stat().st_mtime < 5: continue
                raise
            if value.get('stage') == 'full_source_local_discovery_complete': complete.add(value['source'])
        if not complete <= expected: raise ValueError('Unregistered source graph output')
        pending = sorted(complete - set(results), key=lambda ds: (reader.archive_partitions[ds]['n_archive_rows'], ds))
        for ds in pending:
            storage_gate(root, cfg)
            meta_path = producer / 'outputs' / (ds + '_native_local_groups_v1.json'); array_path = meta_path.with_suffix('.npz')
            meta = json.loads(meta_path.read_text())
            record = {'source': ds, 'n_archive_rows': meta['n_archive_rows'], 'arrays': array_path.name,
                      'arrays_sha256': sha256_file(array_path), 'metadata': meta_path.name, 'metadata_sha256': sha256_file(meta_path)}
            results[ds] = source_evidence(root, cfg, reader, producer, record, reference, refs[ds])
            save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'source': ds, 'completed_sources': len(results), 'wall_seconds': time.monotonic() - began})
            reader._fields.clear(); gc.collect()
        if set(results) != expected:
            terminal = producer / 'provenance/supervisor_terminal_v1.json'
            if terminal.exists():
                value = json.loads(terminal.read_text())
                if value['actual_child_exit_code'] or value['wall_timeout']: raise RuntimeError('Upstream graph producer failed; partial evidence preserved')
                if not pending and not (expected <= complete): raise ValueError('Successful graph producer has incomplete source checkpoints')
            elif not exact_process_live(producer): raise RuntimeError('Upstream graph process identity unavailable; partial evidence preserved')
            if time.monotonic() - began > cfg['maximum_producer_wait_seconds']: raise TimeoutError('Bounded graph/evidence time exhausted')
            time.sleep(15)
    while not (producer / 'provenance/supervisor_terminal_v1.json').exists():
        if not exact_process_live(producer):
            time.sleep(2)
            if not (producer / 'provenance/supervisor_terminal_v1.json').exists(): raise RuntimeError('Missing final graph producer supervision')
        if time.monotonic() - began > cfg['maximum_producer_wait_seconds']: raise TimeoutError('Bounded final graph supervision wait exhausted')
        time.sleep(5)
    full_graph = index_at(producer, 'full_native_source_local_graph_index_v1.json')
    if set(full_graph) != expected: raise ValueError('Final graph scope differs')
    for ds in sorted(expected):
        if full_graph[ds] != results[ds]['source_graph_record']: raise ValueError('Completed-source snapshot differs from final graph index')
    storage_gate(root, cfg)
    output = root / 'outputs/full_source_local_marker_evidence_index_v1.json'; save(output, [results[ds] for ds in sorted(expected)])
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_NATIVE_LOCAL_MARKER_EVIDENCE_COMPLETE_INDEPENDENT_AUDIT_PENDING',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(results),
        'n_archive_rows': sum(r['n_archive_rows'] for r in results.values()), 'n_fitted_rows': sum(r['n_fitted_rows'] for r in results.values()),
        'n_local_groups': sum(r['n_local_groups'] for r in results.values()), 'index_sha256': sha256_file(output),
        'graph_producer_terminal_sha256': sha256_file(producer / 'provenance/terminal_v1.json'),
        'wall_seconds': time.monotonic() - began, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'new_final_identity': False, 'inferential_DEG_computed': False, 'native_counts_copied': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); root = Path(parser.parse_args().output)
    try: main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
