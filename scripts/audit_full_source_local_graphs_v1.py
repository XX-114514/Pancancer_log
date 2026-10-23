"""Independent axis/native/PCA/graph audit, full or frozen completed-source scope.

Replays one bounded, predetermined branch numerically and checks every saved
model axis. This establishes computational consistency, not identity truth.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_lineage_routing_v1 import FAMILIES
from source_local_graph_audit_rules_v1 import validate_partition, validate_model_axes
from summarize_full_native_annotation_evidence_v1 import read_source


def index_at(root, name):
    sup = json.loads((root / 'provenance/supervisor_terminal_v1.json').read_text())
    terminal = json.loads((root / 'provenance/terminal_v1.json').read_text())
    path = root / 'outputs' / name
    if sup['actual_child_exit_code'] or sup['wall_timeout'] or sha256_file(path) != terminal['index_sha256']:
        raise ValueError('Upstream index not bound to successful producer')
    records = json.loads(path.read_text())
    if len({r['source'] for r in records}) != len(records): raise ValueError('Duplicate source index')
    return {r['source']: r for r in records}


def direct_native_rows(group, rows):
    data, indices, ptr = [], [], [0]
    for row in rows:
        low, high = map(int, group['indptr'][int(row):int(row) + 2])
        value, col = np.asarray(group['data'][low:high]), np.asarray(group['indices'][low:high])
        if not np.isfinite(value).all() or (value < 0).any() or not np.equal(value, np.floor(value)).all():
            raise ValueError('Invalid native count domain')
        data.append(value); indices.append(col); ptr.append(ptr[-1] + len(value))
    return sp.csr_matrix((np.concatenate(data), np.concatenate(indices), np.asarray(ptr)),
                         shape=(len(rows), int(group.attrs['shape'][1])))


def replay_branch(reader, ds, rows, member, genes, columns, expected_totals, stored, params, cfg):
    native = reader.store(ds)
    key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
    raw = direct_native_rows(native[key], rows[member])
    totals = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
    np.testing.assert_array_equal(totals, expected_totals)
    x = raw[:, columns].astype(np.float32).tocsr()
    scale = np.asarray(10000. / totals, np.float32)
    x = x.multiply(scale[:, None]).tocsr(); x.data = np.log1p(x.data)
    mean = np.asarray(x.mean(axis=0)).ravel()
    var = np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel() - mean ** 2, 0)
    detected = np.asarray((x > 0).sum(axis=0)).ravel()
    variable = np.flatnonzero((var > 0) & (detected >= 5))
    # Recompute after slicing, matching float32 sparse summation domain.
    v = x[:, variable]
    vm = np.asarray(v.mean(axis=0)).ravel()
    vv = np.maximum(np.asarray(v.multiply(v).mean(axis=0)).ravel() - vm ** 2, 0)
    score = vv / (vm + 1e-3)
    # Average ranks followed by lexicographic tie-breaking are equivalent to
    # descending dispersion then the original feature index for one source.
    chosen = variable[np.lexsort((np.arange(len(variable)), -score))[:min(cfg['n_features'], len(variable))]]
    np.testing.assert_array_equal(columns[chosen], stored['native_columns'])
    dense = x[:, chosen].toarray().astype(np.float32, copy=False)
    model = PCA(n_components=min(params['components'], len(member) - 1, dense.shape[1] - 1),
                svd_solver='randomized', random_state=cfg['seed'], iterated_power=params['pca_iterations'], copy=False)
    embedding = model.fit_transform(dense).astype(np.float32)
    sign = np.sign(np.sum(model.components_ * stored['components'], axis=1)); sign[sign == 0] = 1
    np.testing.assert_allclose(model.mean_, stored['mean'], rtol=2e-5, atol=2e-5)
    np.testing.assert_allclose(model.components_ * sign[:, None], stored['components'], rtol=2e-4, atol=2e-4)
    np.testing.assert_allclose(embedding * sign, stored['fit_embedding'], rtol=2e-4, atol=2e-3)
    np.testing.assert_allclose(model.explained_variance_, stored['explained_variance'], rtol=2e-4, atol=2e-4)
    # A fresh graph call uses the saved fit-transform coordinates. Randomized
    # fit_transform is not assumed identical to a later model.transform.
    import anndata as ad
    import scanpy as sc
    import numba
    sc.settings.n_jobs = params['threads']; numba.set_num_threads(params['threads'])
    obj = ad.AnnData(np.zeros((len(member), 1), np.float32)); obj.obsm['X_pca'] = stored['fit_embedding']
    sc.pp.neighbors(obj, n_neighbors=min(params['neighbors'], len(member) - 1), use_rep='X_pca',
                    metric='euclidean', random_state=cfg['seed'], transformer='pynndescent')
    sc.tl.leiden(obj, resolution=cfg['resolution'], random_state=cfg['seed'], flavor='igraph',
                 directed=False, n_iterations=2)
    ari = adjusted_rand_score(stored['source_local_group_code'], obj.obs['leiden'].to_numpy())
    if ari != 1.: raise ValueError('Independent graph partition not reproduced')
    return {'source': ds, 'n_full_branch_rows_replayed': len(member), 'graph_ARI': ari,
            'max_embedding_absolute_error': float(np.max(np.abs(embedding * sign - stored['fit_embedding']))),
            'native_feature_selection_reproduced': True, 'not_independent_biological_identity_truth': True}


def main(root):
    start = time.monotonic()
    spec = json.loads((root / 'config/local_graph_audit_spec_v1.json').read_text())
    producer = Path(spec['producer']); cfg_path = producer / 'config/full_local_graph_spec_v1.json'
    if sha256_file(cfg_path) != spec['producer_config_sha256']: raise ValueError('Producer configuration changed')
    cfg = json.loads(cfg_path.read_text())
    request = json.loads((producer / 'config/supervisor_request_v1.json').read_text())
    for path, digest in request['script_SHA256'].items():
        if sha256_file(producer / path) != digest: raise ValueError('Producer deployed code changed')
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']: raise ValueError('Producer input binding changed')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_source_local_graph_audit_rules_v1.py')],
                          capture_output=True, text=True, timeout=90)
    (root / 'provenance/numerical_tests.stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Source graph audit counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 4, 'no_skips': True})
    inputs, reference = Path(cfg['local_inputs']), Path(cfg['reference'])
    pools = index_at(inputs, 'full_native_local_discovery_input_index_v1.json')
    refs = index_at(reference, 'full_archive_reference_index_v1.json')
    full = spec['scope'] == 'full_completed_archive'
    if full:
        sources = index_at(producer, 'full_native_source_local_graph_index_v1.json')
    else:
        sources = {r['source']: r for r in spec['completed_sources']}
        if len(sources) != len(spec['completed_sources']) or not sources: raise ValueError('Empty or duplicate snapshot')
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    if set(pools) != set(refs) or set(pools) != set(reader.archive_partitions) or not set(sources) <= set(pools):
        raise ValueError('Source scope mismatch')
    if full and set(sources) != set(pools): raise ValueError('Full audit missing sources')
    params = json.loads(Path(cfg['original_model_spec']).read_text())
    reports, replay_candidates = [], []
    for ds in sorted(sources):
        item = sources[ds]; meta, local = read_source(producer, item)
        pm, pool = read_source(inputs, pools[ds]); rm, evidence = read_source(reference, refs[ds])
        rows = reader.rows(ds); n = len(rows)
        for arrays in [local, pool, evidence]: np.testing.assert_array_equal(rows, arrays['native_row'])
        if any(m['native_rows_sha256'] != row_axis_sha256(rows) or m['native_ids_sha256'] != text_axis_sha256(reader.native_ids(ds)) for m in [meta, pm, rm]):
            raise ValueError('Native source axis digest differs')
        if meta['source'] != ds or meta['provisional_family_code_dictionary'] != list(FAMILIES):
            raise ValueError('Source or provisional family dictionary differs')
        for name in ['fit_pool_candidate', 'provisional_family_code']:
            np.testing.assert_array_equal(local[name], pool[name])
        # Rebuild fit eligibility independently from stored component flags.
        expected = (pool['positive_native_RNA'] & pool['in_frozen_context'] & pool['declared_expression_unit_available'] &
                    ~pool['carried_doublet_positive'] & ~pool['carried_primary_exclusion_positive'])
        np.testing.assert_array_equal(expected, pool['fit_pool_candidate'])
        np.testing.assert_array_equal(pool['positive_native_RNA'], evidence['positive_RNA'])
        np.testing.assert_array_equal(pool['in_frozen_context'], reader.membership(ds)['in_frozen_context'])
        validate_partition(rows, local['source_local_group_code'], local['fine_graph_fitted'], expected,
                           local['provisional_family_code'], meta['n_source_local_groups'])
        if meta['n_archive_rows'] != n or meta['n_fit_pool_rows'] != int(expected.sum()) or meta['n_fine_graph_rows'] != int(local['fine_graph_fitted'].sum()):
            raise ValueError('Source denominator metadata differs')
        feature_path = inputs / 'outputs' / pm['feature_candidates']
        if sha256_file(feature_path) != pm['feature_candidates_sha256']: raise ValueError('Feature candidates changed')
        features = json.loads(feature_path.read_text()); genes = features['genes']; columns = np.asarray(features['native_columns'], np.int64)
        lookup = reader.feature_lookup(ds)
        if any(lookup.get(g, set()) != {int(c)} for g, c in zip(genes, columns)): raise ValueError('Native measured-feature axes changed')
        code_by_gene = dict(zip(genes, columns)); observed_models = 0; seen_families = []
        sampled = set(np.linspace(0, n - 1, min(n, 32), dtype=np.int64).tolist())
        for branch in meta['branches']:
            family = branch['provisional_family_code']; seen_families.append(family)
            branch_mask = local['provisional_family_code'] == family
            member = np.flatnonzero(branch_mask & expected); sampled.update(member[:2].tolist())
            if branch['n_archive_rows'] != int(branch_mask.sum()) or branch['n_fit_pool_rows'] != len(member):
                raise ValueError('Branch denominators differ')
            assigned = np.unique(local['source_local_group_code'][member]); status = branch['status']
            if status == 'full_qualified_source_branch_graph':
                if len(member) < params['minimum_fit_group'] or not local['fine_graph_fitted'][member].all(): raise ValueError('Fine branch fit scope differs')
                for field in ['model', 'model_metadata']:
                    if sha256_file(producer / branch[field]) != branch[field + '_sha256']: raise ValueError('Model digest changed')
                mm = json.loads((producer / branch['model_metadata']).read_text())
                with np.load(producer / branch['model'], allow_pickle=False) as f: model = {k: f[k] for k in f.files}
                chosen = np.asarray([code_by_gene[g] for g in mm['native_feature_names']], np.int64)
                if len(set(chosen)) != len(chosen) or set(chosen) & set(features['heldout_physical_alias_columns_excluded']):
                    raise ValueError('Model feature alias or RNA readout leak')
                if mm['source'] != ds or mm['branch_id'] != branch['branch_id'] or mm['fit_axis_sha256'] != row_axis_sha256(rows[member]):
                    raise ValueError('Model identity or row hash differs')
                if mm['source_local_group_offset'] != int(assigned.min()) or len(assigned) != branch['n_groups']:
                    raise ValueError('Model group boundaries differ')
                validate_model_axes(model, member, rows[member], local['source_local_group_code'][member], chosen)
                observed_models += 1
                if len(member) <= spec['maximum_replay_branch_rows']:
                    replay_candidates.append((len(member), ds, branch['branch_id'], branch))
                del model
            elif status in ['small_branch_coarse_only', 'too_few_variable_features_coarse_only']:
                if len(assigned) != 1 or local['fine_graph_fitted'][member].any(): raise ValueError('Coarse-only branch status differs')
                if status == 'small_branch_coarse_only' and not 0 < len(member) < params['minimum_fit_group']:
                    raise ValueError('Small branch count differs')
                if status == 'too_few_variable_features_coarse_only' and branch['n_native_variable_candidates'] >= 3:
                    raise ValueError('Feature-limited coarse branch differs')
            elif status == 'no_qualified_fit_rows_archive_kept':
                if len(member): raise ValueError('Empty-fit status changed')
            else: raise ValueError('Unrecognized branch status')
        if sorted(seen_families) != sorted(np.unique(local['provisional_family_code']).tolist()):
            raise ValueError('Source branch list incomplete or duplicated')
        sampled = np.array(sorted(sampled), np.int64)
        native = reader.store(ds); key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        raw = direct_native_rows(native[key], rows[sampled])
        np.testing.assert_array_equal(np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel(), evidence['native_total_RNA'][sampled])
        reports.append({'source': ds, 'n_archive_rows': n, 'n_models_checked': observed_models, 'n_native_rows_directly_recomputed': len(sampled),
                        'source_arrays_sha256': item['arrays_sha256'], 'source_metadata_sha256': item['metadata_sha256']})
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                  'source': ds, 'completed_sources': len(reports)})
        reader._fields.clear(); del local, pool, evidence, raw; gc.collect()
    if not replay_candidates: raise ValueError('No bounded full branch available for independent replay')
    _, ds, branch_id, branch = sorted(replay_candidates, key=lambda r: r[:3])[0]
    pm, pool = read_source(inputs, pools[ds]); rm, evidence = read_source(reference, refs[ds])
    features = json.loads((inputs / 'outputs' / pm['feature_candidates']).read_text())
    with np.load(producer / branch['model'], allow_pickle=False) as f: model = {k: f[k] for k in f.files}
    member = model['source_row_positions']
    replay = replay_branch(reader, ds, reader.rows(ds), member, features['genes'], np.asarray(features['native_columns'], np.int64),
                           evidence['native_total_RNA'][member], model, params, cfg)
    replay['branch_id'] = branch_id
    output = root / 'outputs/source_local_graph_audit_v1.json'
    save(output, {'scope': spec['scope'], 'sources': reports, 'bounded_independent_branch_replay': replay,
                  'all_graphs_recomputed': False, 'new_final_identity': False})
    terminal = {'status': 'FULL_LOCAL_GRAPH_AUDIT_COMPLETE' if full else 'COMPLETED_SOURCE_GRAPH_PREFLIGHT_ONLY',
                'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'scope': spec['scope'],
                'n_sources': len(reports), 'n_archive_rows': sum(r['n_archive_rows'] for r in reports),
                'n_models_checked': sum(r['n_models_checked'] for r in reports),
                'n_native_rows_directly_recomputed': sum(r['n_native_rows_directly_recomputed'] for r in reports),
                'replayed_branch_rows': len(member), 'output_sha256': sha256_file(output),
                'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                'not_identity_truth': True}
    if full: terminal['producer_terminal_sha256'] = sha256_file(producer / 'provenance/terminal_v1.json')
    save(root / 'provenance/terminal_v1.json', terminal)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); root = Path(parser.parse_args().output)
    try: main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
