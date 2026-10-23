"""Qualify representative fitting against frozen all-pilot source-local graphs.

The only discovery factor changed is source/route fit population size. Fixed
feature selection is transductive. These are conditional RNA-structure and
mapping-coverage diagnostics, not independent identity accuracy or source holdout.
"""
import argparse
import datetime
import gzip
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from annotation_representation_pilot_v1 import pca_fit, save
from annotation_pilot_sampling_v1 import balanced_cell_selection, stable_seed
from compare_annotation_granularity_v1 import graph_grid
from annotation_granularity_metrics_v1 import granularity, unit_out_state, gamma_group_retention
from paired_granularity_review_v1 import paired_scores, direct_unit_out
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from load_annotation_pilot_native_v1 import load_pilot
from representative_local_mapping_v1 import weighted_vote, without_self, support_mask


def save_npz(path, **arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as f: np.savez_compressed(f, **arrays)
    with np.load(path, allow_pickle=False) as f:
        if set(f.files) != set(arrays): raise ValueError('Serialized field axis changed')
        for key, values in arrays.items(): np.testing.assert_array_equal(f[key], values)


def audit_neighbors(reference, query, labels, d, ix, prediction, confidence):
    """Direct float64 Euclidean distances and dictionary votes; ties are reported."""
    checked = 0
    boundary_ties = 0
    maximum_distance_error = 0.
    for row in np.unique(np.linspace(0, len(query) - 1, min(32, len(query)), dtype=int)):
        distances = np.sqrt(np.sum((reference.astype(float) - query[row].astype(float)) ** 2, axis=1))
        k = ix.shape[1]
        order = np.argsort(distances, kind='stable')
        maximum_distance_error = max(maximum_distance_error, float(np.max(np.abs(distances[ix[row]] - d[row]))))
        np.testing.assert_allclose(distances[ix[row]], d[row], rtol=2e-5, atol=2e-5)
        if len(order) > k and np.isclose(distances[order[k - 1]], distances[order[k]], rtol=2e-5, atol=2e-5):
            boundary_ties += 1
            continue
        chosen = order[:k]
        if set(chosen) != set(ix[row]): raise ValueError('Direct and fitted nearest neighbors differ')
        weights = {}
        zeros = distances[chosen] <= 1e-12
        for j, pos in enumerate(chosen):
            weight = float(zeros[j]) if zeros.any() else 1. / (distances[pos] + 1e-8)
            label = int(labels[pos])
            weights[label] = weights.get(label, 0.) + weight
        winner = min(weights, key=lambda label: (-weights[label], label))
        if winner != prediction[row]: raise ValueError('Direct mapping vote differs')
        np.testing.assert_allclose(confidence[row], weights[winner] / sum(weights.values()), rtol=2e-5, atol=2e-5)
        checked += 1
    return {'n_direct_query_checks': checked, 'n_distance_boundary_ties_not_called_errors': boundary_ties,
            'maximum_direct_distance_error': maximum_distance_error}


def main(root):
    started = time.monotonic()
    cfg = json.loads((root / 'config/representative_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']: raise ValueError('Frozen qualification input changed')
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'),
             {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'stage': stage,
              'wall_seconds': time.monotonic() - started, **extra})
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_representative_local_mapping_v1.py')],
                          capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stdout.txt').write_text(test.stdout)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 6 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Mapping counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 6, 'no_skips': True})
    pilot, discovery, original = (Path(cfg[k]) for k in ['pilot', 'discovery', 'source_lineage'])
    data = load_pilot(pilot, cfg['manifest'], checkpoint=checkpoint)
    selected = json.loads((discovery / 'outputs/discovery_feature_selection_v1.json').read_text())
    features = selected['selected_feature_indices']
    if [data['genes'][i] for i in features] != selected['selected_feature_names']: raise ValueError('Feature axis differs')
    x = data.pop('x')[:, features]
    params = json.loads((discovery / 'config/model_spec_v1.json').read_text())
    with np.load(original / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        source, route, positive = f['source'], f['provisional_native_route'], f['modeled']
    with np.load(original / ('outputs/candidate_partitions_' + str(cfg['seed']) + '.npz'), allow_pickle=False) as f:
        baseline = f['A']
    if not np.array_equal(source, data['source']) or not np.array_equal(positive, data['native_total'] > 0):
        raise ValueError('Qualification cell axes differ')
    if not np.array_equal(baseline >= 0, positive): raise ValueError('Baseline modeling domain differs')
    reader = FullContextAtlas(cfg['manifest'])
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    units = []
    for item in inventory:
        selector = pilot / item['selector_file']
        if sha256_file(selector) != item['selector_sha256']: raise ValueError('Pilot selector changed')
        positions = np.load(selector, allow_pickle=False)
        units.extend(reader.frozen_expression_units(item['dataset'])[positions].astype(str))
    units = np.asarray(units)
    if len(units) != len(source): raise ValueError('Expression-unit axis differs')
    del reader
    proposed = np.full(len(source), -1, np.int32)
    fitted = np.zeros(len(source), bool)
    supported = np.zeros(len(source), bool)
    vote = np.full(len(source), np.nan, np.float32)
    distance = np.full(len(source), np.nan, np.float32)
    records = []
    offset = 0
    for ds in np.unique(source):
        for parent in np.unique(route[source == ds]):
            member = np.flatnonzero((source == ds) & (route == parent) & positive)
            if not len(member): continue
            record = {'source': str(ds), 'route_proposal_not_truth': str(parent), 'n_cells': len(member)}
            branch_id = 'branch_%03d' % len(records)
            if len(member) < params['minimum_fit_group']:
                proposed[member] = offset
                offset += 1
                record.update(status='coarse_only_no_fine_fit', n_fit=0, n_query=0)
            else:
                selected_rows = balanced_cell_selection(units[member], data['fine'][member], np.ones(len(member), bool),
                    cfg['fit_cap_per_source_route'], stable_seed(cfg['seed'], ds, parent), cfg['rare_reserve'])
                train = member[selected_rows]
                query = member[~np.isin(np.arange(len(member)), selected_rows)]
                embedding, model = pca_fit(x[train], cfg['seed'], params['components'], params['pca_iterations'])
                labels = graph_grid(embedding, cfg['seed'], params, [cfg['resolution']])[str(cfg['resolution'])]
                proposed[train] = labels + offset
                fitted[train] = True
                record.update(status='representative_fit_with_query_mapping' if len(query) else 'full_branch_fit_reproduction',
                              n_fit=len(train), n_query=len(query))
                save_npz(root / 'models' / (branch_id + '.npz'), feature_indices=np.asarray(features, np.int32),
                    mean=model.mean_, components=model.components_, explained_variance=model.explained_variance_,
                    fit_pilot_rows=train, fit_embedding=embedding, fit_group_code=labels)
                if not len(query):
                    ari = float(adjusted_rand_score(baseline[member], proposed[member]))
                    if abs(ari - 1.) > 1e-12: raise ValueError('Uncapped branch failed original partition reproduction')
                    record['uncapped_partition_ARI'] = ari
                else:
                    k = min(cfg['mapping_neighbors'], len(train) - 1)
                    nn = NearestNeighbors(n_neighbors=k + 1, algorithm='brute', metric='euclidean', n_jobs=4).fit(embedding)
                    reference_mean_distance = np.empty(len(train), float)
                    for start in range(0, len(train), cfg['query_block_size']):
                        end = min(start + cfg['query_block_size'], len(train))
                        d, ix = nn.kneighbors(embedding[start:end], n_neighbors=k + 1)
                        reference_mean_distance[start:end] = without_self(d, ix, np.arange(start, end), k).mean(axis=1)
                    cutoff = float(np.quantile(reference_mean_distance, cfg['distance_quantile']))
                    audit = []
                    for start in range(0, len(query), cfg['query_block_size']):
                        rows = query[start:start + cfg['query_block_size']]
                        qe = model.transform(x[rows].toarray()).astype(np.float32)
                        d, ix = nn.kneighbors(qe, n_neighbors=k)
                        lab, confidence = weighted_vote(d, labels[ix])
                        proposed[rows] = lab + offset
                        vote[rows], distance[rows] = confidence, d.mean(axis=1)
                        supported[rows] = support_mask(d.mean(axis=1), confidence, cutoff, cfg['minimum_vote'])
                        if start == 0: audit.append(audit_neighbors(embedding, qe, labels, d, ix, lab, confidence))
                    record.update(reference_self_excluded_distance_quantile=cfg['distance_quantile'], distance_cutoff=cutoff,
                                  n_query_supported=int(supported[query].sum()), query_supported_fraction=float(supported[query].mean()),
                                  independent_mapping_checks=audit, query_partition_ARI=float(adjusted_rand_score(baseline[query], proposed[query])))
                    gamma = data['gamma_RNA_candidate']
                    gamma_query = query[np.isfinite(gamma[query]) & (gamma[query] > 0)]
                    record.update(n_gamma_RNA_candidates_fit=int(np.sum(gamma[train] > 0)),
                                  n_gamma_RNA_candidates_query=len(gamma_query),
                                  gamma_query_supported_fraction=float(supported[gamma_query].mean()) if len(gamma_query) else None)
                offset += len(np.unique(labels))
            records.append(record)
            checkpoint('source_local_branch_complete', branch_id=branch_id, **record)
    del x
    if not np.array_equal(proposed >= 0, positive): raise ValueError('Native archive modeling eligibility changed')
    for code in np.unique(proposed[positive]):
        rows = proposed == code
        if len(np.unique(source[rows])) != 1 or len(np.unique(route[rows])) != 1:
            raise ValueError('Proposed local group crosses a source/route boundary')
    save_npz(root / 'outputs/representative_partitions_v1.npz', source=source, route=route, positive_native_RNA=positive,
             baseline_group=baseline, proposed_group=proposed, participated_in_fit=fitted, query_support=supported,
             query_vote_fraction=vote, query_mean_neighbor_distance=distance)
    profiles = []
    n_direct_predictions = 0
    max_direct_prediction_error = 0.
    audit_source = min(inventory, key=lambda item: (item['n_selected_cells'], item['dataset']))['dataset']
    for record in records:
        ds, parent = record['source'], record['route_proposal_not_truth']
        member = np.flatnonzero((source == ds) & (route == parent) & positive)
        if len(member) < params['minimum_fit_group']: continue
        query = ~fitted[member]
        u, depth = units[member], np.log1p(data['native_total'][member])
        row = {**record, 'n_expression_units_not_patients': len(np.unique(u)),
               'baseline_granularity': granularity(baseline[member]), 'representative_granularity': granularity(proposed[member]),
               'baseline_gamma_RNA': gamma_group_retention(baseline[member], data['gamma_RNA_candidate'][member]),
               'representative_gamma_RNA': gamma_group_retention(proposed[member], data['gamma_RNA_candidate'][member]),
               'state_comparisons_all_cells': {}, 'state_comparisons_unfitted_query': {}}
        for state, values in data['state'].items():
            results = []
            for labels in [baseline, proposed]:
                result = unit_out_state(labels[member], u, depth, values[member], cfg['minimum_other_unit_group_cells'], True)
                if ds == audit_source:
                    direct, direct_base = direct_unit_out(labels[member], u, depth, values[member], cfg['minimum_other_unit_group_cells'])
                    for field, expected in [('prediction', direct), ('baseline', direct_base)]:
                        np.testing.assert_allclose(result[field], expected, rtol=1e-8, atol=1e-8, equal_nan=True)
                        valid = np.isfinite(expected)
                        if valid.any(): max_direct_prediction_error = max(max_direct_prediction_error, float(np.max(np.abs(result[field][valid] - expected[valid]))))
                    n_direct_predictions += int(np.isfinite(direct).sum())
                results.append(result)
            row['state_comparisons_all_cells'][state] = paired_scores(values[member], u, *results)
            if query.any():
                row['state_comparisons_unfitted_query'][state] = paired_scores(values[member][query], u[query],
                    *[{k: r[k][query] for k in ['prediction', 'baseline']} for r in results])
        profiles.append(row)
        checkpoint('conditional_RNA_diagnostic_complete', source=ds, route_proposal_not_truth=parent)
    output = root / 'outputs/representative_local_qualification_v1.json.gz'
    with gzip.open(output, 'xt') as f:
        json.dump({'profiles': profiles, 'all_branch_records': records, 'state_gene_sets': data['state_genes'],
                   'feature_selection_transductive': True, 'clinical_or_CM_inputs_used': False,
                   'partition_ARI_not_identity_accuracy': True, 'full_archive_model_acceptance': False,
                   'unsupported_queries_preserved_with_proposal_and_separate_flag': True,
                   'RNA_comparison_training_domain_includes_all_positive_branch_rows': True}, f, allow_nan=False)
    written = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if written > cfg['output_cap_bytes']: raise ValueError('Bounded qualification output cap exceeded')
    save(root / 'provenance/terminal_v1.json', {'status': 'REPRESENTATIVE_LOCAL_QUALIFICATION_COMPLETE_NOT_PRODUCTION_ACCEPTANCE',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_tests_actual_no_skips': 6,
        'n_branch_records': len(records), 'n_profiles': len(profiles), 'n_direct_predictions': n_direct_predictions,
        'maximum_direct_prediction_error': max_direct_prediction_error, 'output_sha256': sha256_file(output),
        'partitions_sha256': sha256_file(root / 'outputs/representative_partitions_v1.npz'),
        'wall_seconds': time.monotonic() - started, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'written_bytes_before_terminal': written, 'native_counts_copied': False, 'new_final_identity': False})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--output', required=True)
    root = Path(p.parse_args().output)
    try: main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
