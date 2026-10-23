"""Audit same-input discovery partitions without treating old labels as truth."""
import argparse
import collections
import datetime
import gzip
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import adjusted_rand_score

from full_context_atlas_v1 import FullContextAtlas, sha256_file


def load(path):
    return json.loads(path.read_text())


def reconcile_partition(labels, source, route, oldfine, tables, method, modeled, minimum):
    if labels.shape != modeled.shape or labels.dtype.kind not in 'iu' or (labels < -1).any():
        raise ValueError('Invalid candidate partition axis/domain')
    expected = modeled.copy()
    grouping = source if method == 'A' else route
    if method != 'B':
        for group in np.unique(grouping):
            mask = modeled & (grouping == group)
            if mask.sum() < minimum:expected[mask] = False
    if not np.array_equal(labels >= 0, expected):
        raise ValueError('Small/zero-RNA/not-modeled cells silently assigned or dropped')
    if len(tables) != len(set(labels)):
        raise ValueError('Partition table size differs')
    for row in tables:
        code = row['partition_code']
        member = labels == code
        if int(member.sum()) != row['n_cells']:
            raise ValueError('Partition denominator differs')
        for observed, values in [('source_counts', source), ('inherited_route_counts', route),
                                 ('inherited_fine_counts_not_truth', oldfine)]:
            if dict(collections.Counter(values[member])) != row[observed]:
                raise ValueError('Partition source/route/old-fine crosswalk differs')
        if code >= 0 and method != 'B' and len(set(grouping[member])) != 1:
            raise ValueError('Source-local or routed partition crosses declared boundary')


def main(args):
    output, run, pilot = Path(args.output), Path(args.run), Path(args.pilot)
    supervisor = load(run / 'provenance/supervisor_terminal_v1.json')
    if supervisor['actual_child_exit_code'] != 0 or supervisor['wall_timeout']:
        raise ValueError('Discovery process did not complete successfully')
    terminal = load(run / 'provenance/terminal_v1.json')
    if terminal['status'] != 'CONTROLLED_DISCOVERY_PARTITIONS_COMPLETE_NOT_FINAL_ANNOTATIONS_OR_HOLDOUT':
        raise ValueError('Unexpected discovery completion status')
    request = load(run / 'config/supervisor_request_v1.json')
    for name, digest in request['script_SHA256'].items():
        if sha256_file(run / name) != digest:raise ValueError('Executed code changed')
    cfg = load(run / 'config/model_spec_v1.json')
    for name, entry in cfg['bindings'].items():
        if sha256_file(entry['path']) != entry['sha256']:raise ValueError('Frozen model input changed: ' + name)
    reader = FullContextAtlas(args.manifest)
    inventory = load(pilot / 'outputs/pilot_input_inventory_v1.json')
    expected_source, expected_route, expected_fine = [], [], []
    # Independent literal mapping audit. Unresolved and unsupported author labels
    # stay unresolved; no speculative T/NK or malignancy promotion.
    aliases = {
        'T_NK': ['T cell', 'T/NK lymphoid', 'T', 'NK', 'T_NK_boundary'],
        'B': ['B cell', 'B'], 'Plasma': ['Plasma cell', 'Plasma'],
        'Myeloid': ['Myeloid', 'Monocyte', 'Macrophage', 'DC', 'Dendritic cell'],
        'Neutrophil': ['Neutrophil'], 'Erythroid': ['Erythroid'], 'Mast': ['Mast cell', 'Mast'],
        'Epithelial': ['Epithelial cell', 'Epithelial'], 'Endothelial': ['Endothelial cell', 'Endothelial'],
        'Fibroblast': ['Fibroblast', 'Fibroblast cell'],
        'Mural': ['Pericyte', 'Smooth muscle cell', 'Mural', 'Mural cell', 'Pericyte_SMC'],
        'Glial': ['Oligodendrocyte']}
    mapping = {name: parent for parent, names in aliases.items() for name in names}
    for item in inventory:
        ds = item['dataset'];path = pilot / item['selector_file']
        if sha256_file(path) != item['selector_sha256']:raise ValueError('Pilot selector changed')
        positions = np.load(path, allow_pickle=False)
        frozen = reader.frozen_labels(ds)
        expected_source.extend([ds] * len(positions))
        expected_route.extend([mapping.get(str(v), 'unresolved_context') for v in frozen['major'][positions]])
        expected_fine.extend(frozen['fine'][positions].astype(str))
    source, route, fine = map(np.asarray, [expected_source, expected_route, expected_fine])
    with np.load(run / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        modeled = f['modeled'].copy()
        if not np.array_equal(source, f['source']) or not np.array_equal(route, f['inherited_route']):
            raise ValueError('Derived discovery metadata axes do not match exact input order')
    if modeled.shape != source.shape or modeled.dtype != bool:
        raise ValueError('Invalid model eligibility mask')
    # Original native stage checked complete RNA >0 for every frozen archive row.
    native_terminal = load(Path(args.native_evidence) / 'provenance/terminal_v1.json')
    if native_terminal['n_native_RNA_evaluable'] == native_terminal['n_archive_cells'] and not modeled.all():
        raise ValueError('Pilot zero-RNA mask conflicts with audited full native archive')
    plan = load(pilot / 'config/model_features_and_heldouts_v1.json')
    selected = load(run / 'outputs/discovery_feature_selection_v1.json')
    indices = np.asarray(selected['selected_feature_indices'])
    genes = selected['selected_feature_names']
    candidates = plan['common_model_feature_candidates']
    if (indices.dtype.kind not in 'iu' or len(indices) != cfg['n_features'] or
            len(set(indices)) != len(indices) or indices.min() < 0 or indices.max() >= len(candidates) or
            genes != [candidates[j] for j in indices]):
        raise ValueError('Frozen discovery feature axis differs')
    heldout = {g for panels in plan['heldout_RNA_readouts'].values() for panel in panels.values() for g in panel}
    for ds in set(source):
        lookup = reader.feature_lookup(ds)
        physical = {next(iter(lookup[g])) for g in genes if len(lookup.get(g, set())) == 1}
        held = {next(iter(lookup[g])) for g in heldout if len(lookup.get(g, set())) == 1}
        if len(physical) != len(genes) or physical & held:raise ValueError('Physical feature/held-out leakage')
    partitions, scientific_summaries = {}, []
    for seed in cfg['model_seeds']:
        with np.load(run / 'outputs' / ('candidate_partitions_' + str(seed) + '.npz'), allow_pickle=False) as f:
            labels = {method: f[method].copy() for method in ['A', 'B', 'C']}
        tables = load(run / 'outputs' / ('partition_crosswalk_diagnostic_' + str(seed) + '.json'))
        for method, values in labels.items():
            reconcile_partition(values, source, route, fine, tables[method], method, modeled, cfg['minimum_fit_group'])
            for ds in np.unique(source):
                mask = (source == ds) & (values >= 0) & modeled
                path = run / 'outputs' / ('native_partition_panels_' + str(seed) + '_' + method + '_' + ds + '.json.gz')
                with gzip.open(path, 'rt') as f:panels = json.load(f)
                if {r['partition_code']: r['n_cells'] for r in panels} != dict(collections.Counter(values[mask])):
                    raise ValueError('Native panel group denominators differ')
                for row in panels:
                    for panel in row['lineage_panel_evidence'].values():
                        n = panel['n_measured']
                        mean, detected = panel['mean_native_logCP10K'], panel['mean_detection_fraction']
                        if n == 0:
                            if mean is not None or detected is not None:raise ValueError('Unmeasured panel represented as zero')
                        elif not (0 <= n <= panel['n_declared'] and mean is not None and np.isfinite(mean) and
                                  0 <= mean <= np.log1p(10000) + 1e-5 and detected is not None and 0 <= detected <= 1):
                            raise ValueError('Native panel evidence out of bounds')
        diagnostics = load(run / 'outputs' / ('representation_diagnostics_' + str(seed) + '.json'))
        rows = diagnostics['conditional_RNA_diagnostics']
        axes = [(r['source'], r['inherited_route'], r['representation']) for r in rows]
        if len(axes) != len(set(axes)):raise ValueError('Repeated RNA diagnostic evidence')
        for row in rows:
            if row['status'] == 'representation_unavailable_small_group':continue
            member = modeled & (source == row['source']) & (route == row['inherited_route'])
            if row.get('n_cells') != min(int(member.sum()), cfg['diagnostic_cells_per_stratum']):
                raise ValueError('RNA diagnostic denominator changed')
            for name, value in row['states'].items():
                if value['status'] == 'estimated':
                    for key in ['depth_adjusted_neighbor_Pearson', 'depth_adjusted_neighbor_Spearman']:
                        if not np.isfinite(value[key]) or not -1 <= value[key] <= 1:
                            raise ValueError('Invalid RNA diagnostic correlation')
                    if row['inherited_route'] == 'T_NK':
                        scientific_summaries.append({'seed': seed, 'source': row['source'], 'state': name,
                                                    'representation': row['representation'], **value})
        partitions[seed] = labels
    stability = []
    seeds = cfg['model_seeds']
    for method in ['A', 'B', 'C']:
        for i, seed1 in enumerate(seeds):
            for seed2 in seeds[i + 1:]:
                a, b = partitions[seed1][method], partitions[seed2][method]
                included = (a >= 0) & (b >= 0)
                for ds in np.unique(source):
                    mask = included & (source == ds)
                    if mask.sum():
                        stability.append({'method': method, 'seed1': seed1, 'seed2': seed2, 'source': str(ds),
                                          'n_cells': int(mask.sum()),
                                          'partition_seed_ARI_not_identity_accuracy': float(adjusted_rand_score(a[mask], b[mask]))})
    output_bytes = sum(p.stat().st_size for p in run.rglob('*') if p.is_file())
    if output_bytes > cfg['output_cap_bytes']:raise ValueError('Bounded discovery output cap exceeded')
    manifest = [{'file': str(p.relative_to(run)), 'sha256': sha256_file(p), 'bytes': p.stat().st_size}
                for p in sorted((run / 'outputs').rglob('*')) if p.is_file()]
    (output / 'outputs').mkdir(exist_ok=True)
    for name, value in [('audited_output_file_manifest_v1.json', manifest),
                        ('private_T_NK_RNA_diagnostics_v1.json', scientific_summaries),
                        ('private_partition_seed_stability_v1.json', stability)]:
        with (output / 'outputs' / name).open('x') as f:json.dump(value, f, indent=2, allow_nan=False);f.write('\n')
    receipt = {'status': 'DISCOVERY_PARTITION_AXES_AND_BOUNDS_AUDIT_PASS_NOT_FINAL_ABC_IDENTITY',
               'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'all_seeds_and_method_partitions_checked': True,
               'all_crosswalk_native_panel_denominators_and_RNA_diagnostic_bounds_pass': True,
               'physical_feature_holdout_disjointness_pass': True,
               'frozen_pilot_order_masks_and_small_group_handling_pass': True,
               'independent_native_panel_value_replay': 'pending_separate_replay_not_claimed',
               'RNA_neighborhood_numeric_replay': 'pending_separate_replay_not_claimed',
               'source_holdout_or_identity_accuracy_or_method_winner_accepted': False,
               'output_manifest_SHA256': sha256_file(output / 'outputs/audited_output_file_manifest_v1.json'),
               'original_science_actual_exit_code': supervisor['actual_child_exit_code'],
               'original_terminal_SHA256': sha256_file(run / 'provenance/terminal_v1.json'),
               'n_scientific_seeds': len(seeds), 'n_cells': len(source), 'output_bytes': output_bytes}
    with (output / 'provenance/terminal_v1.json').open('x') as f:json.dump(receipt, f, indent=2);f.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['output', 'run', 'pilot', 'manifest', 'native-evidence']:parser.add_argument('--' + name, required=True)
    main(parser.parse_args())
