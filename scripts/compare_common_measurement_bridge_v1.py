"""Compare original/common measurement references on the same query cells.

Uses four-way matched finite coverage, frozen native targets and exactly matched
training anchors. Source-level diagnostics remain conditional, not patient tests.
"""
import argparse
import datetime
import json
from pathlib import Path
import time

import numpy as np

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from lineage_holdout_diagnostics_v1 import VARIANTS
from run_annotation_source_holdout_v1 import conditional_scores


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/comparison_spec_v1.json').read_text())
    for entry in cfg['bindings'].values():
        if sha256_file(entry['path']) != entry['sha256']:
            raise ValueError('Frozen comparison input changed')
    runs = {name: Path(cfg[name]) for name in ['original', 'common']}
    for path in [*runs.values(), Path(cfg['original_audit']), Path(cfg['common_audit'])]:
        terminal = json.loads((path / 'provenance/supervisor_terminal_v1.json').read_text())
        if terminal['actual_child_exit_code'] != 0 or terminal['wall_timeout']:
            raise ValueError('Producer and independent audit completion required')
    indices = {name: json.loads((path / 'outputs/lineage_source_holdout_index_v1.json').read_text()) for name, path in runs.items()}
    by_key = {name: {(r['fold'], r['route']): r for r in records} for name, records in indices.items()}
    if set(by_key['original']) != set(by_key['common']):
        raise ValueError('Route/fold coverage differs')
    definitions = {name: json.loads((path / 'outputs/readout_definition_v1.json').read_text()) for name, path in runs.items()}
    if definitions['original'] != definitions['common']:
        raise ValueError('Native readout definition differs')
    programs = definitions['original']['programs']; profiles = []; missing = []
    for key in sorted(by_key['original']):
        entries = {name: table[key] for name, table in by_key.items()}
        if any(entry['status'] != 'complete' for entry in entries.values()):
            missing.append({'fold': key[0], 'route': key[1], 'status': 'unavailable_in_at_least_one_feature_universe'})
            continue
        folders = {name: runs[name] / entry['relative_folder'] for name, entry in entries.items()}
        axes, predictions, support, scopes = {}, {}, {}, {}
        for name, folder in folders.items():
            for filename, digest in entries[name]['files_SHA256'].items():
                if sha256_file(folder / filename) != digest:
                    raise ValueError('Frozen bridge output changed')
            with np.load(folder / 'axes_and_readouts.npz', allow_pickle=False) as f:
                axes[name] = {k: f[k] for k in f.files}
            coverage = json.loads((folder / 'models_and_coverage.json').read_text())
            scopes[name] = [r['scope'] for r in coverage['scopes']]
            for variant in VARIANTS:
                with np.load(folder / (variant + '_coordinates.npz'), allow_pickle=False) as f:
                    support[name + '__' + variant] = f['query_supported']
        if scopes['original'] != scopes['common'] or set(axes['original']) != set(axes['common']):
            raise ValueError('Identity scope or stored axes differ')
        for name, values in axes['original'].items():
            np.testing.assert_array_equal(values, axes['common'][name])
        axis = axes['original']; source = axis['query_source']; direct = axis['query_direct_candidate']
        target = axis['query_RNA_readouts']; totals = axis['query_native_totals']
        for scope in scopes['original']:
            for name, folder in folders.items():
                with np.load(folder / (scope + '_predictions.npz'), allow_pickle=False) as f:
                    for variant in VARIANTS:
                        predictions[name + '__' + variant] = f[variant]
            wanted = np.ones(len(source), bool) if scope == 'all_route_cells_not_pure_type' else direct == scope
            for ds in np.unique(source[wanted]):
                subset = wanted & (source == ds)
                for j, program in enumerate(programs):
                    matched = subset & np.isfinite(target[:, j]) & np.isfinite(totals)
                    for condition, pred in predictions.items():
                        matched &= support[condition] & np.isfinite(pred[:, j])
                    metrics = {condition: conditional_scores(target[matched, j], pred[matched, j], totals[matched], cfg['minimum_query_evaluation'])
                               for condition, pred in predictions.items()}
                    differences = {}
                    for variant in VARIANTS:
                        first, second = metrics['original__' + variant], metrics['common__' + variant]
                        if first.get('depth_adjusted_structure_status') == second.get('depth_adjusted_structure_status') == 'estimated':
                            differences[variant] = second['depth_adjusted_Pearson'] - first['depth_adjusted_Pearson']
                    profiles.append({'fold': key[0], 'route': key[1], 'source': str(ds), 'identity_scope': scope,
                        'program': program, 'n_query_scope_cells': int(subset.sum()), 'n_four_way_matched_queries': int(matched.sum()),
                        'four_way_matched_fraction': float(matched.sum() / subset.sum()), 'metrics': metrics,
                        'delta_common_minus_original_depth_adjusted_Pearson': differences})
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'), {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'fold': key[0], 'route': key[1], 'n_profiles': len(profiles)})
    path = root / 'outputs/common_measurement_paired_profiles_v1.json'
    save(path, {'profiles': profiles, 'unavailable_route_folds': missing,
        'anchors_and_native_readouts_exactly_equal': True, 'all_four_conditions_share_query_denominator': True,
        'new_identity_accuracy_or_independent_patient_claim': False, 'production_winner_selected': False})
    save(root / 'provenance/terminal_v1.json', {'status': 'COMMON_MEASUREMENT_BRIDGE_PAIRED_SOURCE_COMPARISON_COMPLETE',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'output_sha256': sha256_file(path),
        'n_profiles': len(profiles), 'n_unavailable_route_folds': len(missing), 'wall_seconds': time.monotonic() - start,
        'native_counts_read_or_copied': False, 'new_authoritative_labels': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
