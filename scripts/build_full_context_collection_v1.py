"""Create and audit a bounded full-context index; never copy expression arrays."""
import argparse
import collections
import datetime
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

import numpy as np
from backed_atlas_collection_v1 import read_csr_rows
from full_context_atlas_v1 import FullContextAtlas, sha256_file


MISSING = {'', 'nan', 'None', '__SOURCE_VALUE_MISSING__', '__MISSING__', 'NA'}
PANEL = ['CD3D', 'CD3E', 'TRAC', 'TRDC', 'TRGC1', 'NKG7', 'GNLY', 'MS4A1',
         'CD79A', 'LST1', 'LYZ', 'PECAM1', 'VWF', 'EPCAM', 'KRT8', 'COL1A1',
         'DCN', 'FOXP3', 'IL2RA', 'GZMB', 'IFIT1', 'MKI67', '__MISSING_FEATURE_FIXTURE__']


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write('\n')


def check_size(root, limit):
    actual = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if actual > limit:
        raise RuntimeError('Metadata byte cap exceeded')
    return actual


def main(args):
    output = Path(args.output)
    if (output / 'provenance/terminal_v1.json').exists():
        raise ValueError('Preserve completed attempt; choose a new output')
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    signal.alarm(args.wall_seconds)
    test = subprocess.run([sys.executable, '-B', str(output / 'tests/test_full_context_atlas_v1.py'), '-v'],
                          capture_output=True, text=True, timeout=120)
    save(output / 'provenance/synthetic_tests_v1.json',
         {'actual_exit_code': test.returncode, 'stdout': test.stdout, 'stderr': test.stderr})
    if test.returncode != 0 or 'Ran 6 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise RuntimeError('Actual scientific numerical fixtures did not pass')
    projection_path = Path(args.projection)
    projection = json.loads(projection_path.read_text())
    primary_path = Path(args.primary_manifest)
    descriptor_path = Path(args.legacy_descriptor)
    descriptor = json.loads(descriptor_path.read_text())
    for path, expected in [(primary_path, projection['primary_manifest_sha256']),
                           (Path(args.legacy_sample_metadata), projection['legacy_sample_metadata_sha256'])]:
        if sha256_file(path) != expected:
            raise ValueError('Changed existing frozen scope input')
    bindings = []
    for role, path in [('full scope', projection_path), ('primary collection', primary_path),
                       ('legacy descriptor', descriptor_path), ('legacy roster', Path(args.legacy_roster)),
                       ('legacy sample scope', Path(args.legacy_sample_metadata))]:
        bindings.append({'role': role, 'path': str(path), 'sha256': sha256_file(path)})
    manifest = {'schema': 'full_context_atlas_collection_v1', 'created_UTC': now(),
                'status': 'readable_archive_index_not_new_annotation_or_joint_embedding',
                'primary_manifest': str(primary_path), 'bindings': bindings,
                'legacy': {'canonical_path': descriptor['canonical_path'], 'counts_key': 'X',
                           'native_source_root': descriptor['native_source_root'],
                           'roster_path': args.legacy_roster,
                           'sample_metadata_path': args.legacy_sample_metadata},
                'partitions': projection['parts'],
                'expected_context_cells': projection['expected_context_cells'],
                'expected_context_expression_units': projection['expected_context_expression_units'],
                'cross_study_alias_dedup_closed': False,
                'independent_patient_count_inferred': False,
                'source_local_labels_are_carried_not_biological_validation': True,
                'large_matrix_copied': False}
    manifest_path = output / 'config/full_context_collection_v1.json'
    save(manifest_path, manifest)
    reader = FullContextAtlas(manifest_path)
    primary = reader.primary.manifest
    profiles = []
    for number, (dataset, partition) in enumerate(reader.partitions.items(), 1):
        rows = reader.rows(dataset)
        labels = reader.frozen_labels(dataset)
        majors = collections.Counter(str(x) for x in labels['major'])
        if sum(majors.values()) != len(rows):
            raise ValueError('Label/source axis mismatch')
        unit_field = partition.get('unit_field', 'sampleID')
        units = reader.field(dataset, unit_field)
        valid = np.asarray([str(v) not in MISSING for v in units])
        if not valid.all() or len(set(units)) != partition['n_context_units']:
            raise ValueError('Original expression-unit scope mismatch')
        lookup = reader.feature_lookup(dataset)
        positions = [len(rows) - 1, 0, len(rows) // 2, 0]
        observed = reader.read_panel(dataset, positions, PANEL)
        z = reader.store(dataset)
        key = reader.legacy_config['counts_key'] if reader.is_legacy(dataset) else reader.primary.partitions[dataset]['counts_key']
        actual_matrix = z[key]
        expected = np.zeros(observed['counts'].shape, dtype=np.float64)
        independent_totals = []
        sampled_raw = []
        for row_index, native_row in enumerate(rows[positions]):
            lo, hi = map(int, actual_matrix['indptr'][int(native_row):int(native_row) + 2])
            indices = np.asarray(actual_matrix['indices'][lo:hi])
            values = np.asarray(actual_matrix['data'][lo:hi])
            sampled_raw.append(values)
            independent_totals.append(sum(float(v) for v in values))
            for column_index, gene in enumerate(PANEL):
                if len(lookup.get(gene, set())) == 1:
                    physical = next(iter(lookup[gene]))
                    expected[row_index, column_index] = sum(float(v) for i, v in zip(indices, values) if i == physical)
        if not np.array_equal(observed['counts'].toarray(), expected):
            raise ValueError('Independent native CSR panel replay failed')
        if not np.array_equal(observed['native_total_counts'], independent_totals):
            raise ValueError('Native denominator changed')
        if observed['measured'][-1] or observed['counts'][:, -1].nnz:
            raise ValueError('Missing-feature contract failed')
        data = np.concatenate(sampled_raw)
        sampled_count_eligible = bool(np.isfinite(data).all() and (data >= 0).all() and np.equal(data, np.floor(data)).all())
        fields = reader.field_names(dataset)
        profile = {'dataset': dataset, 'provider': partition['provider'],
                   'n_context_cells': len(rows), 'unit_field': unit_field,
                   'n_expression_units': len(set(units)),
                   'frozen_major_counts': dict(majors),
                   'n_frozen_fine_labels': len(set(str(v) for v in labels['fine'])),
                   'carried_annotation_status': labels['status'],
                   'n_unique_resolvable_features': len({next(iter(v)) for v in lookup.values() if len(v) == 1}),
                   'diagnostic_panel_measured': {gene: bool(flag) for gene, flag in zip(PANEL, observed['measured'])},
                   'metadata_fields': fields,
                   'sampled_raw_count_domain_eligible': sampled_count_eligible,
                   'sampled_count_domain_not_full_source_certificate': True,
                   'independent_panel_and_denominator_replay_pass': True,
                   'source_cell_order_and_feature_membership_checked': True,
                   'cross_study_alias_status': 'pending_paper_and_native_library_review',
                   'new_biological_identity_promoted': False,
                   'elapsed_seconds': time.monotonic() - started}
        if not reader.is_legacy(dataset):
            prior = reader.primary.partitions[dataset]
            profile['source_design'] = {k: prior[k] for k in ['cancer_type', 'capture_scope', 'therapy_background', 'role']}
        else:
            profile['source_design'] = {'status': 'frozen_source_metadata_and_external_registry_bridge_pending'}
        profiles.append(profile)
        save(output / ('checkpoints/source_%03d.json' % number),
             {'observed_UTC': now(), 'n_sources_audited': number, 'last_source': dataset,
              'completed_source_cells': sum(p['n_context_cells'] for p in profiles),
              'last_profile': profile})
        check_size(output, args.metadata_cap_bytes)
    if sum(p['n_context_cells'] for p in profiles) != manifest['expected_context_cells']:
        raise ValueError('Full-context cell count reconciliation failed')
    if sum(p['n_expression_units'] for p in profiles) != manifest['expected_context_expression_units']:
        raise ValueError('Full-context expression-unit reconciliation failed')
    save(output / 'outputs/source_structure_and_annotation_coverage_v1.json', profiles)
    terminal = {'status': 'FULL_CONTEXT_NATIVE_COLLECTION_STRUCTURE_AND_CARRIED_LABELS_AUDITED',
                'completed_UTC': now(), 'actual_exit_code': 0,
                'n_source_keys': len(profiles),
                'n_context_cells': sum(p['n_context_cells'] for p in profiles),
                'n_expression_units_not_patients': sum(p['n_expression_units'] for p in profiles),
                'all_sources_independent_CSR_panel_replay': all(p['independent_panel_and_denominator_replay_pass'] for p in profiles),
                'all_sampled_raw_count_domains_eligible': all(p['sampled_raw_count_domain_eligible'] for p in profiles),
                'metadata_bytes_before_terminal': check_size(output, args.metadata_cap_bytes),
                'max_RSS_bytes_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                'wall_seconds': time.monotonic() - started,
                'cross_study_alias_dedup_closed': False,
                'new_coarse_annotation_or_A_B_C_fitted': False,
                'single_materialized_Zarr_created': False,
                'large_matrix_copied_or_downloaded': False,
                'biological_identity_or_clinical_release_accepted': False,
                'manifest_sha256': sha256_file(manifest_path),
                'coverage_sha256': sha256_file(output / 'outputs/source_structure_and_annotation_coverage_v1.json')}
    save(output / 'provenance/terminal_v1.json', terminal)
    print(json.dumps(terminal))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for flag in ['projection', 'primary-manifest', 'legacy-descriptor', 'legacy-roster',
                 'legacy-sample-metadata', 'output']:
        parser.add_argument('--' + flag, required=True)
    parser.add_argument('--metadata-cap-bytes', type=int, default=16 * 1024 ** 2)
    parser.add_argument('--wall-seconds', type=int, default=5400)
    args = parser.parse_args()
    try:
        main(args)
    except Exception as error:
        failure = Path(args.output) / 'provenance/failure_v1.json'
        if not failure.exists():
            save(failure, {'observed_UTC': now(), 'status': 'FAILED_NOT_ACCEPTED',
                           'error_class': type(error).__name__, 'error': str(error)})
        raise
