"""Stream full native marker evidence by inherited group without relabeling."""
import argparse
import datetime
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import (GENES, LINEAGE_PANELS, STATE_PANELS,
                                       IDENTITY_CANDIDATE_PANELS, MarkerEvidence,
                                       native_marker_block)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, ensure_ascii=False, separators=(',', ':'), allow_nan=False)
        handle.write('\n')


def check_size(root, limit):
    size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size > limit:
        raise RuntimeError('Bounded evidence directory byte cap exceeded')
    return size


def compact_evidence(e):
    out = []
    for i, group in enumerate(e.group_names):
        n = int(e.n_evaluable[i])
        row = {'inherited_group': group, 'n_archive_cells': int(e.n_cells[i]),
               'n_native_RNA_evaluable': n, 'n_detected': [],
               'mean_native_CP10K': [], 'mean_native_log1p_CP10K': []}
        for j in range(len(e.genes)):
            valid = bool(e.measured[j]) and n > 0
            row['n_detected'].append(int(e.detected[i, j]) if valid else None)
            row['mean_native_CP10K'].append(float(e.cp10k_sum[i, j] / n) if valid else None)
            row['mean_native_log1p_CP10K'].append(float(e.logcp10k_sum[i, j] / n) if valid else None)
        out.append(row)
    return out


def main(args):
    output = Path(args.output)
    started = time.monotonic()
    def timeout(_signal, _frame):
        raise TimeoutError('Frozen stage wall budget exhausted')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(args.wall_seconds)
    stage1 = Path(args.source_run)
    audit = json.loads((stage1 / 'provenance/independent_stage1_audit_v1.json').read_text())
    if audit['status'] != 'INDEPENDENT_STAGE1_STRUCTURE_AND_FROZEN_LABEL_COUNT_AUDIT_PASS':
        raise ValueError('Source stage independent audit not accepted')
    manifest_path = stage1 / 'config/full_context_collection_v1.json'
    terminal = json.loads((stage1 / 'provenance/terminal_v1.json').read_text())
    if sha256_file(manifest_path) != terminal['manifest_sha256']:
        raise ValueError('Changed source collection')
    test = subprocess.run([sys.executable, '-B', str(output / 'tests/test_native_panel_aggregate_v1.py'), '-v'],
                          capture_output=True, text=True, timeout=120)
    save(output / 'provenance/numerical_tests_v1.json',
         {'actual_exit_code': test.returncode, 'stdout': test.stdout, 'stderr': test.stderr})
    if test.returncode != 0 or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise RuntimeError('Native evidence numerical tests not accepted')
    reader = FullContextAtlas(manifest_path)
    save(output / 'config/native_marker_contract_v1.json',
         {'created_UTC': now(), 'source_manifest_sha256': terminal['manifest_sha256'],
          'genes': GENES, 'lineage_panels': LINEAGE_PANELS, 'state_panels': STATE_PANELS,
          'identity_candidate_panels': IDENTITY_CANDIDATE_PANELS,
          'normalization': 'per-cell native counts / complete native RNA total * 10000; then log1p',
          'zero_RNA_policy': 'archive cells retained; excluded from expression means and detection denominators',
          'missing_gene_policy': 'not measured remains null, never biological zero',
          'outputs_are_evidence_not_new_identity_or_state_calls': True,
          'unresolved_author_groups_split_for_review_only': True,
          'new_effect_P_or_FDR': False})
    outputs = []
    progress_number = 0
    total_cells = 0
    for source_number, (dataset, partition) in enumerate(reader.partitions.items(), 1):
        rows = reader.rows(dataset)
        labels = reader.frozen_labels(dataset)
        group_values = np.asarray([str(major) if str(major) != 'unresolved'
                                   else 'unresolved::author::' + str(fine)
                                   for major, fine in zip(labels['major'], labels['fine'])], dtype=object)
        codes, groups = pd.factorize(group_values, sort=True)
        lookup = reader.feature_lookup(dataset)
        measured = np.asarray([len(lookup.get(gene, set())) == 1 for gene in GENES])
        evidence = MarkerEvidence(groups.tolist(), GENES, measured)
        nnz = 0
        for start in range(0, len(rows), args.block_size):
            stop = min(len(rows), start + args.block_size)
            block = native_marker_block(reader, dataset, np.arange(start, stop))
            if not np.array_equal(block['measured'], measured):
                raise ValueError('Marker measurement mask changed within source')
            evidence.update(block['counts'], block['native_total_counts'], codes[start:stop])
            nnz += block['full_native_nnz_read']
            if start == 0 or stop == len(rows) or (start // args.block_size) % 16 == 0:
                progress_number += 1
                save(output / ('checkpoints/progress_%06d.json' % progress_number),
                     {'observed_UTC': now(), 'source_number': source_number, 'dataset': dataset,
                      'source_cells_processed': stop, 'source_cells': len(rows),
                      'full_archive_cells_processed': total_cells + stop,
                      'completed_sources': len(outputs), 'wall_seconds': time.monotonic() - started})
                check_size(output, args.output_cap_bytes)
        if int(evidence.n_cells.sum()) != len(rows):
            raise ValueError('Full native source evidence count differs from archive')
        result = {'dataset': dataset, 'n_archive_cells': len(rows),
                  'n_native_RNA_evaluable': int(evidence.n_evaluable.sum()),
                  'genes': GENES, 'measured': measured.tolist(),
                  'native_nnz_checked': nnz, 'all_selected_native_count_values_checked': True,
                  'group_definition': 'inherited coarse identity; author subgroup only where unresolved',
                  'groups': compact_evidence(evidence),
                  'source_annotation_promoted': False, 'state_or_malignancy_promoted': False}
        destination = output / ('outputs/source_%03d_evidence.json' % source_number)
        save(destination, result)
        outputs.append({'dataset': dataset, 'file': str(destination.relative_to(output)),
                        'sha256': sha256_file(destination), 'n_archive_cells': len(rows),
                        'n_native_RNA_evaluable': int(evidence.n_evaluable.sum()),
                        'n_groups': len(groups), 'n_measured_markers': int(measured.sum()),
                        'native_nnz_checked': nnz})
        total_cells += len(rows)
        check_size(output, args.output_cap_bytes)
    if total_cells != reader.manifest['expected_context_cells']:
        raise ValueError('Full archive evidence reconciliation failed')
    save(output / 'outputs/source_evidence_index_v1.json', outputs)
    save(output / 'provenance/terminal_v1.json',
         {'status': 'FULL_NATIVE_COARSE_MARKER_EVIDENCE_MEASURED_NOT_NEW_ANNOTATION',
          'completed_UTC': now(), 'n_archive_cells': total_cells, 'n_sources': len(outputs),
          'n_native_RNA_evaluable': sum(x['n_native_RNA_evaluable'] for x in outputs),
          'native_nnz_checked': sum(x['native_nnz_checked'] for x in outputs),
          'all_selected_native_count_values_checked': True,
          'new_identity_state_malignancy_or_A_B_C_fitted': False,
          'output_bytes_before_terminal': check_size(output, args.output_cap_bytes),
          'wall_seconds': time.monotonic() - started,
          'max_RSS_bytes_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
          'output_index_sha256': sha256_file(output / 'outputs/source_evidence_index_v1.json')})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-run', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--block-size', type=int, default=2048)
    parser.add_argument('--wall-seconds', type=int, default=10800)
    parser.add_argument('--output-cap-bytes', type=int, default=16 * 1024 ** 2)
    args = parser.parse_args()
    if args.block_size < 1 or args.block_size > 4096:
        raise ValueError('Use 1..4096 selected rows per block')
    try:
        main(args)
    except Exception as error:
        failure = Path(args.output) / 'provenance/failure_v1.json'
        if not failure.exists():
            save(failure, {'observed_UTC': now(), 'status': 'FAILED_NOT_ACCEPTED',
                           'error_class': type(error).__name__, 'error': str(error)})
        raise
