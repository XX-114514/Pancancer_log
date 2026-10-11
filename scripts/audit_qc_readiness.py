#!/usr/bin/env python3
"""Read QC ledgers and Zarr headers; publish only aggregate readiness evidence.

Does not load expression matrices, modify QC, or grant final source acceptance.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def shape(path):
    attrs = path / '.zattrs'
    if attrs.exists():
        value = json.loads(attrs.read_text())
        if 'shape' in value:
            return value['shape']
    array = path / '.zarray'
    if array.exists():
        return json.loads(array.read_text())['shape']
    raise ValueError('Missing readable Zarr shape metadata')


def audit(run):
    result = {'schema': 'qc_ledger_readiness_v1',
              'observed_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'all_registered_pending_QC_shards', 'shards': [],
              'final_source_acceptance_granted': False,
              'limits': ['Header and ledger audit; numeric matrix readback is the prior producer proof.',
                         'Source aliases, library overlap and clinical linkage remain separate gates.']}
    for proof_path in sorted((run / 'manifests').glob('accepted_*.json')):
        proof = json.loads(proof_path.read_text())
        name = proof['gse_id']
        if proof['status'] != 'qc_shard_written_pending_final_source_acceptance':
            raise ValueError('Unexpected proof status; review schema before publication')
        ledger = run / 'manifests' / (name + '_cell_qc_and_filter_reason.tsv.gz')
        thresholds = run / 'manifests' / (name + '_sample_qc_thresholds.tsv')
        reasons, statuses = Counter(), Counter()
        raw = kept = 0
        with gzip.open(ledger, 'rt', newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                raw += 1
                if row['keep'].lower() in ('true', '1'):
                    kept += 1
                    statuses[row['doublet_status']] += 1
                else:
                    reasons[row['filter_reason'] or 'unspecified'] += 1
        with thresholds.open(newline='') as handle:
            libraries = list(csv.DictReader(handle, delimiter='\t'))
        zarr = Path(proof['output'])
        # Evidence paths must stay within the explicitly selected run.
        zarr.relative_to(run)
        x_shape = shape(zarr / 'X')
        counts_shape = shape(zarr / 'layers' / 'counts')
        checks = {
            'ledger_raw_matches_proof': raw == proof['n_raw'],
            'ledger_keep_matches_proof': kept == proof['n_qc_nonduplicate'],
            'library_raw_sum_matches_ledger': sum(int(r['n_raw']) for r in libraries) == raw,
            'X_and_counts_shape_match': x_shape == counts_shape,
            'Zarr_rows_match_keep': x_shape[0] == kept,
            'every_library_doublet_computed': all(r['doublet_status'] == 'computed' for r in libraries),
            'every_library_mt_measured': all(r['mt_metric_status'] == 'measured' for r in libraries),
            'no_retained_predicted_doublet': statuses.get('predicted_doublet', 0) == 0,
        }
        result['shards'].append({
            'dataset': name, 'raw_ledger_rows': raw, 'retained_ledger_rows': kept,
            'filtered_reason_counts': dict(reasons), 'retained_doublet_status_counts': dict(statuses),
            'n_libraries': len(libraries), 'Zarr_X_shape': x_shape,
            'Zarr_counts_shape': counts_shape, 'checks': checks,
            'clinical_sidecar_present': bool(proof.get('clinical_sidecar_sha256')),
            'evidence_sha256': {p.name: sha256(p) for p in (proof_path, ledger, thresholds)},
            'evidence_directory': '${RUN_ROOT}/manifests',
        })
    if not result['shards']:
        raise ValueError('No registered QC shards found')
    result['total_retained_rows'] = sum(r['retained_ledger_rows'] for r in result['shards'])
    result['total_libraries'] = sum(r['n_libraries'] for r in result['shards'])
    result['all_checked_conditions_pass'] = all(all(r['checks'].values()) for r in result['shards'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.run.resolve())
    # Immutable observation: never overwrite an earlier audit.
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('total_retained_rows', 'total_libraries', 'all_checked_conditions_pass')}))
    raise SystemExit(0 if result['all_checked_conditions_pass'] else 1)


if __name__ == '__main__':
    main()
