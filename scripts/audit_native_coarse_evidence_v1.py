"""Independently check complete evidence axes and replay one whole native source.

The replay reads CSR arrays directly, not the producer's gather/normalization
helpers. Scientific summaries remain private. Evidence acceptance is separate
from biological annotation acceptance.
"""
import argparse
import collections
import datetime
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from full_context_atlas_v1 import FullContextAtlas, sha256_file


def main(args):
    root, evidence = Path(args.output), Path(args.evidence)
    terminal = json.loads((evidence / 'provenance/terminal_v1.json').read_text())
    supervisor = json.loads((evidence / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervisor['actual_child_exit_code'] != 0 or supervisor['wall_timeout']:
        raise ValueError('Evidence scientific process did not finish successfully')
    if terminal['status'] != 'FULL_NATIVE_COARSE_MARKER_EVIDENCE_MEASURED_NOT_NEW_ANNOTATION':
        raise ValueError('Unexpected evidence stage')
    launch = json.loads((evidence / 'provenance/launch_receipt_v1.json').read_text())
    for name, digest in launch['script_SHA256'].items():
        if sha256_file(evidence / name) != digest:
            raise ValueError('Changed executed evidence script')
    index_path = evidence / 'outputs/source_evidence_index_v1.json'
    if sha256_file(index_path) != terminal['output_index_sha256']:
        raise ValueError('Changed evidence index')
    index = json.loads(index_path.read_text())
    reader = FullContextAtlas(args.manifest)
    if {x['dataset'] for x in index} != set(reader.partitions) or len(index) != len(reader.partitions):
        raise ValueError('Evidence source axis differs')
    checked, missing_cases, zero_cases = [], [], []
    for item in index:
        path = evidence / item['file']
        if sha256_file(path) != item['sha256']:
            raise ValueError('Changed evidence source output')
        data = json.loads(path.read_text())
        dataset = item['dataset']
        n = len(reader.rows(dataset))
        labels = reader.frozen_labels(dataset)
        expected = collections.Counter(str(a) if str(a) != 'unresolved' else 'unresolved::author::' + str(b)
                                       for a, b in zip(labels['major'], labels['fine']))
        observed = {g['inherited_group']: g['n_archive_cells'] for g in data['groups']}
        if observed != dict(expected) or sum(observed.values()) != n or n != item['n_archive_cells']:
            raise ValueError('Evidence inherited-group/count reconciliation failed')
        lookup = reader.feature_lookup(dataset)
        measured = np.array([len(lookup.get(g, set())) == 1 for g in data['genes']])
        if measured.tolist() != data['measured']:
            raise ValueError('Evidence source measurement mask changed')
        for group in data['groups']:
            evaluable = group['n_native_RNA_evaluable']
            if not 0 <= evaluable <= group['n_archive_cells']:
                raise ValueError('Invalid native RNA denominator')
            for j, gene in enumerate(data['genes']):
                values = [group[k][j] for k in ['n_detected', 'mean_native_CP10K', 'mean_native_log1p_CP10K']]
                if not measured[j] or not evaluable:
                    if values != [None, None, None]:raise ValueError('Missing measurement treated as zero')
                    if len(missing_cases) < 3:missing_cases.append({'dataset': dataset, 'gene': gene})
                elif (any(v is None or not np.isfinite(v) or v < 0 for v in values) or
                      values[0] > evaluable or values[1] > 10000 + 1e-6):
                    raise ValueError('Native expression evidence out of bounds')
                elif values[0] == 0:
                    if values[1:] != [0, 0]:raise ValueError('True zero marker has nonzero mean')
                    if len(zero_cases) < 3:zero_cases.append({'dataset': dataset, 'gene': gene})
        if sum(g['n_native_RNA_evaluable'] for g in data['groups']) != item['n_native_RNA_evaluable']:
            raise ValueError('Evaluable denominator not reconciled')
        checked.append({'dataset': dataset, 'n_cells': n, 'axes_bounds_nulls_pass': True})
    # Replay the smallest complete source with at least two inherited groups.
    candidate = min((x for x in index if x['n_groups'] >= 2), key=lambda x: x['n_archive_cells'])
    ds = candidate['dataset']
    native = reader.rows(ds)
    data = json.loads((evidence / candidate['file']).read_text())
    z = reader.store(ds)
    key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
    csr = z[key]
    shape = tuple(csr.attrs['shape'])
    lookup = reader.feature_lookup(ds)
    columns = [next(iter(lookup[g])) if len(lookup.get(g, set())) == 1 else None for g in data['genes']]
    lab = reader.frozen_labels(ds)
    names = np.asarray([str(a) if str(a) != 'unresolved' else 'unresolved::author::' + str(b)
                        for a, b in zip(lab['major'], lab['fine'])])
    actual = {}
    for group in data['groups']:
        actual[group['inherited_group']] = {'n': 0, 'det': np.zeros(len(columns), np.int64),
                                          'cp': np.zeros(len(columns)), 'log': np.zeros(len(columns))}
    for pos, physical in enumerate(native):
        low, high = np.asarray(csr['indptr'][physical:physical + 2], dtype=np.int64)
        vals = np.asarray(csr['data'][low:high], dtype=np.float64)
        cols = np.asarray(csr['indices'][low:high])
        total = vals.sum()
        if not total:continue
        # Sum possible repeated CSR column entries before detection/logging.
        row = sp.csr_matrix((vals, cols, [0, len(vals)]), shape=(1, shape[1]))
        row.sum_duplicates()
        vector = np.array([float(row[0, c]) if c is not None else 0 for c in columns])
        a = actual[names[pos]];a['n'] += 1
        a['det'] += vector > 0
        a['cp'] += vector / total * 10000
        a['log'] += np.log1p(vector / total * 10000)
    for group in data['groups']:
        a = actual[group['inherited_group']]
        if a['n'] != group['n_native_RNA_evaluable']:raise ValueError('Independent native total differs')
        for j, c in enumerate(columns):
            if c is not None and a['n']:
                if a['det'][j] != group['n_detected'][j]:raise ValueError('Independent detection differs')
                np.testing.assert_allclose(a['cp'][j] / a['n'], group['mean_native_CP10K'][j], rtol=1e-10, atol=1e-10)
                np.testing.assert_allclose(a['log'][j] / a['n'], group['mean_native_log1p_CP10K'][j], rtol=1e-10, atol=1e-10)
    if sum(x['n_cells'] for x in checked) != terminal['n_archive_cells']:
        raise ValueError('Full source counts differ')
    result = {'status': 'INDEPENDENT_NATIVE_MARKER_EVIDENCE_AUDIT_PASS_NOT_IDENTITY',
              'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'n_checked_sources': len(checked), 'n_checked_archive_cells': terminal['n_archive_cells'],
              'all_source_group_masks_nulls_and_bounds_pass': True,
              'whole_source_direct_CSR_replay': {'dataset': ds, 'n_cells': len(native), 'passed': True},
              'true_zero_examples_found': len(zero_cases), 'unmeasured_examples_found': len(missing_cases),
              'zero_RNA_real_examples': terminal['n_archive_cells'] - terminal['n_native_RNA_evaluable'],
              'zero_RNA_synthetic_counterexample_already_passed': True,
              'terminal_SHA256': sha256_file(evidence / 'provenance/terminal_v1.json'),
              'index_SHA256': sha256_file(index_path), 'identity_or_state_accepted': False}
    with (root / 'provenance/terminal_v1.json').open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False);f.write('\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['output', 'evidence', 'manifest']:p.add_argument('--' + name, required=True)
    main(p.parse_args())
