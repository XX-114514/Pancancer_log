"""Native evidence for rows restored to the archive, without promoting labels."""
import argparse
import collections
import datetime
import json
from pathlib import Path
import resource
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection
from native_panel_aggregate_v1 import GENES, LINEAGE_PANELS, STATE_PANELS, IDENTITY_CANDIDATE_PANELS, MarkerEvidence, native_marker_block


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:json.dump(value, f, indent=2, allow_nan=False);f.write('\n')


def main(args):
    begin = time.monotonic();root, archive = Path(args.output), Path(args.archive)
    audit = json.loads((archive / 'provenance/terminal_v1.json').read_text())
    supervisor = json.loads((archive / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervisor['actual_child_exit_code'] != 0 or supervisor['wall_timeout']:raise ValueError('Native archive not accepted')
    manifest = archive / 'config/native_archive_collection_v1.json'
    if sha256_file(manifest) != audit['manifest_sha256']:raise ValueError('Archive manifest changed')
    reader = NativeArchiveCollection(manifest)
    summaries = []
    smallest, smallest_n = None, None
    for ds in reader.partitions:
        n = int((~reader.membership(ds)['in_frozen_context']).sum())
        if n and (smallest_n is None or n < smallest_n):smallest, smallest_n = ds, n
    contract = {'archive_manifest_sha256': sha256_file(manifest), 'genes': GENES,
                'lineage_panels': LINEAGE_PANELS, 'state_panels': STATE_PANELS,
                'identity_candidate_panels': IDENTITY_CANDIDATE_PANELS,
                'scope': 'all selected-source native rows outside the frozen analysis context',
                'group_definition': 'archival membership status and inherited coarse label, not a new cluster',
                'normalization': 'native counts / all stored native RNA * 10000; per-cell log1p',
                'zero_RNA_excluded_from_expression_denominator_but_retained_in_archive': True,
                'unmeasured_is_null_not_zero': True,
                'independent_direct_CSR_replay_source_selected_by_smallest_supplementary_count': smallest,
                'new_identity_state_or_malignancy_calls': False}
    save(root / 'config/supplementary_native_evidence_contract_v1.json', contract)
    for ds in reader.partitions:
        membership = reader.membership(ds)
        selected = np.flatnonzero(~membership['in_frozen_context'])
        if not len(selected):continue
        inherited = reader.frozen_labels(ds)['major']
        keys = [json.dumps([str(membership['membership_status'][i]), str(inherited[i])]) for i in selected]
        codes, names = pd.factorize(keys, sort=True)
        lookup = reader.feature_lookup(ds)
        mask = np.asarray([len(lookup.get(g, set())) == 1 for g in GENES])
        evidence = MarkerEvidence(names.tolist(), GENES, mask)
        nnz = 0
        for offset in range(0, len(selected), 2048):
            block = native_marker_block(reader, ds, selected[offset:offset + 2048])
            evidence.update(block['counts'], block['native_total_counts'], codes[offset:offset + 2048])
            nnz += block['full_native_nnz_read']
        if evidence.n_cells.sum() != len(selected):raise ValueError('Supplementary evidence denominator differs')
        records = evidence.records()
        # Independent scalar replay, using direct native CSR rather than the
        # producer's gather/normalization/aggregation helpers.
        if ds == smallest:
            actual = {name: {'n': 0, 'det': np.zeros(len(GENES), int), 'cp': np.zeros(len(GENES)), 'log': np.zeros(len(GENES))} for name in names}
            csr = reader.store(ds)[reader.archive_partitions[ds]['count_origin']]
            cols = [next(iter(lookup[g])) if mask[j] else None for j, g in enumerate(GENES)]
            for index, physical in enumerate(reader.rows(ds)[selected]):
                low, high = np.asarray(csr['indptr'][int(physical):int(physical) + 2], dtype=np.int64)
                row = sp.csr_matrix((np.asarray(csr['data'][low:high]), np.asarray(csr['indices'][low:high]), [0, high - low]), shape=(1, csr.attrs['shape'][1]))
                row.sum_duplicates();total = float(row.sum())
                if total <= 0:continue
                values = np.array([float(row[0, c]) if c is not None else 0. for c in cols])
                group = actual[keys[index]];group['n'] += 1;group['det'] += values > 0
                group['cp'] += values * 10000 / total;group['log'] += np.log1p(values * 10000 / total)
            for record in records:
                group = actual[record['frozen_group']]
                if group['n'] != record['n_native_RNA_evaluable']:raise ValueError('Independent supplemental denominator differs')
                for j, gene in enumerate(GENES):
                    value = record['marker_evidence'][gene]
                    if not mask[j] or group['n'] == 0:
                        if any(value[k] is not None for k in ['n_detected','detection_fraction','mean_native_CP10K','mean_native_log1p_CP10K']):
                            raise ValueError('Supplementary missing measurement treated as zero')
                    elif (value['n_detected'] != group['det'][j] or
                          not np.isclose(value['mean_native_CP10K'], group['cp'][j] / group['n'], rtol=1e-9, atol=1e-9) or
                          not np.isclose(value['mean_native_log1p_CP10K'], group['log'][j] / group['n'], rtol=1e-9, atol=1e-9)):
                        raise ValueError('Independent supplementary marker numeric replay differs')
        path = root / 'outputs' / (ds + '_supplementary_evidence.json')
        save(path, {'source': ds, 'genes': GENES, 'measured': mask.tolist(), 'records': records,
                     'n_archive_supplementary_rows': len(selected), 'native_nnz_checked': nnz,
                     'previous_exclusion_reason_counts': dict(collections.Counter(str(v) for v in reader.field(ds, 'expression_v1_exclusion_reason')[selected])) if reader.is_legacy(ds) else {},
                     'new_label_or_model_eligibility_assigned': False})
        summaries.append({'source': ds, 'file': str(path.relative_to(root)), 'sha256': sha256_file(path),
                           'n_rows': len(selected), 'n_groups': len(names), 'n_native_RNA_evaluable': int(evidence.n_evaluable.sum()),
                           'native_nnz_checked': nnz})
        save(root / 'checkpoints' / (ds + '.json'), {'source': ds, 'stage': 'supplementary_native_evidence_complete',
             'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat()})
    if sum(s['n_rows'] for s in summaries) != audit['n_supplementary_rows']:raise ValueError('Archive/evidence scope differs')
    index = root / 'outputs/supplementary_evidence_index_v1.json';save(index, summaries)
    save(root / 'provenance/terminal_v1.json', {
        'status': 'SUPPLEMENTARY_NATIVE_MARKER_EVIDENCE_COMPLETE_WITH_SIZE_SELECTED_DIRECT_REPLAY_NOT_NEW_IDENTITY',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'n_supplementary_rows': sum(s['n_rows'] for s in summaries), 'n_sources_with_supplement': len(summaries),
        'n_native_RNA_evaluable': sum(s['n_native_RNA_evaluable'] for s in summaries),
        'n_native_nnz_checked': sum(s['native_nnz_checked'] for s in summaries),
        'index_sha256': sha256_file(index), 'direct_CSR_numeric_replay_source': smallest,
        'direct_replay_scope': 'all supplemental rows in the size-selected source, not all sources',
        'all_other_source_evidence_independent_numeric_replay': False,
        'new_identity_or_model_eligibility_assigned': False,
        'wall_seconds': time.monotonic() - begin,
        'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})
    if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > 8 * 1024 ** 2:
        raise ValueError('Supplementary evidence output cap exceeded')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--archive', required=True);parser.add_argument('--output', required=True)
    main(parser.parse_args())
