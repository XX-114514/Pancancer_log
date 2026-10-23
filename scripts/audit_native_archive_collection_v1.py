"""Independent membership audit and native readback for a selected-source archive.

Every native/source/cell/unit axis is checked. All newly restored rows receive
a count-domain scan; per-source representative panels are read directly from
CSR and compared with the public reader. Existing context count scans are a
separate bound evidence stage, not rerun or overstated here.
"""
import argparse
import collections
import datetime
import hashlib
import json
from pathlib import Path
import resource
import time

import numpy as np
import scipy.sparse as sp
from anndata._io.specs import read_elem

from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False);handle.write('\n')


def direct_csr(group, start, end):
    ptr = np.asarray(group['indptr'][int(start):int(end) + 1], dtype=np.int64)
    low, high = int(ptr[0]), int(ptr[-1])
    return sp.csr_matrix((np.asarray(group['data'][low:high]), np.asarray(group['indices'][low:high]), ptr - low),
                         shape=(end - start, group.attrs['shape'][1]))


def main(args):
    start = time.monotonic();root = Path(args.output)
    reader = NativeArchiveCollection(args.manifest)
    context = FullContextAtlas(reader.archive_manifest['frozen_context']['path'])
    context_evidence = Path(args.context_evidence)
    previous = json.loads((context_evidence / 'provenance/supervisor_terminal_v1.json').read_text())
    if previous['actual_child_exit_code'] != 0 or previous['wall_timeout']:
        raise ValueError('Previous context evidence not successful')
    expected_total, context_total, extra_total, checked_values = 0, 0, 0, 0
    sources = []
    reason_fields = ['authoritative_l1_exclusion_reason', 'expression_v1_exclusion_reason',
                     'cell_class_status', 'doublet_status', 'v9_sample_flag',
                     'annotation_exclude_from_primary', 'residual_doublet_flag',
                     'preintegration_cell_eligible', 'downstream_release_eligible']
    for dataset in reader.partitions:
        z = reader.store(dataset);entry = reader.archive_partitions[dataset]
        raw_source = np.asarray(read_elem(z['obs/' + entry['source_field']]), dtype=object)
        expected = np.flatnonzero(raw_source == dataset)
        rows = reader.rows(dataset)
        if not np.array_equal(expected, rows):raise ValueError('Archive omitted or reordered native source rows')
        ids = reader.native_ids(dataset);units = reader.frozen_expression_units(dataset)
        if len(set(ids)) != len(ids):raise ValueError('Repeated native ID within source namespace')
        old = context.rows(dataset)
        if not np.array_equal(old, reader.context_rows(dataset)):raise ValueError('Historical context changed')
        where_old = np.searchsorted(rows, old)
        if not np.array_equal(rows[where_old], old):raise ValueError('Invalid historical membership')
        if not np.array_equal(units[where_old], context.frozen_expression_units(dataset)):
            raise ValueError('Frozen expression unit keys changed')
        old_labels, labels = context.frozen_labels(dataset), reader.frozen_labels(dataset)
        for kind in ['major', 'fine']:
            if not np.array_equal(labels[kind][where_old], old_labels[kind]):
                raise ValueError('Inherited context label changed')
        membership = reader.membership(dataset)
        if membership['in_frozen_context'].sum() != len(old):raise ValueError('Context flag denominator differs')
        if reader.is_legacy(dataset):
            native_roster = np.isin(rows, context._positions)
            if not np.array_equal(native_roster, membership['in_legacy_roster']):
                raise ValueError('Roster membership differs')
        lookup = reader.feature_lookup(dataset)
        # Deterministic measured genes and one guaranteed unmeasured feature.
        available = sorted(g for g, cols in lookup.items() if len(cols) == 1)
        genes, used = [], set()
        for gene in available:
            col = next(iter(lookup[gene]))
            if col not in used:genes.append(gene);used.add(col)
            if len(genes) == 6:break
        absent = '__ARCHIVE_UNMEASURED_SENTINEL__'
        if absent in lookup:raise ValueError('Unexpected measurement sentinel collision')
        genes.append(absent)
        selected = np.asarray(sorted(set([0, len(rows) // 2, len(rows) - 1])), dtype=np.int64)
        # Include a duplicate to prove API order preservation in real sources.
        selected = np.r_[selected, selected[0]]
        panel = reader.read_panel(dataset, selected, genes)
        csr = z[entry['count_origin']]
        native = sp.vstack([direct_csr(csr, int(r), int(r) + 1) for r in rows[selected]], format='csr')
        cols = [next(iter(lookup[g])) for g in genes[:-1]]
        wanted = sp.hstack([native[:, cols], sp.csr_matrix((len(selected), 1))], format='csr')
        if (panel['counts'] != wanted).nnz:raise ValueError('Representative native panel readback differs')
        if not np.array_equal(panel['native_total_counts'], np.asarray(native.sum(axis=1, dtype=np.float64)).ravel()):
            raise ValueError('Native denominator readback differs')
        if panel['measured'].tolist() != [True] * len(cols) + [False]:raise ValueError('Missing feature lost its mask')
        extra = rows[~membership['in_frozen_context']]
        domains = {'n_rows': len(extra), 'n_zero_RNA_rows': 0, 'n_nonfinite_or_negative_rows': 0,
                   'n_noninteger_rows': 0, 'n_stored_nnz_checked': 0,
                   'integer_domain_is_not_UMI_provenance_acceptance': True}
        digest = hashlib.sha256()
        # Read bounded native spans, then subset; avoids one Zarr request per
        # isolated row and does not densify the full feature axis.
        for block in np.unique(extra // 2048):
            physical = extra[extra // 2048 == block]
            lo, hi = int(physical.min()), int(physical.max()) + 1
            x = direct_csr(csr, lo, hi)[physical - lo].tocsr()
            domains['n_stored_nnz_checked'] += len(x.data)
            totals = np.asarray(x.sum(axis=1, dtype=np.float64)).ravel()
            domains['n_zero_RNA_rows'] += int((totals == 0).sum())
            for a, b in zip(x.indptr[:-1], x.indptr[1:]):
                values = x.data[a:b]
                domains['n_nonfinite_or_negative_rows'] += int((~np.isfinite(values) | (values < 0)).any())
                domains['n_noninteger_rows'] += int((np.isfinite(values) & (values != np.floor(values))).any())
            for a in [physical, x.indptr, x.indices, x.data]:
                digest.update(str(a.dtype).encode());digest.update(str(a.shape).encode());digest.update(a.tobytes())
        domains['supplementary_native_CSR_digest_sha256'] = digest.hexdigest()
        reason_summary = {}
        fields = reader.field_names(dataset)
        for field in reason_fields:
            if field in fields:
                values = reader.field(dataset, field)
                reason_summary[field] = {
                    str(group): dict(collections.Counter(str(v) for v in values[membership['membership_status'] == group]))
                    for group in sorted(set(membership['membership_status']))}
        sources.append({'source': dataset, 'n_archive_rows': len(rows), 'n_context_rows': len(old),
                        'n_expression_unit_keys': len(set(units)),
                        'native_row_sha256': row_axis_sha256(rows), 'native_id_sha256': text_axis_sha256(ids),
                        'native_unit_sha256': text_axis_sha256(units),
                        'membership_counts': dict(collections.Counter(membership['membership_status'])),
                        'carried_reason_flag_counts': reason_summary,
                        'supplementary_count_domain': domains,
                        'all_native_axes_and_representative_CSR_readback_pass': True,
                        'new_identity_or_model_eligibility_accepted': False})
        expected_total += len(rows);context_total += len(old);extra_total += len(extra)
        checked_values += domains['n_stored_nnz_checked']
        save(root / 'checkpoints' / (dataset + '.json'), {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'stage': 'native_archive_source_readback_complete', 'source': dataset})
    output = root / 'outputs/private_source_archive_audit_v1.json'
    save(output, {'sources': sources, 'archive_scope': reader.archive_manifest['archive_scope']})
    terminal = {'status': 'NATIVE_SELECTED_SOURCE_ARCHIVE_MEMBERSHIP_AND_READBACK_PASS_NOT_NEW_MODEL_ELIGIBILITY',
                'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'manifest_sha256': sha256_file(args.manifest), 'source_summary_sha256': sha256_file(output),
                'n_sources': len(sources), 'n_archive_rows': expected_total, 'n_frozen_context_rows': context_total,
                'n_supplementary_rows': extra_total, 'n_supplementary_stored_nnz_checked': checked_values,
                'all_source_axes_checked': True, 'all_new_rows_count_domain_scanned': True,
                'count_readback_scope': 'all supplementary rows plus per-source representative native panels; previous frozen context evidence bound separately',
                'previous_context_evidence_terminal_sha256': sha256_file(context_evidence / 'provenance/terminal_v1.json'),
                'upstream_QC_removed_droplets_reconstructed': False, 'standalone_Zarr': False,
                'new_identity_or_model_or_clinical_eligibility_accepted': False,
                'wall_seconds': time.monotonic() - start,
                'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024}
    save(root / 'provenance/terminal_v1.json', terminal)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--context-evidence', required=True)
    parser.add_argument('--output', required=True)
    main(parser.parse_args())
