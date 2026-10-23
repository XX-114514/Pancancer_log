"""Summarize audited native evidence before promoting routes or identities.

This reads compact evidence sidecars, not expression matrices. Multi-family
support is not a doublet call; missing reference coverage is not an unknown type.
"""
import argparse
import datetime
import json
from pathlib import Path

import numpy as np

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file


def accepted(producer, audit, index_name):
    for root in [producer, audit]:
        result = json.loads((root / 'provenance/supervisor_terminal_v1.json').read_text())
        if result['actual_child_exit_code'] != 0 or result['wall_timeout']:
            raise ValueError('Both producer and independent auditor must succeed')
    terminal = json.loads((producer / 'provenance/terminal_v1.json').read_text())
    checked = json.loads((audit / 'provenance/terminal_v1.json').read_text())
    if checked['producer_terminal_sha256'] != sha256_file(producer / 'provenance/terminal_v1.json'):
        raise ValueError('Audit belongs to a different producer terminal')
    path = producer / 'outputs' / index_name
    if sha256_file(path) != terminal['index_sha256']:
        raise ValueError('Accepted index changed')
    records = json.loads(path.read_text())
    if len({r['source'] for r in records}) != len(records):
        raise ValueError('Duplicate source in evidence index')
    return {r['source']: r for r in records}, terminal


def read_source(root, record):
    for field in ['arrays', 'metadata']:
        if sha256_file(root / 'outputs' / record[field]) != record[field + '_sha256']:
            raise ValueError('Accepted source artifact changed')
    meta = json.loads((root / 'outputs' / record['metadata']).read_text())
    with np.load(root / 'outputs' / record['arrays'], allow_pickle=False) as stored:
        values = {key: stored[key] for key in stored.files}
    if any(len(v) != record['n_archive_rows'] for v in values.values()):
        raise ValueError('Sidecar row denominator changed')
    return meta, values


def main(args):
    refs, rt = accepted(args.reference, args.reference_audit, 'full_archive_reference_index_v1.json')
    cues, ct = accepted(args.supplemental, args.supplemental_audit, 'supplemental_native_index_v1.json')
    if set(refs) != set(cues) or rt['n_archive_rows'] != ct['n_archive_rows']:
        raise ValueError('Full native evidence source scope differs')
    reports = []
    for ds in sorted(refs):
        rm, rv = read_source(args.reference, refs[ds]); cm, cv = read_source(args.supplemental, cues[ds])
        if rm['native_ids_sha256'] != cm['native_ids_sha256']:
            raise ValueError('Native identity axes differ')
        np.testing.assert_array_equal(rv['native_row'], cv['native_row'])
        np.testing.assert_array_equal(rv['native_total_RNA'], cv['native_total_RNA'])
        positive = rv['positive_RNA']; bits = rv['native_family_bitmask']
        families = rm['native_family_bit_order']
        n_support = np.zeros(len(bits), np.uint8)
        family_counts = {}
        for j, family in enumerate(families):
            mask = (bits & (1 << j)) != 0
            family_counts[family] = int(mask.sum()); n_support += mask
        if np.any(bits.astype(np.uint64) >> len(families)):
            raise ValueError('Out-of-range family evidence bit')
        reference_counts = {}
        for model in rm['source_reference_models']:
            scope = rv['immune_reference_scope'] if model.startswith('immune_') else rv['organ_reference_scope']
            n = rv[model + '__n_positive_labels']
            if np.any((n >= 0) & ~scope):
                raise ValueError('Reference predictions outside frozen scope')
            reference_counts[model] = {
                'n_in_scope': int(scope.sum()),
                'n_in_scope_unscored': int((scope & (n < 0)).sum()),
                'n_scored_no_positive': int((n == 0).sum()),
                'n_scored_unique_positive': int((n == 1).sum()),
                'n_scored_multiple_positive': int((n > 1).sum()),
                'feature_coverage': rm['reference_feature_coverage'][model]}
        extra = {}
        for j, name in enumerate(cm['panel_order']):
            cue = cv['cue_support'][:, j]
            extra[name] = {
                'n_positive': int(cue.sum()),
                'n_positive_with_no_base_family_support': int((cue & (n_support == 0)).sum()),
                'n_positive_with_one_base_family_support': int((cue & (n_support == 1)).sum()),
                'n_positive_with_multiple_base_family_support': int((cue & (n_support > 1)).sum())}
        reports.append({'source': ds, 'n_archive_rows': len(bits), 'n_positive_RNA': int(positive.sum()),
            'n_nonmissing_author_labels': rm['n_nonmissing_author_labels'],
            'native_support_counts_not_identity_labels': family_counts,
            'positive_RNA_support_multiplicity': {str(k): int((positive & (n_support == k)).sum()) for k in np.unique(n_support[positive])},
            'reference_proposals_not_truth': reference_counts,
            'n_positive_RNA_without_qualified_reference_scope': int((positive & ~rv['immune_reference_scope'] & ~rv['organ_reference_scope']).sum()),
            'supplemental_cues_not_identity_or_malignancy': extra,
            'new_identity_assigned': False})
        save(args.output / 'checkpoints' / (ds + '.json'), {
            'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'source': ds,
            'completed_sources': len(reports)})
    if sum(r['n_archive_rows'] for r in reports) != rt['n_archive_rows']:
        raise ValueError('Native summary denominator differs')
    result = args.output / 'outputs/native_annotation_evidence_coverage_v1.json'
    save(result, {'sources': reports, 'native_expression_read': False, 'new_authoritative_labels': False,
        'multiple_family_support_is_not_doublet_call': True,
        'reference_unavailable_is_not_final_unknown_identity': True,
        'identity_malignancy_and_model_eligibility_still_require_review': True})
    save(args.output / 'provenance/terminal_v1.json', {
        'status': 'AUDITED_NATIVE_EVIDENCE_COVERAGE_SUMMARIZED_NOT_FINAL_ANNOTATION',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'output_sha256': sha256_file(result), 'n_sources': len(reports),
        'new_authoritative_labels': False, 'native_expression_read': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['reference', 'reference-audit', 'supplemental', 'supplemental-audit', 'output']:
        parser.add_argument('--' + name, required=True, type=Path)
    args = parser.parse_args()
    try:
        main(args)
    except Exception as error:
        save(args.output / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
