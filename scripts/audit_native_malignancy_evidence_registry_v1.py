"""Inventory carried CNV/malignancy evidence without interpreting cell lineage.

Field presence, finite scores and negative flags are not interchangeable with a
completed CNV evaluation. Original fields/version/method remain authoritative.
"""
import argparse
import collections
import datetime
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256


PAIRS = [('ev_v7_infercnv_evaluable', 'ev_v7_infercnv_positive'),
         ('legacy_sample_infercnv_evaluable_v7', 'legacy_sample_infercnv_positive_v7'),
         ('cnv_evaluable', 'infercnv_positive')]
STATUS = ['legacy_cnv_role_v7', 'legacy_malignancy_call_v7', 'legacy_malignancy_tier_v7',
    'legacy_is_malignant_primary_v7', 'legacy_is_malignant_high_confidence_v7', 'legacy_is_malignant_sensitive_v7',
    'legacy_strict_malignant_downstream_eligible_v7', 'legacy_infercnv_run_status_v7', 'legacy_cnv_skip_reason_v7',
    'legacy_malignancy_evidence_version_v7', 'malignant_status', 'CNV_method_status', 'CNV_candidate_basis',
    'CNV_input', 'CNV_evidence_tier', 'malignancy_scope', 'cnv_condition']
SCORES = ['ev_v7_cnv_burden', 'legacy_sample_cnv_burden_v7', 'cnv_score', 'cnv_cors', 'cnv_cors_max',
          'cnv_ref_score', 'cnv_score_norm', 'cnv_score_norm_norm', 'cnv_score_norm_norm2']


def boolean_codes(values):
    values = np.char.lower(np.char.strip(np.asarray(values).astype(str)))
    codes = np.full(len(values), -2, np.int8)
    codes[np.isin(values, ['', 'nan', 'none', 'na', 'n/a', '<na>', 'unknown', 'not_assessed', 'not_evaluated'])] = -1
    codes[np.isin(values, ['true', '1', '1.0'])] = 1
    codes[np.isin(values, ['false', '0', '0.0'])] = 0
    return codes


def pair_counts(evaluable, positive):
    ev, pos = boolean_codes(evaluable), boolean_codes(positive)
    if ev.shape != pos.shape:
        raise ValueError('CNV flag axes differ')
    result = np.full(len(ev), 'evaluation_unknown_result_uninterpreted', dtype='<U70')
    result[ev == 0] = 'explicitly_not_evaluable_not_negative'
    result[(ev == 1) & (pos == 1)] = 'evaluable_positive_flag'
    result[(ev == 1) & (pos == 0)] = 'evaluable_negative_flag'
    result[(ev == 1) & (pos < 0)] = 'evaluable_result_unknown'
    result[(ev < 0) & (pos == 1)] = 'positive_flag_evaluation_status_unknown'
    result[(ev == 0) & (pos == 1)] = 'conflict_positive_flag_while_not_evaluable'
    result[(ev == -2) | (pos == -2)] = 'unrecognized_boolean_encoding_requires_review'
    return dict(collections.Counter(result))


def main(root):
    start = time.monotonic(); cfg = json.loads((root / 'config/malignancy_registry_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Native malignancy input changed')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_native_malignancy_registry_v1.py')], capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Native malignancy flag counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 3, 'no_skips': True})
    reader = NativeArchiveCollection(cfg['archive_manifest']); records = []
    for ds in reader.archive_partitions:
        rows = reader.rows(ds); n = len(rows); names = set(reader.field_names(ds))
        row = {'source': ds, 'n_archive_rows': n, 'native_rows_sha256': row_axis_sha256(rows),
            'native_ids_sha256': text_axis_sha256(reader.native_ids(ds)), 'flag_pair_evidence': [],
            'carried_statuses': {}, 'numeric_score_coverage': {}, 'new_malignancy_call': 'not_assigned',
            'lineage_used_to_determine_malignancy': False, 'prior_evidence_provenance_not_replaced': True}
        for ev, pos in PAIRS:
            if ev not in names or pos not in names:
                row['flag_pair_evidence'].append({'evaluation_field': ev, 'result_field': pos,
                    'status': 'missing_flag_pair_not_negative', 'available_fields': [f for f in [ev, pos] if f in names]})
                continue
            a, b = reader.field(ds, ev), reader.field(ds, pos)
            counts = pair_counts(a, b)
            if sum(counts.values()) != n:
                raise ValueError('CNV flag pair denominator differs')
            # Scalar cross-tab checks the vectorized tally using the same rule.
            pairs = collections.Counter(zip(np.asarray(a).astype(str), np.asarray(b).astype(str)))
            scalar = collections.Counter()
            for (left, right), number in pairs.items():
                scalar.update({next(iter(pair_counts([left], [right]))): number})
            if dict(scalar) != counts:
                raise ValueError('CNV category cross-tab differs')
            row['flag_pair_evidence'].append({'evaluation_field': ev, 'result_field': pos,
                'status': 'carried_flags_separated_method_scope_not_requalified', 'counts': counts})
        for field in STATUS:
            if field not in names:
                continue
            values = reader.field(ds, field).astype(str)
            row['carried_statuses'][field] = dict(collections.Counter(values))
            if sum(row['carried_statuses'][field].values()) != n:
                raise ValueError('Carried malignancy status denominator differs')
        for field in SCORES:
            if field not in names:
                continue
            import pandas as pd
            values = pd.to_numeric(pd.Series(reader.field(ds, field)), errors='coerce').to_numpy(dtype=float)
            row['numeric_score_coverage'][field] = {'n_finite': int(np.isfinite(values).sum()),
                'n_nonfinite_or_non_numeric': int((~np.isfinite(values)).sum()),
                'finite_value_does_not_certify_evaluable_CNV': True}
        records.append(row)
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source': ds, 'completed_sources': len(records)})
        reader._fields.clear()
    output = root / 'outputs/native_malignancy_evidence_registry_v1.json'
    save(output, {'sources': records, 'n_sources': len(records), 'native_expression_read': False,
        'new_CNV_analysis': False, 'individual_records_exported': False, 'new_authoritative_labels': False,
        'method_scope_gene_order_controls_and_original_authority_review_still_required': True})
    if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
        raise ValueError('Malignancy metadata cap exceeded')
    save(root / 'provenance/terminal_v1.json', {'status': 'CARRIED_NATIVE_MALIGNANCY_FLAG_MEASUREMENT_REGISTRY_COMPLETE_NOT_RECLASSIFICATION',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(records),
        'output_sha256': sha256_file(output), 'native_expression_read': False, 'new_authoritative_labels': False,
        'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
