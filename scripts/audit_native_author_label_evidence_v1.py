"""Inventory carried author labels without assigning new cell identities.

Only explicit class names in a frozen vocabulary receive a coarse context.
All remaining original values stay in source-private vocabulary reports. An
author malignancy label is an evidence axis separate from the lineage context.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

from annotation_representation_pilot_v1 import ROUTES, save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, text_axis_sha256
from audit_supplemental_native_context_v1 import context_counts


MISSING = {'', 'na', 'n/a', 'nan', 'none', '<na>', 'unknown', 'unassigned', 'not available',
    'not_available', 'not provided', '__source_value_missing__', '__missing__',
    'author_label_missing', 'source_missing'}
VOCABULARY = {
    'T_NK': ['t', 't cell', 't cells', 't-cell', 't-cells', 't lymphocyte', 't lymphocytes',
        'cd4 t', 'cd4 t cell', 'cd4 t cells', 'cd8 t', 'cd8 t cell', 'cd8 t cells',
        'regulatory t cell', 'regulatory t cells', 'treg', 'tregs', 'gamma-delta t cell',
        'gamma-delta t cells', 'gd t', 'nk', 'nk cell', 'nk cells', 'natural killer cell',
        'natural killer cells', 't/nk', 't/nk cells', 't and nk cells',
        'cd4-positive, alpha-beta t cell', 'cd8-positive, alpha-beta t cell'],
    'B_Plasma': ['b', 'b cell', 'b cells', 'b-cell', 'b-cells', 'b lymphocyte', 'b lymphocytes',
        'plasma', 'plasma cell', 'plasma cells', 'plasmablast', 'plasmablasts'],
    'Myeloid': ['myeloid', 'myeloid cell', 'myeloid cells', 'macrophage', 'macrophages',
        'monocyte', 'monocytes', 'microglia', 'microglial cell', 'dendritic cell',
        'dendritic cells', 'dc', 'dcs', 'neutrophil', 'neutrophils',
        'classical monocyte', 'non-classical monocyte', 'conventional dendritic cell',
        'plasmacytoid dendritic cell', 'alveolar macrophage', 'kupffer cell'],
    'Mast': ['mast', 'mast cell', 'mast cells'],
    'Endothelial': ['endothelial', 'endothelial cell', 'endothelial cells',
        'vascular endothelial cell', 'lymphatic endothelial cell'],
    'Fibroblast': ['fibroblast', 'fibroblasts', 'fibroblast cell', 'fibroblast cells',
        'cancer associated fibroblast', 'cancer-associated fibroblast'],
    'Mural': ['mural', 'mural cell', 'mural cells', 'pericyte', 'pericytes',
        'smooth muscle cell', 'smooth muscle cells', 'vascular smooth muscle cell'],
    'Epithelial': ['epithelial', 'epithelial cell', 'epithelial cells',
        'luminal epithelial cell', 'basal epithelial cell', 'type i pneumocyte',
        'type ii pneumocyte', 'club cell', 'goblet cell', 'ciliated epithelial cell'],
    'Glial': ['oligodendrocyte', 'oligodendrocytes', 'oligodendrocyte precursor cell',
        'oligodendrocyte progenitor cell', 'glial cell', 'glial cells'],
    'Astrocyte': ['astrocyte', 'astrocytes'],
    'Neuron': ['neuron', 'neurons', 'neuronal cell', 'excitatory neuron', 'inhibitory neuron'],
    'Erythroid': ['erythroid', 'erythrocyte', 'erythrocytes', 'erythroid cell', 'erythroid cells'],
    'Hepatocyte': ['hepatocyte', 'hepatocytes'],
    'Melanocyte': ['melanocyte', 'melanocytes'],
    'Mesothelial': ['mesothelial', 'mesothelial cell', 'mesothelial cells'],
    'Adipocyte': ['adipocyte', 'adipocytes'],
    'Schwann': ['schwann cell', 'schwann cells'],
    'Neuroendocrine': ['neuroendocrine cell', 'neuroendocrine cells'],
    'Megakaryocyte_platelet': ['megakaryocyte', 'megakaryocytes', 'platelet', 'platelets']}
MALIGNANCY_LABELS = {'tumor', 'tumour', 'tumor cell', 'tumor cells', 'tumour cell',
    'tumour cells', 'malignant', 'malignant cell', 'malignant cells', 'cancer cell', 'cancer cells'}
LOOKUP = {name: family for family, names in VOCABULARY.items() for name in names}
if len(LOOKUP) != sum(map(len, VOCABULARY.values())):
    raise ValueError('Author vocabulary contains duplicate classes')


def author_context(value):
    """Do not strip subtype/marker suffixes or infer a lineage from malignancy."""
    value = ' '.join(str(value).strip().lower().split())
    if value in MISSING:
        return 'missing_or_explicit_unknown'
    if value in MALIGNANCY_LABELS:
        return 'author_malignancy_label_lineage_unspecified'
    return LOOKUP.get(value, 'unreviewed_author_value')


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/author_evidence_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Author evidence input changed')
    producer, audit = Path(cfg['producer']), Path(cfg['audit'])
    for path in [producer, audit]:
        terminal = json.loads((path / 'provenance/supervisor_terminal_v1.json').read_text())
        if terminal['actual_child_exit_code'] != 0 or terminal['wall_timeout']:
            raise ValueError('Accepted supplemental producer and independent audit required')
    checked = json.loads((audit / 'provenance/terminal_v1.json').read_text())
    if checked['producer_terminal_sha256'] != sha256_file(producer / 'provenance/terminal_v1.json'):
        raise ValueError('Supplemental audit belongs to another producer')
    terminal = json.loads((producer / 'provenance/terminal_v1.json').read_text())
    index_path = producer / 'outputs/supplemental_native_index_v1.json'
    if sha256_file(index_path) != terminal['index_sha256']:
        raise ValueError('Accepted supplemental index changed')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_native_author_label_evidence_v1.py')],
        capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 5 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Author evidence counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 5, 'no_skips': True})
    index = json.loads(index_path.read_text()); reader = NativeArchiveCollection(cfg['archive_manifest'])
    if len(index) != len(reader.archive_partitions) or {r['source'] for r in index} != set(reader.archive_partitions):
        raise ValueError('Native author evidence source scope differs')
    reports = []
    for record in index:
        ds = record['source']; rows = reader.rows(ds); n = len(rows)
        for name in ['arrays', 'metadata']:
            if sha256_file(producer / 'outputs' / record[name]) != record[name + '_sha256']:
                raise ValueError('Supplemental artifact changed')
        meta = json.loads((producer / 'outputs' / record['metadata']).read_text())
        if meta['native_ids_sha256'] != text_axis_sha256(reader.native_ids(ds)) or n != record['n_archive_rows']:
            raise ValueError('Native author evidence identity axis differs')
        with np.load(producer / 'outputs' / record['arrays'], allow_pickle=False) as f:
            np.testing.assert_array_equal(f['native_row'], rows)
            support = f['cue_support']
        fields = set(reader.field_names(ds))
        field = 'ev_author_label' if reader.is_legacy(ds) else (
            reader.primary.partitions[ds].get('annotation', {}).get('field') or 'source_author_label')
        original = reader.field(ds, field).astype(str) if field in fields else np.full(n, '', dtype='<U1')
        values, inverse, counts = np.unique(original, return_inverse=True, return_counts=True)
        contexts = np.asarray([author_context(v) for v in values])[inverse]
        inherited = np.asarray([ROUTES.get(str(v), 'unresolved_context') for v in reader.frozen_labels(ds)['major']])
        recovered = (inherited == 'unresolved_context') & np.isin(contexts, list(VOCABULARY))
        vocab_path = root / 'outputs' / (ds + '_private_author_vocabulary.json')
        save(vocab_path, {'source': ds, 'native_author_field': field if field in fields else None,
            'original_label_vocabulary_not_individual_rows': [
                {'original_value': str(value), 'n_archive_rows': int(count), 'coarse_context': author_context(value)}
                for value, count in zip(values, counts)],
            'all_original_values_source_private_not_exportable_as_safe_summary': True})
        summaries = context_counts(contexts, support, meta['panel_order'])
        if sum(x['n_archive_rows'] for x in summaries) != n:
            raise ValueError('Author context denominator differs')
        reports.append({'source': ds, 'n_archive_rows': n, 'native_author_field': field if field in fields else None,
            'n_original_label_values': len(values), 'author_coarse_contexts': summaries,
            'n_original_unresolved_with_recognized_author_lineage': int(recovered.sum()),
            'recovered_context_counts_not_new_labels': {
                str(c): int((recovered & (contexts == c)).sum()) for c in np.unique(contexts[recovered])},
            'private_vocabulary_sha256': sha256_file(vocab_path),
            'n_author_labels_nonmissing': int((contexts != 'missing_or_explicit_unknown').sum()),
            'lineage_and_author_malignancy_not_inferred_from_each_other': True,
            'new_identity_assigned': False})
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source': ds, 'completed_sources': len(reports)})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
            raise ValueError('Author evidence output cap exceeded')
        reader._fields.clear(); del original, contexts, inherited, support; gc.collect()
    if sum(x['n_archive_rows'] for x in reports) != terminal['n_archive_rows']:
        raise ValueError('Author evidence archive denominator differs')
    path = root / 'outputs/native_author_label_evidence_v1.json'
    save(path, {'sources': reports, 'native_expression_read': False, 'new_authoritative_labels': False,
        'missing_value_vocabulary': sorted(MISSING),
        'merged_source_missing_placeholders_are_not_author_classes': True,
        'author_label_is_not_independent_identity_truth': True,
        'unreviewed_value_is_not_a_biologically_unknown_identity': True,
        'supplemental_cue_overlap_is_descriptive_not_validation': True})
    save(root / 'provenance/terminal_v1.json', {'status': 'NATIVE_AUTHOR_EVIDENCE_INVENTORY_COMPLETE_NOT_FINAL_ANNOTATION',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(reports),
        'output_sha256': sha256_file(path), 'native_expression_read': False, 'new_authoritative_labels': False,
        'wall_seconds': time.monotonic() - start,
        'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
