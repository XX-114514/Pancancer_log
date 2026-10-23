"""Route all audited archive rows for source-local discovery, preserving conflicts.

The frozen pilot policy is reused with no invented local group evidence. These
are provisional analysis branches, not final identities or model eligibility.
Supplemental differentiation cues and original labels remain separate evidence.
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

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_lineage_routing_v1 import FAMILIES, STATUS, reference_family, route_cells
from summarize_full_native_annotation_evidence_v1 import accepted, read_source


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/full_native_routing_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Frozen full-native routing input changed')
    refs, rt = accepted(Path(cfg['reference']), Path(cfg['reference_audit']), 'full_archive_reference_index_v1.json')
    cues, ct = accepted(Path(cfg['supplemental']), Path(cfg['supplemental_audit']), 'supplemental_native_index_v1.json')
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    if set(refs) != set(cues) or set(refs) != set(reader.archive_partitions) or rt['n_archive_rows'] != ct['n_archive_rows']:
        raise ValueError('Accepted evidence source or row scope differs')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_native_lineage_routing_v1.py')],
        capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 6 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Frozen provisional routing counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 6, 'no_skips': True})
    reference_root = Path(cfg['reference'])
    models_path = reference_root / 'outputs/reference_models_v1.json'
    model_meta = json.loads(models_path.read_text())
    maps = {name: np.asarray([FAMILIES.index(family) if family is not None else -1
        for family in [reference_family(name, label) for label in model['labels']]], np.int16)
        for name, model in model_meta.items()}
    save(root / 'config/routing_dictionary_v1.json', {'families': FAMILIES, 'reasons': STATUS,
        'no_local_group_fit_available': True, 'author_labels_not_substituted_for_local_group_evidence': True,
        'supplemental_cues_do_not_override_coarse_family_conflicts': True,
        'reference_models_sha256': sha256_file(models_path)})
    index = []
    for ds in sorted(refs):
        rm, rv = read_source(reference_root, refs[ds])
        cm, cv = read_source(Path(cfg['supplemental']), cues[ds])
        rows = reader.rows(ds); n = len(rows)
        np.testing.assert_array_equal(rows, rv['native_row'])
        np.testing.assert_array_equal(rows, cv['native_row'])
        np.testing.assert_array_equal(rv['native_total_RNA'], cv['native_total_RNA'])
        identity_sha = text_axis_sha256(reader.native_ids(ds))
        if rm['native_ids_sha256'] != identity_sha or cm['native_ids_sha256'] != identity_sha:
            raise ValueError('Provisional routing native identity axes differ')
        if rm['native_family_bit_order'] != FAMILIES:
            raise ValueError('Native family vocabulary order differs')
        native = ((rv['native_family_bitmask'][:, None] >> np.arange(len(FAMILIES))) & 1).astype(bool)
        reference = np.zeros_like(native)
        for name in rm['source_reference_models']:
            codes = rv[name + '__unique_code']; count = rv[name + '__n_positive_labels']
            np.testing.assert_array_equal(codes >= 0, count == 1)
            if np.any(codes >= len(maps[name])):
                raise ValueError('Reference class index outside frozen vocabulary')
            valid = codes >= 0
            families = np.full(n, -1, np.int16); families[valid] = maps[name][codes[valid]]
            for j in range(len(FAMILIES)):
                reference[:, j] |= families == j
        # No full-source group has been fitted at this stage. Old author labels
        # and old immune-only labels cannot be used as newly fitted local groups.
        route, reason = route_cells(native, reference, np.full(n, -1, np.int16), rv['positive_RNA'])
        native_n, reference_n = native.sum(axis=1), reference.sum(axis=1)
        conflict = (native_n > 1) | (reference_n > 1)
        conflict |= (native_n == 1) & (reference_n == 1) & (native.argmax(axis=1) != reference.argmax(axis=1))
        if np.any(route[conflict | ~rv['positive_RNA']] >= 0):
            raise ValueError('Conflict or zero-RNA record was routed')
        if np.any((route >= 0) & (native_n != 1)):
            raise ValueError('Reference-only or author-only route invented without local group evidence')
        values = {'native_row': rows, 'provisional_family_code': route, 'routing_reason_code': reason,
            'native_family_bitmask': rv['native_family_bitmask'],
            'reference_family_bitmask': (reference.astype(np.uint32) * (1 << np.arange(len(FAMILIES), dtype=np.uint32))).sum(axis=1).astype(np.uint32),
            'supplemental_cue_bitmask': (cv['cue_support'].astype(np.uint8) * (1 << np.arange(len(cm['panel_order']), dtype=np.uint8))).sum(axis=1).astype(np.uint8),
            'positive_RNA': rv['positive_RNA']}
        output = root / 'outputs' / (ds + '_provisional_routes_v1.npz')
        with output.open('xb') as f:
            np.savez_compressed(f, **values)
        with np.load(output, allow_pickle=False) as saved:
            for name, value in values.items():
                np.testing.assert_array_equal(saved[name], value)
        meta = {'source': ds, 'n_archive_rows': n, 'native_rows_sha256': row_axis_sha256(rows),
            'native_ids_sha256': identity_sha, 'family_code_dictionary': FAMILIES, 'unresolved_code': -1,
            'routing_reason_dictionary': STATUS, 'supplemental_cue_bit_order': cm['panel_order'],
            'provisional_family_counts': {('unresolved' if j < 0 else FAMILIES[j]): int((route == j).sum()) for j in np.unique(route)},
            'routing_reason_counts': {STATUS[j]: int((reason == j).sum()) for j in np.unique(reason)},
            'n_unresolved_with_supplemental_cue': int(((route < 0) & cv['cue_support'].any(axis=1)).sum()),
            'reference_evidence_arrays_sha256': refs[ds]['arrays_sha256'],
            'supplemental_evidence_arrays_sha256': cues[ds]['arrays_sha256'],
            'local_group_fit_available': False, 'native_expression_read': False,
            'author_and_original_labels_unchanged_in_source': True,
            'model_composition_and_clinical_eligibility_assigned': False,
            'new_authoritative_labels': False}
        if sum(meta['provisional_family_counts'].values()) != n or sum(meta['routing_reason_counts'].values()) != n:
            raise ValueError('Routing summary denominator differs')
        meta_path = output.with_suffix('.json'); save(meta_path, meta)
        index.append({'source': ds, 'n_archive_rows': n, 'arrays': output.name, 'arrays_sha256': sha256_file(output),
            'metadata': meta_path.name, 'metadata_sha256': sha256_file(meta_path)})
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source': ds, 'completed_sources': len(index)})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
            raise ValueError('Full-native provisional routing output cap exceeded')
        reader._fields.clear(); del rv, cv, native, reference, route, reason, values; gc.collect()
    if sum(r['n_archive_rows'] for r in index) != rt['n_archive_rows']:
        raise ValueError('Full-native routing archive denominator differs')
    path = root / 'outputs/full_native_provisional_routing_index_v1.json'; save(path, index)
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_NATIVE_PROVISIONAL_ROUTING_COMPLETE_SOURCE_LOCAL_DISCOVERY_PENDING',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(index),
        'index_sha256': sha256_file(path), 'wall_seconds': time.monotonic() - start,
        'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'native_expression_read': False, 'new_authoritative_labels': False,
        'model_and_downstream_eligibility_assigned': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
