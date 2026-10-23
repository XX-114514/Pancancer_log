"""Measure uncovered lineage/differentiation cues without assigning identity.

All selected archive rows are retained. Negative cues do not exclude a lineage;
positive cues do not imply malignancy, purity or a calibrated cell-type label.
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
from native_panel_aggregate_v1 import gather_csr_rows


PANELS = {'hepatocyte_related': ['ASGR1', 'HNF4A', 'ALB', 'TTR'],
    'melanocytic_differentiation_related': ['MITF', 'TYR', 'PMEL', 'MLANA'],
    'neural_neuroendocrine_differentiation_related': ['INSM1', 'ASCL1', 'NEUROD1', 'CHGA', 'CHGB', 'SYP']}
GENES = sorted({g for panel in PANELS.values() for g in panel} | {'AXL', 'NGFR'})


def cue_support(detected, measured, genes=GENES):
    if detected.ndim != 2 or detected.shape[1] != len(genes) or len(measured) != len(genes):
        raise ValueError('Supplementary marker axes differ')
    lookup = {g: i for i, g in enumerate(genes)}
    def has(g):
        return detected[:, lookup[g]] & bool(measured[lookup[g]])
    hepatocyte = (has('ASGR1') | has('HNF4A')) & (has('ALB') | has('TTR'))
    pigment = has('TYR').astype(int) + has('PMEL') + has('MLANA')
    melanocytic = (pigment + has('MITF') >= 2) & (pigment >= 1)
    secretory = has('CHGA') | has('CHGB') | has('SYP')
    neural_ne = (has('INSM1') | has('ASCL1') | has('NEUROD1')) & secretory
    return np.column_stack([hepatocyte, melanocytic, neural_ne])


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/supplemental_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Supplemental native input changed')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_supplemental_lineage_cues_v1.py')], capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Supplemental cue counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'n_tests': 4, 'actual_exit_code': 0, 'no_skips': True})
    reader = NativeArchiveCollection(cfg['archive_manifest']); index = []; total_nnz = 0; total_rows = 0
    for ds in cfg['source_order']:
        rows = reader.rows(ds); n = len(rows); lookup = reader.feature_lookup(ds)
        measured = np.array([len(lookup.get(g, set())) == 1 for g in GENES])
        columns = np.array([next(iter(lookup[g])) if measured[j] else 0 for j, g in enumerate(GENES)], np.int64)
        if len(set(columns[measured])) != int(measured.sum()):
            raise ValueError('Physical aliases in supplemental marker axis')
        arrays = {'native_row': rows, 'native_total_RNA': np.zeros(n, np.float64),
            'gene_detection_packbits': np.zeros((n, (len(GENES) + 7) // 8), np.uint8),
            'cue_support': np.zeros((n, len(PANELS)), bool),
            'panel_logCP10K': np.full((n, len(PANELS)), np.nan, np.float32)}
        z = reader.store(ds); key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        coverage = {name: [g for g in panel if measured[GENES.index(g)]] for name, panel in PANELS.items()}
        nnz = 0
        for low in range(0, n, cfg['read_block_rows']):
            sl = slice(low, min(n, low + cfg['read_block_rows']))
            raw = gather_csr_rows(z[key], rows[sl]); nnz += raw.nnz
            if not np.isfinite(raw.data).all() or (raw.data < 0).any() or not np.equal(raw.data, np.floor(raw.data)).all():
                raise ValueError('Native count domain changed')
            total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel(); positive = total > 0
            values = raw[:, columns].toarray(); values[:, ~measured] = 0
            detected = values > 0
            arrays['native_total_RNA'][sl] = total
            arrays['gene_detection_packbits'][sl] = np.packbits(detected, axis=1, bitorder='little')
            arrays['cue_support'][sl] = cue_support(detected, measured)
            scales = np.divide(10000., total, out=np.zeros_like(total), where=positive)
            log = np.log1p(values * scales[:, None])
            for j, (name, genes) in enumerate(coverage.items()):
                if len(genes) >= 3:
                    score = log[:, [GENES.index(g) for g in genes]].mean(axis=1); score[~positive] = np.nan
                    arrays['panel_logCP10K'][sl, j] = score
            if low == 0 or sl.stop == n or (low // cfg['read_block_rows']) % 64 == 0:
                save(root / 'checkpoints' / (str(time.time_ns()) + '.json'), {
                    'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'stage': 'supplemental_native_block',
                    'source': ds, 'n_source_rows_processed': sl.stop, 'completed_sources': len(index), 'completed_archive_rows': total_rows + sl.stop})
        output = root / 'outputs' / (ds + '_supplemental_native.npz')
        with output.open('xb') as f:
            np.savez_compressed(f, **arrays)
        with np.load(output, allow_pickle=False) as f:
            for name, value in arrays.items():
                np.testing.assert_array_equal(f[name], value)
        meta = {'source': ds, 'n_archive_rows': n, 'native_rows_sha256': row_axis_sha256(rows),
            'native_ids_sha256': text_axis_sha256(reader.native_ids(ds)), 'genes': GENES, 'measured': measured.tolist(),
            'panel_order': list(PANELS), 'panel_coverage': coverage, 'packbits_order': 'little',
            'cue_positive_counts': {name: int(arrays['cue_support'][:, j].sum()) for j, name in enumerate(PANELS)},
            'n_positive_RNA': int((arrays['native_total_RNA'] > 0).sum()), 'native_nnz': nnz,
            'negative_cue_excludes_lineage': False, 'positive_cue_assigns_identity': False,
            'neural_and_neuroendocrine_shared_expression_requires_context': True,
            'AXL_NGFR_detection_is_context_not_melanoma_identity': True, 'malignancy_inferred': False,
            'new_authoritative_labels': False, 'native_count_origin': key}
        meta_path = output.with_suffix('.json'); save(meta_path, meta)
        record = {'source': ds, 'n_archive_rows': n, 'arrays': output.name, 'arrays_sha256': sha256_file(output),
            'metadata': meta_path.name, 'metadata_sha256': sha256_file(meta_path)}
        save(root / 'outputs' / (ds + '_completion.json'), record); index.append(record)
        total_rows += n; total_nnz += nnz
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
            raise ValueError('Supplementary compact sidecar cap exceeded')
        del arrays, raw, values, detected, log
        reader._fields.clear(); gc.collect()
    if set(cfg['source_order']) != set(reader.archive_partitions) or total_rows != sum(p['n_archive_rows'] for p in reader.archive_partitions.values()):
        raise ValueError('Supplemental archive scope incomplete')
    path = root / 'outputs/supplemental_native_index_v1.json'; save(path, index)
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_ARCHIVE_SUPPLEMENTAL_NATIVE_CUES_COMPLETE_AUDIT_PENDING',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'index_sha256': sha256_file(path),
        'n_sources': len(index), 'n_archive_rows': total_rows, 'n_native_nnz_read': total_nnz,
        'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'new_authoritative_labels': False, 'native_counts_copied': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
