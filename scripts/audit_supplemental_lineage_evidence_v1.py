"""Replay supplemental cue axes, logic and sampled native count measurements."""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import time

import numpy as np

from annotation_representation_pilot_v1 import save
from audit_source_lineage_native_review_v1 import direct_rows
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_supplemental_lineage_evidence_v1 import PANELS, GENES


def main(root, producer):
    start = time.monotonic()
    supervision = json.loads((producer / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervision['actual_child_exit_code'] != 0 or supervision['wall_timeout']:
        raise ValueError('Supplemental producer completion required')
    terminal = json.loads((producer / 'provenance/terminal_v1.json').read_text())
    cfg = json.loads((producer / 'config/supplemental_spec_v1.json').read_text())
    for entry in cfg['bindings'].values():
        if sha256_file(entry['path']) != entry['sha256']:
            raise ValueError('Supplemental input changed')
    index_path = producer / 'outputs/supplemental_native_index_v1.json'
    if sha256_file(index_path) != terminal['index_sha256']:
        raise ValueError('Supplemental index changed')
    index = json.loads(index_path.read_text()); reader = NativeArchiveCollection(cfg['archive_manifest'])
    if len(index) != len(reader.archive_partitions) or {r['source'] for r in index} != set(reader.archive_partitions):
        raise ValueError('Full native source scope differs')
    summaries = []; total = replayed = 0; largest_error = 0.
    for record in index:
        ds = record['source']; rows = reader.rows(ds); n = len(rows)
        for field in ['arrays', 'metadata']:
            if sha256_file(producer / 'outputs' / record[field]) != record[field + '_sha256']:
                raise ValueError('Supplemental artifact changed')
        meta = json.loads((producer / 'outputs' / record['metadata']).read_text())
        with np.load(producer / 'outputs' / record['arrays'], allow_pickle=False) as f:
            side = {key: f[key] for key in f.files}
        if any(len(v) != n for v in side.values()) or meta['n_archive_rows'] != n:
            raise ValueError('Supplemental row denominator differs')
        np.testing.assert_array_equal(side['native_row'], rows)
        if row_axis_sha256(rows) != meta['native_rows_sha256'] or text_axis_sha256(reader.native_ids(ds)) != meta['native_ids_sha256']:
            raise ValueError('Supplemental native identity differs')
        if meta['genes'] != GENES or meta['panel_order'] != list(PANELS):
            raise ValueError('Frozen panel vocabulary differs')
        lookup = reader.feature_lookup(ds); measured = np.array([len(lookup.get(g, set())) == 1 for g in GENES])
        np.testing.assert_array_equal(measured, meta['measured'])
        counts = np.zeros(3, np.int64)
        for low in range(0, n, 100000):
            sl = slice(low, min(n, low + 100000))
            detection = np.unpackbits(side['gene_detection_packbits'][sl], axis=1, bitorder='little')
            if detection[:, len(GENES):].any() or detection[:, :len(GENES)][:, ~measured].any():
                raise ValueError('Missing or padding gene became positive')
            d = dict(zip(GENES, detection[:, :len(GENES)].astype(bool).T))
            liver = (d['ASGR1'] | d['HNF4A']) & (d['ALB'] | d['TTR'])
            pigment = d['TYR'].astype(np.int64) + d['PMEL'] + d['MLANA']
            melan = ((pigment + d['MITF']) >= 2) & (pigment > 0)
            neural = (d['INSM1'] | d['ASCL1'] | d['NEUROD1']) & (d['CHGA'] | d['CHGB'] | d['SYP'])
            expected = np.column_stack([liver, melan, neural])
            np.testing.assert_array_equal(expected, side['cue_support'][sl]); counts += expected.sum(axis=0)
        for j, panel in enumerate(PANELS):
            if int(counts[j]) != meta['cue_positive_counts'][panel]:
                raise ValueError('Cue summary count differs')
            available = [g for g in PANELS[panel] if measured[GENES.index(g)]]
            if available != meta['panel_coverage'][panel]:
                raise ValueError('Panel measurement coverage differs')
            finite = (side['native_total_RNA'] > 0) if len(available) >= 3 else np.zeros(n, bool)
            np.testing.assert_array_equal(finite, np.isfinite(side['panel_logCP10K'][:, j]))
        chosen = set(np.linspace(0, n - 1, min(n, 48), dtype=np.int64).tolist())
        for j in range(3):
            for value in [False, True]:
                chosen.update(np.flatnonzero(side['cue_support'][:, j] == value)[:4].tolist())
        chosen.update(np.flatnonzero(side['native_total_RNA'] == 0)[:4].tolist()); chosen = np.array(sorted(chosen))
        raw = direct_rows(reader.store(ds)[meta['native_count_origin']], rows[chosen])
        totals = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
        np.testing.assert_array_equal(totals, side['native_total_RNA'][chosen])
        values = np.zeros((len(chosen), len(GENES)), np.float64)
        for j, gene in enumerate(GENES):
            if measured[j]:
                values[:, j] = raw[:, next(iter(lookup[gene]))].toarray().ravel()
        np.testing.assert_array_equal(np.packbits(values > 0, axis=1, bitorder='little'), side['gene_detection_packbits'][chosen])
        scales = np.divide(10000., totals, out=np.zeros_like(totals), where=totals > 0)
        log = np.log1p(values * scales[:, None])
        for j, panel in enumerate(PANELS):
            cols = [GENES.index(g) for g in PANELS[panel] if measured[GENES.index(g)]]
            expected = np.full(len(chosen), np.nan)
            if len(cols) >= 3:
                expected = log[:, cols].mean(axis=1); expected[totals == 0] = np.nan
            actual = side['panel_logCP10K'][chosen, j]
            np.testing.assert_allclose(actual, expected, atol=2e-6, rtol=2e-6, equal_nan=True)
            finite = np.isfinite(expected)
            if finite.any():
                largest_error = max(largest_error, float(np.max(np.abs(actual[finite] - expected[finite]))))
        total += n; replayed += len(chosen)
        summaries.append({'source': ds, 'n_archive_rows': n, 'n_direct_native_rows': len(chosen)})
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source': ds, 'completed_sources': len(summaries), 'completed_archive_rows': total})
        del side, raw, values, log
        reader._fields.clear(); gc.collect()
    if total != terminal['n_archive_rows']:
        raise ValueError('Supplemental total scope differs')
    report = {'status': 'SUPPLEMENTAL_ALL_AXES_CUE_LOGIC_AND_STRATIFIED_NATIVE_REPLAY_ACCEPTED',
        'source_checks': summaries, 'n_archive_rows': total, 'n_direct_native_rows': replayed,
        'max_native_panel_score_error': largest_error, 'biological_identity_truth_verified': False,
        'full_native_counts_second_pass': False, 'new_authoritative_labels': False}
    path = root / 'outputs/supplemental_native_audit_v1.json'; save(path, report)
    save(root / 'provenance/terminal_v1.json', {**{k: v for k, v in report.items() if k != 'source_checks'},
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'output_sha256': sha256_file(path),
        'producer_terminal_sha256': sha256_file(producer / 'provenance/terminal_v1.json'),
        'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); parser.add_argument('--producer', required=True)
    args = parser.parse_args(); root = Path(args.output)
    try:
        main(root, Path(args.producer))
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
