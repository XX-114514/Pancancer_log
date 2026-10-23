"""Audit full-archive evidence axes and independently sampled native values.

This verifies computation, scope and missingness. It is not biological truth
or a full second pass of reference predictions over every archived cell.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import time

import anndata as ad
from anndata._io.specs import read_elem
import celltypist
from celltypist.models import Model
import numpy as np
import pandas as pd

from annotation_representation_pilot_v1 import ROUTES, save
from audit_source_lineage_native_review_v1 import direct_rows
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_annotation_review_v1 import LINEAGE_PANELS, TYPE_PANELS, STATE_PANELS


def independent_panels(detected, genes):
    values = dict(zip(genes, detected.T))
    zero = np.zeros(len(detected), bool)
    def has(gene):
        return values.get(gene, zero)
    calls = {name: sum((has(g).astype(np.uint8) for g in panel)) >= 2
             for name, panel in {**LINEAGE_PANELS, **TYPE_PANELS}.items()}
    calls['B'] &= has('CD79A') | has('MS4A1')
    calls['Treg_related'] = has('FOXP3') & (has('IL2RA') | has('CTLA4'))
    calls['gamma_delta_T_related'] = has('TRDC') & (has('TRGC1') | has('TRGC2')) & has('CD3D')
    calls['CD8_T_related'] = has('CD8A') & has('CD8B') & (has('CD3D') | has('CD3E'))
    calls['Erythroid_related'] &= has('AHSP') | has('GYPA') | has('GYPB') | has('ALAS2') | has('SLC4A1')
    calls['Neutrophil_related'] &= has('CSF3R') | has('FCGR3B')
    return calls


def independent_families(calls):
    mapping = {'T_NK': ['T', 'NK_related'], 'B_Plasma': ['B', 'Plasma'],
        'Myeloid': ['Myeloid', 'DC_related', 'Neutrophil_related'], 'Mast': ['Mast'],
        'Endothelial': ['Endothelial', 'Lymphatic_related'], 'Fibroblast': ['Fibroblast'],
        'Mural': ['Mural'], 'Epithelial': ['Epithelial'], 'Glial': ['Glial_Oligodendrocyte_related'],
        'Astrocyte': ['Astrocyte_related'], 'Neuron': ['Neuron_related'], 'Erythroid': ['Erythroid_related']}
    return {name: np.logical_or.reduce([calls[p] for p in panels]) for name, panels in mapping.items()}


def unpack_named(bits, names):
    if len(names) > 32 or np.any(bits.astype(np.uint64) >> len(names)):
        raise ValueError('Out-of-range evidence bit')
    return {name: (bits & (1 << j)) != 0 for j, name in enumerate(names)}


def sampled_positions(side, meta):
    """Fixed coverage plus edge cases, independent of prediction correctness."""
    n = len(side['native_row'])
    positions = set(np.linspace(0, n - 1, min(48, n), dtype=np.int64).tolist())
    for key in ['positive_RNA', 'immune_reference_scope', 'organ_reference_scope']:
        for value in [False, True]:
            positions.update(np.flatnonzero(side[key] == value)[:4].tolist())
    for code in range(len(meta['archive_membership']['dictionary'])):
        positions.update(np.flatnonzero(side['archive_membership_code'] == code)[:4].tolist())
    for name in meta['source_reference_models']:
        positive = side[name + '__n_positive_labels']
        for mask in [positive < 0, positive == 0, positive == 1, positive > 1]:
            positions.update(np.flatnonzero(mask)[:4].tolist())
    return np.array(sorted(positions), np.int64)


def main(root, producer, completed_sources=None):
    start = time.monotonic()
    terminal = None
    if completed_sources is None:
        terminal = json.loads((producer / 'provenance/terminal_v1.json').read_text())
        supervision = json.loads((producer / 'provenance/supervisor_terminal_v1.json').read_text())
        if supervision['actual_child_exit_code'] != 0 or supervision['wall_timeout']:
            raise ValueError('Completed producer required')
    cfg = json.loads((producer / 'config/full_archive_reference_spec_v1.json').read_text())
    for entry in cfg['bindings'].values():
        if sha256_file(entry['path']) != entry['sha256']:
            raise ValueError('Frozen input changed')
    if completed_sources is None:
        index_path = producer / 'outputs/full_archive_reference_index_v1.json'
        if sha256_file(index_path) != terminal['index_sha256']:
            raise ValueError('Evidence index changed')
        index = json.loads(index_path.read_text())
    else:
        if not completed_sources or len(set(completed_sources)) != len(completed_sources):
            raise ValueError('Unique completed source list required')
        index = [json.loads((producer / 'outputs' / (ds + '_completion.json')).read_text()) for ds in completed_sources]
        if [r['source'] for r in index] != completed_sources:
            raise ValueError('Completed source snapshot differs')
        save(root / 'provenance/completed_source_snapshot_v1.json', index)
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    expected_sources = set(reader.archive_partitions) if completed_sources is None else set(completed_sources)
    if len(index) != len(expected_sources) or {r['source'] for r in index} != expected_sources or not expected_sources <= set(reader.archive_partitions):
        raise ValueError('Archive source scope differs')
    model_meta = json.loads((producer / 'outputs/reference_models_v1.json').read_text())
    reports = []; total_rows = direct_native_rows = official_predictions = 0
    max_error = 0.
    for record in index:
        ds = record['source']; rows = reader.rows(ds); n = len(rows)
        for field in ['arrays', 'metadata']:
            if sha256_file(producer / 'outputs' / record[field]) != record[field + '_sha256']:
                raise ValueError('Source evidence changed')
        meta = json.loads((producer / 'outputs' / record['metadata']).read_text())
        with np.load(producer / 'outputs' / record['arrays'], allow_pickle=False) as f:
            side = {name: f[name] for name in f.files}
        if any(len(value) != n for value in side.values()) or record['n_archive_rows'] != n:
            raise ValueError('Evidence row count differs')
        np.testing.assert_array_equal(rows, side['native_row'])
        if row_axis_sha256(rows) != meta['native_rows_sha256'] or text_axis_sha256(reader.native_ids(ds)) != meta['native_ids_sha256']:
            raise ValueError('Native identity axis differs')
        if not np.isfinite(side['native_total_RNA']).all() or np.any(side['native_total_RNA'] < 0):
            raise ValueError('Invalid RNA total')
        np.testing.assert_array_equal(side['positive_RNA'], side['native_total_RNA'] > 0)
        membership = np.asarray(meta['archive_membership']['dictionary'])[side['archive_membership_code']]
        np.testing.assert_array_equal(membership, reader.membership(ds)['membership_status'])
        np.testing.assert_array_equal(np.bincount(side['archive_membership_code'], minlength=len(meta['archive_membership']['dictionary'])), meta['archive_membership']['counts'])
        fields = set(reader.field_names(ds))
        for name, status in meta['carried_status_fields'].items():
            if name not in fields:
                if status['status'] != 'field_absent_not_false' or 'carried_status__' + name in side:
                    raise ValueError('Missing status fabricated')
                continue
            codes = side['carried_status__' + name]
            np.testing.assert_array_equal(np.asarray(status['dictionary'])[codes], reader.field(ds, name).astype(str))
            np.testing.assert_array_equal(np.bincount(codes, minlength=len(status['dictionary'])), status['counts'])
        expected_author = 'ev_author_label' if reader.is_legacy(ds) else (reader.primary.partitions[ds].get('annotation', {}).get('field') or 'source_author_label')
        if meta['author_field'] != (expected_author if expected_author in fields else None):
            raise ValueError('Author field origin differs')
        author = reader.field(ds, expected_author).astype(str) if expected_author in fields else np.full(n, '', dtype='<U1')
        if int((~np.isin(np.char.lower(author), ['', 'nan', 'none', 'na', 'unknown', 'unassigned'])).sum()) != meta['n_nonmissing_author_labels']:
            raise ValueError('Author label coverage differs')
        genes = meta['review_genes']; measured = np.asarray(meta['gene_measured'], bool)
        lookup = reader.feature_lookup(ds)
        np.testing.assert_array_equal(measured, [len(lookup.get(g, set())) == 1 for g in genes])
        inherited = np.asarray([ROUTES.get(str(v), 'unresolved_context') for v in reader.frozen_labels(ds)['major']])
        for low in range(0, n, 100000):
            sl = slice(low, min(n, low + 100000))
            packed = np.unpackbits(side['native_gene_detection_packbits'][sl], axis=1, bitorder='little')
            if packed[:, len(genes):].any() or packed[:, :len(genes)][:, ~measured].any():
                raise ValueError('Unmeasured or padding gene detection')
            calls = independent_panels(packed[:, :len(genes)].astype(bool), genes)
            stored = unpack_named(side['native_panel_bitmask'][sl], meta['native_panel_bit_order'])
            if set(stored) != set(calls):
                raise ValueError('Panel vocabulary differs')
            for name, value in calls.items():
                np.testing.assert_array_equal(stored[name], value)
            families = independent_families(calls)
            stored_families = unpack_named(side['native_family_bitmask'][sl], meta['native_family_bit_order'])
            if set(stored_families) != set(families):
                raise ValueError('Family vocabulary differs')
            for name, value in families.items():
                np.testing.assert_array_equal(stored_families[name], value)
            immune = np.logical_or.reduce([families[k] for k in ['T_NK', 'B_Plasma', 'Myeloid', 'Mast']])
            immune |= np.isin(inherited[sl], ['T_NK', 'B', 'Plasma', 'Myeloid', 'Neutrophil', 'Mast'])
            nonimmune = np.logical_or.reduce([families[k] for k in ['Endothelial', 'Fibroblast', 'Mural', 'Epithelial', 'Glial', 'Astrocyte', 'Neuron']])
            np.testing.assert_array_equal(side['immune_reference_scope'][sl], side['positive_RNA'][sl] & immune)
            np.testing.assert_array_equal(side['organ_reference_scope'][sl], side['positive_RNA'][sl] & nonimmune & (ds in cfg['source_organ_model']))
        for name in meta['source_reference_models']:
            unique, best, positive, probability = [side[name + '__' + key] for key in ['unique_code', 'best_code', 'n_positive_labels', 'maximum_probability']]
            scope = side['immune_reference_scope'] if name.startswith('immune_') else side['organ_reference_scope']
            evaluable = scope & (meta['reference_feature_coverage'][name]['feature_fraction'] >= cfg['minimum_model_feature_fraction'])
            np.testing.assert_array_equal(positive >= 0, evaluable)
            np.testing.assert_array_equal(np.isfinite(probability), evaluable)
            np.testing.assert_array_equal(unique >= 0, positive == 1)
            if np.any(best[~evaluable] != -1) or np.any(unique[~evaluable] != -1):
                raise ValueError('Out-of-scope reference prediction')
            if np.any(best[evaluable] < 0) or np.any(best >= len(model_meta[name]['labels'])) or np.any(unique >= len(model_meta[name]['labels'])):
                raise ValueError('Reference label bounds differ')
            if np.any((probability[evaluable] < 0) | (probability[evaluable] > 1)):
                raise ValueError('Reference probability outside bounds')
        chosen = sampled_positions(side, meta)
        z = reader.store(ds); key = meta['native_count_origin']
        raw = direct_rows(z[key], rows[chosen])
        totals = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
        np.testing.assert_array_equal(totals, side['native_total_RNA'][chosen])
        values = np.zeros((len(chosen), len(genes)), np.float64)
        for j, gene in enumerate(genes):
            if measured[j]:
                values[:, j] = raw[:, next(iter(lookup[gene]))].toarray().ravel()
        expected_detection = np.packbits(values > 0, axis=1, bitorder='little')
        np.testing.assert_array_equal(expected_detection, side['native_gene_detection_packbits'][chosen])
        scales = np.divide(10000., totals, out=np.zeros_like(totals), where=totals > 0)
        marker_log = np.log1p(values * scales[:, None])
        for j, state in enumerate(meta['state_order']):
            positions = [genes.index(g) for g in STATE_PANELS[state] if measured[genes.index(g)]]
            expected = np.full(len(chosen), np.nan)
            if len(positions) >= 3:
                expected = marker_log[:, positions].mean(axis=1); expected[totals == 0] = np.nan
            np.testing.assert_allclose(side['state_native_logCP10K'][chosen, j], expected, atol=2e-6, rtol=2e-6, equal_nan=True)
            valid = side['positive_RNA'] if len(positions) >= 3 else np.zeros(n, bool)
            np.testing.assert_array_equal(np.isfinite(side['state_native_logCP10K'][:, j]), valid)
        # Official API receives the complete native count denominator and an
        # independently constructed physical feature axis, including nonmodels.
        var = read_elem(z['var']); names = np.array(['NATIVE_NONMODEL_COLUMN_' + str(i) for i in range(raw.shape[1])], dtype=object)
        for j in range(len(var)):
            gene = str(var.iloc[j]['gene_symbol']) if 'gene_symbol' in var else str(var.index[j])
            if lookup.get(gene, set()) == {j}:
                names[j] = gene
        log = raw.astype(np.float64).multiply(scales[:, None]).tocsr(); log.data = np.log1p(log.data)
        source_predictions = 0
        for name in meta['source_reference_models']:
            if sha256_file(cfg['models'][name]['path']) != cfg['models'][name]['sha256']:
                raise ValueError('Reference file changed')
            ok = side[name + '__n_positive_labels'][chosen] >= 0
            if not ok.any():
                continue
            official = celltypist.annotate(ad.AnnData(log[ok], var=pd.DataFrame(index=names)),
                model=Model.load(cfg['models'][name]['path']), majority_voting=False, mode='prob match', p_thres=.5)
            probabilities = official.probability_matrix.to_numpy()
            positive = (probabilities > .5).sum(axis=1)
            best = official.decision_matrix.to_numpy().argmax(axis=1)
            for field, expected in [('best_code', best), ('n_positive_labels', positive), ('unique_code', np.where(positive == 1, best, -1))]:
                np.testing.assert_array_equal(side[name + '__' + field][chosen][ok], expected)
            actual = side[name + '__maximum_probability'][chosen][ok]
            np.testing.assert_allclose(actual, probabilities.max(axis=1), atol=3e-6, rtol=3e-6)
            max_error = max(max_error, float(np.max(np.abs(actual - probabilities.max(axis=1)))))
            source_predictions += int(ok.sum())
        counts = {'n_positive_RNA': int(side['positive_RNA'].sum()), 'n_immune_reference_scope': int(side['immune_reference_scope'].sum()),
            'n_organ_reference_scope': int(side['organ_reference_scope'].sum()),
            'n_cells_without_reference_scope': int((side['positive_RNA'] & ~side['immune_reference_scope'] & ~side['organ_reference_scope']).sum())}
        if any(meta[k] != v for k, v in counts.items()):
            raise ValueError('Source summary denominator differs')
        total_rows += n; direct_native_rows += len(chosen); official_predictions += source_predictions
        reports.append({'source': ds, 'n_archive_rows_checked': n, 'n_direct_native_rows': len(chosen),
                        'n_official_reference_predictions': source_predictions})
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source': ds, 'completed_sources': len(reports), 'completed_archive_rows': total_rows})
        del side, raw, log, values, marker_log, inherited, author
        reader._fields.clear(); gc.collect()
    if terminal is not None and total_rows != terminal['n_archive_rows']:
        raise ValueError('Archive total differs')
    status = ('ALL_SOURCE_AXES_SCOPES_AND_STRATIFIED_NATIVE_OFFICIAL_REFERENCE_REPLAY_ACCEPTED' if terminal is not None
              else 'COMPLETED_SOURCE_PREFLIGHT_REPLAY_ACCEPTED_NOT_FULL_ARCHIVE_ACCEPTANCE')
    output = {'status': status, 'full_archive_audit': terminal is not None,
        'source_checks': reports, 'n_archive_rows_checked': total_rows, 'n_direct_native_rows': direct_native_rows,
        'n_official_reference_predictions': official_predictions, 'max_reference_probability_error': max_error,
        'all_native_counts_and_predictions_recomputed': False, 'independent_biological_validation': False,
        'new_authoritative_labels': False}
    path = root / 'outputs/full_archive_evidence_audit_v1.json'; save(path, output)
    save(root / 'provenance/terminal_v1.json', {**{k: v for k, v in output.items() if k != 'source_checks'},
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'output_sha256': sha256_file(path),
        'producer_terminal_sha256': sha256_file(producer / 'provenance/terminal_v1.json') if terminal is not None else None,
        'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); parser.add_argument('--producer', required=True)
    parser.add_argument('--completed-source', action='append')
    args = parser.parse_args(); root = Path(args.output)
    try:
        main(root, Path(args.producer), args.completed_source)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
