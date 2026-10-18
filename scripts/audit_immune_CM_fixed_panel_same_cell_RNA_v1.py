"""Reconstruct original T membership and independently replay all selected RNA."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 8 * 1024**2)
    run = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    terminal = json.loads((run / 'private_terminal.json').read_text())
    assert terminal['source_sha256'] == '920fad637b6c583e0eac04662b72cbe20e9de4fbb17ae756c3cac18bbc2514d0'
    for n, v in terminal['files'].items(): assert sha(run / n) == v['sha256']
    plan = json.loads((run / 'plan_before_native_same_cell_RNA_read.json').read_text())
    out = root / 'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    (out / 'plan_before_full_original_T_RNA_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        producer_terminal_sha256=sha(run / 'private_terminal.json'),
        scope='All original T cells, every34 gene and full-native RNA/depth. Independently rebuild cells/units from original obs and refined catalog, source physical maps from original var/masks and aliases, CM/pre membership from original metadata. No sampling, cell-ID export, new labels or clinical model.'), indent=2) + '\n')
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    assert sha(atlas / 'manifests/immune_authority_release_public.json') == plan['atlas_release_sha256']
    catalog_file = atlas / 'tables/subcluster_catalog_refined_private.tsv'
    assert sha(catalog_file) == plan['refined_catalog_sha256']
    catalog = pd.read_csv(catalog_file, sep='\t', keep_default_na=False)
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    metadata = {k: np.asarray(read_elem(z['obs'][k])) for k in ['dataset', 'library_id', 'node_id', 'major']}
    use = np.isin(metadata['node_id'], catalog.loc[catalog.major.eq('T'), 'node_id'].to_numpy())
    assert np.array_equal(use, metadata['major'] == 'T')
    original_rows = np.flatnonzero(use)
    table = pd.DataFrame({k: v[use] for k, v in metadata.items()})
    units = pd.read_csv(run / 'private_original_T_unit_scope.tsv', sep='\t', keep_default_na=False)
    keys = ['dataset', 'library_id', 'node_id', 'major']
    independent = table.groupby(keys, sort=True).size().reset_index(name='n_cells')
    assert independent[keys].equals(units[keys]) and np.array_equal(independent.n_cells, units.n_cells)
    lookup = pd.MultiIndex.from_frame(independent[keys])
    unit_index = lookup.get_indexer(pd.MultiIndex.from_frame(table[keys])).astype(np.int32)
    assert (unit_index >= 0).all()
    original_meta = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    assert sha(original_meta) == plan['original_CM_material_metadata_sha256']
    m = pd.read_csv(original_meta, sep='\t', keep_default_na=False)
    sample_keys = independent.dataset.astype(str) + '::' + independent.library_id.astype(str)
    assert np.array_equal(sample_keys.isin(m.sampleID), units.original_CM_scope_library)
    assert np.array_equal(sample_keys.isin(m.loc[m.timepoint == 'pre', 'sampleID']), units.original_CM_scope_pre_library)
    cache = np.load(run / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False)
    assert np.array_equal(original_rows, cache['global_original_T_row_index'])
    assert np.array_equal(unit_index, cache['original_T_unit_index'])
    panel = cache['raw_gene_counts']
    genes = cache['genes'].tolist()
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed_path) == plan['fixed_panel_contract_sha256']
    fixed = json.loads(fixed_path.read_text())
    assert genes == sorted(set(g for gg in fixed['panels'].values() for g in gg))
    source_order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    symbol = np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    aliases = np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    alias_sets = [set(a.split('|')) if st == 'stable_Ensembl_ID' else set() for a, st in zip(aliases, identity)]
    ambiguous = np.array([c.startswith('AMBIGUOUS::') or st == 'ambiguous_source_feature' for c, st in zip(canonical, identity)])
    native_maps = []
    for si, source in enumerate(source_order):
        mapping = []
        for gene in genes:
            candidates = np.flatnonzero(measured[:, si] & ~ambiguous & ((symbol == gene) | (canonical == gene) |
                np.array([gene in v for v in alias_sets])))
            mapping.append(int(candidates[0]) if len(candidates) == 1 else -1)
        collision = {v for v in mapping if v >= 0 and mapping.count(v) > 1}
        native_maps.append(np.array([-1 if v in collision else v for v in mapping], np.int64))
    source_contract = json.loads((run / 'private_source_gene_measurement_contract.json').read_text())
    assert len(source_contract) == len(source_order)
    for si, source in enumerate(source_order):
        assert source_contract[si]['dataset'] == source
        assert np.array_equal(native_maps[si], source_contract[si]['physical_columns'])
    source_by_unit = np.array([source_order.index(v) for v in independent.dataset])
    source_by_cell = source_by_unit[unit_index]
    c = z['layers/counts'];nc, nf = map(int, c.attrs['shape'])
    expected_total = cache['full_native_total_counts']
    expected_detected = cache['full_native_detected_genes']
    done = 0
    for start in range(0, nc, 16384):
        stop = min(nc, start + 16384)
        mask = use[start:stop]
        if mask.any():
            ptr = np.asarray(c['indptr'][start:stop + 1], np.int64)
            x = csr_matrix((np.asarray(c['data'][int(ptr[0]):int(ptr[-1])]), np.asarray(c['indices'][int(ptr[0]):int(ptr[-1])]), ptr - ptr[0]), shape=(stop - start, nf))[mask]
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert (x.data > 0).all() and np.isfinite(x.data).all() and (x.data == np.floor(x.data)).all()
            count = x.shape[0];s = source_by_cell[done:done + count]
            assert measured[x.indices, np.repeat(s, np.diff(x.indptr))].all()
            totals = np.add.reduceat(x.data.astype(np.int64), x.indptr[:-1])
            assert np.array_equal(totals, expected_total[done:done + count])
            assert np.array_equal(np.diff(x.indptr), expected_detected[done:done + count])
            raw = np.full((count, len(genes)), -1, np.int32)
            for si in np.unique(s):
                selected = np.flatnonzero(s == si);present = np.flatnonzero(native_maps[si] >= 0)
                raw[np.ix_(selected, present)] = x[selected][:, native_maps[si][present]].toarray()
            assert np.array_equal(raw, panel[done:done + count])
            done += count
        if start // 16384 % 16 == 0 or stop == nc:
            (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), original_rows_replayed=stop, T_cells_independently_replayed=done)) + '\n')
    assert done == plan['original_T_cells'] == len(original_rows)
    assert sha(atlas / 'manifests/immune_authority_release_public.json') == plan['atlas_release_sha256']
    result = dict(status='ACCEPTED_ALL_ORIGINAL_T_CELLS_FIXED34_RNA_FULL_DEPTH_AND_SCOPE_REPLAY',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), T_cells=done,
        gene_values_replayed=int(panel.size), original_T_units=len(independent),
        producer_terminal_sha256=sha(run / 'private_terminal.json'), source_native_measurability_and_zero_distinction=True,
        original_obs_unit_and_CM_pre_membership_independently_rebuilt=True,
        input='Same accepted immutable authority object independently re-read; not independent biological collection.',
        no_new_state_or_clinical_or_cell_interaction_claim=True)
    (out / 'private_independent_acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    result['files'] = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], completed_UTC=result['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
