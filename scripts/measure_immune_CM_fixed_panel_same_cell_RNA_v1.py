"""Measure all fixed clinical panels on original native T-cell RNA.

Full source-native depth; exact original cells, units and material membership.
No new cell labels, clinical refit, statistical effects or whole-gene matrix copy.
"""
import json
import os
from collections import defaultdict
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
    budget = resource(root, 128 * 1024**2)
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    scope = root / 'existing_immune_fine_native_measurement_scope_20261016_v1'
    terminal = json.loads((scope / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(scope / name) == digest
    contract = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed = json.loads(contract.read_text())
    genes = sorted(set(g for panel in fixed['panels'].values() for g in panel))
    assert len(genes) == 34 and len(fixed['panels']) == 7
    units = pd.read_csv(scope / 'private_all_frozen_T_myeloid_node_units.tsv', sep='\t', keep_default_na=False)
    assert np.array_equal(units.unit_index, np.arange(len(units)))
    mapping = np.load(scope / 'private_original_cell_to_fine_node_unit.npy', mmap_mode='r')
    t_units = units[units.major.eq('T')].copy().reset_index(drop=True)
    convert = np.full(len(units), -1, np.int32)
    convert[t_units.unit_index.to_numpy()] = np.arange(len(t_units))
    cell_unit = np.full(len(mapping), -1, np.int32)
    use = mapping >= 0
    cell_unit[use] = convert[mapping[use]]
    rows = np.flatnonzero(cell_unit >= 0)
    assert len(rows) == int(t_units.n_cells.sum())
    assert np.array_equal(np.bincount(cell_unit[rows], minlength=len(t_units)), t_units.n_cells)
    original_meta = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    meta = pd.read_csv(original_meta, sep='\t', keep_default_na=False)
    sample_keys = t_units.dataset.astype(str) + '::' + t_units.library_id.astype(str)
    assert meta.sampleID.is_unique and set(meta.sampleID) <= set(sample_keys)
    assert (meta.n_original_libraries == 1).all()
    t_units['original_CM_scope_library'] = sample_keys.isin(meta.sampleID)
    t_units['original_CM_scope_pre_library'] = sample_keys.isin(meta.loc[meta.timepoint.eq('pre'), 'sampleID'])
    catalog = pd.read_csv(atlas / 'tables/subcluster_catalog_refined_private.tsv', sep='\t', keep_default_na=False)
    assert set(t_units.node_id) == set(catalog.loc[catalog.major.eq('T'), 'node_id'])
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    counts = z['layers/counts']
    nc, nf = map(int, counts.attrs['shape'])
    assert len(mapping) == nc
    source_order = np.asarray(read_elem(z['uns/source_order'])).astype(str).tolist()
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    symbols = np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    status = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    aliases = np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    maps = []
    source_contract = []
    for si, source in enumerate(source_order):
        names = defaultdict(set)
        for col in np.flatnonzero(measured[:, si]):
            if canonical[col].startswith('AMBIGUOUS::') or status[col] == 'ambiguous_source_feature':
                continue
            labels = {symbols[col], canonical[col]}
            if status[col] == 'stable_Ensembl_ID': labels.update(aliases[col].split('|'))
            for label in labels - {'', 'nan', 'None', '__SOURCE_VALUE_MISSING__'}:
                names[label].add(int(col))
        columns = [next(iter(names[g])) if len(names[g]) == 1 else -1 for g in genes]
        # Physical collisions across requested symbols remain unavailable.
        collisions = {col for col in columns if col >= 0 and columns.count(col) > 1}
        columns = np.array([-1 if col in collisions else col for col in columns], np.int64)
        maps.append(columns)
        source_contract.append(dict(dataset=source, requested_genes=genes, physical_columns=columns.tolist(),
            complete_panels={k: bool(all(columns[genes.index(g)] >= 0 for g in gg)) for k, gg in fixed['panels'].items()}))
    source_by_unit = np.array([source_order.index(v) for v in t_units.dataset], np.int16)
    out = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        input_scope_terminal_sha256=sha(scope / 'private_terminal.json'),
        atlas_release_sha256=sha(atlas / 'manifests/immune_authority_release_public.json'),
        refined_catalog_sha256=sha(atlas / 'tables/subcluster_catalog_refined_private.tsv'),
        fixed_panel_contract_sha256=sha(contract), original_CM_material_metadata_sha256=sha(original_meta),
        original_T_cells=len(rows), original_T_units=len(t_units), original_T_nodes=int(t_units.node_id.nunique()),
        original_CM_scope_libraries=len(meta), original_CM_scope_T_cells=int(t_units.loc[t_units.original_CM_scope_library, 'n_cells'].sum()),
        original_CM_scope_pre_libraries=int(meta.timepoint.eq('pre').sum()),
        original_CM_scope_pre_T_cells=int(t_units.loc[t_units.original_CM_scope_pre_library, 'n_cells'].sum()),
        source_CSR_stored_entries=int(counts['data'].shape[0]),
        measurement='Read original layers/counts, canonicalize working row blocks; all34 original fixed genes with exact unique source-native physical mapping. Unmeasured/unresolved=-1, measured absent=0. Full native RNA and detected genes from all measured features, never panel-only depth.',
        scopes='All original T cells; exact original CM library membership; pre only when exact original CM timepoint==pre. Unknown timepoints remain unknown. These are overlapping descriptive subsets, not independent cohorts.',
        prospective_description='All7 unchanged panels, all21 panel pairs, original16 nodes and all source/library units; detection thresholds1/2, >=2 distinct panel genes; compare all T versus pan-T>=2-gene capture. FOXP3-supported regulatory capture requires FOXP3 plus>=1 other original regulatory gene. Missing FOXP3 capture is not absence of regulatory biology.',
        interpretation='Descriptive source/measurement attribution after clinical/spatial results, not new clinical effects, DEG, learned expressionMCP or formal Treg/state annotation. No clonotype values or patient-level output exported.',
        output_cap_bytes=128 * 1024**2, no_original_CM_or_annotation_or_clinical_changes=True)
    (out / 'plan_before_native_same_cell_RNA_read.json').write_text(json.dumps(plan, indent=2) + '\n')
    (out / 'private_source_gene_measurement_contract.json').write_text(json.dumps(source_contract, indent=2) + '\n')
    t_units.to_csv(out / 'private_original_T_unit_scope.tsv', sep='\t', index=False)
    panel = np.full((len(rows), len(genes)), -1, np.int32)
    total = np.zeros(len(rows), np.int64)
    detected = np.zeros(len(rows), np.int32)
    selected_units = np.zeros(len(rows), np.int32)
    offset = 0
    for start in range(0, nc, 8192):
        stop = min(nc, start + 8192)
        local = cell_unit[start:stop] >= 0
        if local.any():
            ptr = np.asarray(counts['indptr'][start:stop + 1], np.int64)
            lo, hi = int(ptr[0]), int(ptr[-1])
            x = csr_matrix((np.asarray(counts['data'][lo:hi]), np.asarray(counts['indices'][lo:hi]), ptr - lo), shape=(stop - start, nf))[local].copy()
            x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
            assert np.isfinite(x.data).all() and (x.data > 0).all() and np.equal(x.data, np.floor(x.data)).all()
            ids = cell_unit[start:stop][local]
            sources = source_by_unit[ids]
            assert measured[x.indices, np.repeat(sources, np.diff(x.indptr))].all()
            ns = len(ids)
            assert np.array_equal(rows[offset:offset + ns], start + np.flatnonzero(local))
            for si in np.unique(sources):
                take = np.flatnonzero(sources == si)
                present = np.flatnonzero(maps[si] >= 0)
                v = x[take][:, maps[si][present]].toarray()
                assert v.max(initial=0) <= np.iinfo(np.int32).max
                panel[offset + take[:, None], present[None, :]] = v.astype(np.int32)
            total[offset:offset + ns] = np.asarray(x.sum(axis=1)).ravel().astype(np.int64)
            detected[offset:offset + ns] = np.diff(x.indptr)
            selected_units[offset:offset + ns] = ids
            offset += ns
        if start // 8192 % 16 == 0 or stop == nc:
            (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), original_rows_read=stop, original_rows=nc, T_cells_measured=offset, original_T_cells=len(rows))) + '\n')
    assert offset == len(rows) and (total > 0).all()
    assert np.array_equal(np.bincount(selected_units, minlength=len(t_units)), t_units.n_cells)
    for si in np.unique(source_by_unit):
        take = source_by_unit[selected_units] == si
        assert np.array_equal((panel[take] >= 0).all(axis=0), maps[si] >= 0)
    np.savez_compressed(out / 'private_fixed34_gene_T_counts_and_depth.npz',
        global_original_T_row_index=rows, original_T_unit_index=selected_units, genes=np.asarray(genes),
        raw_gene_counts=panel, full_native_total_counts=total, full_native_detected_genes=detected)
    assert sha(original_meta) == plan['original_CM_material_metadata_sha256']
    assert sha(atlas / 'manifests/immune_authority_release_public.json') == plan['atlas_release_sha256']
    assert sha(contract) == plan['fixed_panel_contract_sha256']
    files = [p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files) < plan['output_cap_bytes']
    result = dict(status='COMPLETE_ALL_ORIGINAL_T_CELLS_FIXED34_GENE_RNA_MEASUREMENT_RAW_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        T_cells=len(rows), T_units=len(t_units), genes=len(genes),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
        independent_raw_replay=False, clinical_or_biological_interaction_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], completed_UTC=result['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
