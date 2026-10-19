"""Read-only cache/reference coverage; no expression scores or clinical fitting."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024**2), b''): h.update(chunk)
    return h.hexdigest()


def label_T(raw, genes, panels, anchor, effector, threshold):
    def captured(gg):
        return np.count_nonzero(raw[:, [genes.index(g) for g in gg]] >= threshold, axis=1) >= 2
    core = captured(panels['pan_T'])
    reg = captured(panels['regulatory_associated'])
    anchored = reg & (raw[:, genes.index(anchor)] >= threshold)
    eff = captured(effector)
    labels = np.full(len(raw), 'T_neither_capture', dtype='U40')
    labels[eff & ~reg] = 'T_effector_only'
    labels[reg & ~anchored & ~eff] = 'T_reg_without_anchor_only'
    labels[reg & ~anchored & eff] = 'T_reg_without_anchor_effector'
    labels[anchored & ~eff] = 'T_anchored_reg_only'
    labels[anchored & eff] = 'T_anchored_reg_effector'
    labels[~core] = 'T_core_not_captured'
    assert len(labels) == len(raw)
    return labels


def main():
    cfg = json.load(sys.stdin)
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed_path) == cfg['fixed_panel_contract_SHA256']
    panels = json.loads(fixed_path.read_text())['panels']
    old = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    manifest = json.loads((old / 'private_terminal.json').read_text())
    for name in ['private_fixed34_gene_T_counts_and_depth.npz', 'private_original_T_unit_scope.tsv', 'private_source_gene_measurement_contract.json']:
        assert sha(old / name) == manifest['files'][name]['sha256']
    cache = dict(np.load(old / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False))
    genes = cache['genes'].tolist()
    assert len(genes) == 34
    mcp = root / 'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    audit_path = root / 'immune_RNA_MCP_native_expression_input_audit_20261018_v1/private_terminal.json'
    audit = json.loads(audit_path.read_text())
    assert audit['producer_terminal_sha256'] == sha(mcp / 'private_terminal.json')
    mm = json.loads((mcp / 'private_terminal.json').read_text())
    name = 'private_frozen_original_cell_selection.tsv'
    assert sha(mcp / name) == mm['files'][name]['sha256']
    v = pd.read_csv(mcp / name, sep='\t', keep_default_na=False)
    assert set(v.view) == {'R_FOXP3_capture', 'E_effector_capture', 'Myeloid'}
    mcp_scope = v.groupby(['dataset', 'view']).size().reset_index(name='selected_cells').to_dict('records')
    legacy_path = root / 'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json'
    legacy = json.loads(legacy_path.read_text())
    legacy_scope = []
    for item in legacy['files']:
        p = Path(item['path']); assert sha(p) == item['sha256']
        obj = ad.read_h5ad(p, backed='r')
        assert 'core_position' in obj.obs
        legacy_scope.append(dict(major=item['major'], cells=obj.n_obs, genes=obj.n_vars,
                                 sources=sorted(set(obj.obs.gse_id.astype(str))),
                                 fixed34_gene_axis_coverage=len(set(obj.var_names) & set(genes)),
                                 same_immune_atlas_native_row_namespace=False))
        obj.file.close()
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    assert sha(atlas / 'manifests/immune_authority_release_public.json') == cfg['release_SHA256']
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    md = pd.DataFrame({k:np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset', 'library_id', 'major']})
    cm_path = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    cm = pd.read_csv(cm_path, sep='\t', keep_default_na=False)
    material = md.dataset + '::' + md.library_id
    use = material.isin(cm.sampleID).to_numpy()
    assert set(cm.sampleID) <= set(material)
    t_rows = cache['global_original_T_row_index']
    np.testing.assert_array_equal(t_rows, np.flatnonzero(md.major.eq('T')))
    source_order = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    present_sources = sorted(set(md.loc[use, 'dataset']))
    source_mask = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    status = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    maps = json.loads((old / 'private_source_gene_measurement_contract.json').read_text())
    excluded = {int(v) for item in maps for v in item['physical_columns'] if v >= 0}
    eligible = source_mask[:, [source_order.index(s) for s in present_sources]].all(axis=1)
    eligible &= ~np.char.startswith(canonical, 'AMBIGUOUS::') & (status != 'ambiguous_source_feature')
    eligible &= ~pd.Index(canonical).duplicated(keep=False)
    eligible[list(excluded)] = False
    coverage = []
    all_classes = (set(md.major) - {'T'}) | {
        'T_neither_capture', 'T_effector_only', 'T_reg_without_anchor_only',
        'T_reg_without_anchor_effector', 'T_anchored_reg_only',
        'T_anchored_reg_effector', 'T_core_not_captured'}
    for threshold in [1, 2]:
        labels = md.major.to_numpy().astype('U80')
        labels[t_rows] = label_T(cache['raw_gene_counts'], genes, panels, cfg['anchor_gene'], cfg['effector_genes'], threshold)
        selected = pd.DataFrame(dict(source=md.dataset[use].to_numpy(), library=md.library_id[use].to_numpy(), candidate=labels[use]))
        classes = sorted(all_classes)
        counts = selected.groupby(['source', 'candidate']).agg(cells=('candidate', 'size'), material_units=('library', 'nunique'))
        for src in present_sources:
            for candidate in classes:
                n, units = counts.loc[(src, candidate)].tolist() if (src, candidate) in counts.index else (0, 0)
                coverage.append(dict(source=src, candidate=candidate, threshold=threshold, cells=int(n), material_units=int(units)))
        assert sum(x['cells'] for x in coverage if x['threshold'] == threshold) == int(use.sum())
    out.mkdir(mode=0o700)
    result = dict(status='COMPLETE_EXISTING_CACHE_AND_MARKER_DISJOINT_REFERENCE_COVERAGE_ONLY', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), fixed_panel_contract_SHA256=sha(fixed_path), immune_release_SHA256=cfg['release_SHA256'], original_CM_material_SHA256=sha(cm_path), original_T_cache_SHA256=sha(old / 'private_fixed34_gene_T_counts_and_depth.npz'),
        existing_three_view_input_terminal_SHA256=sha(mcp / 'private_terminal.json'), existing_three_view_audit_SHA256=sha(audit_path), existing_three_view_source_scope=mcp_scope,
        existing_three_views_cover_NK_or_both_neither_T=False, legacy_all18_input_inventory_SHA256=sha(legacy_path), legacy_all18_scope=legacy_scope,
        new_exact_atlas_source_scope=present_sources, new_exact_atlas_material_cells=int(use.sum()), all_candidate_classes=sorted(all_classes),
        source_class_coverage=coverage, source_common_marker_disjoint_physical_genes=int(eligible.sum()), excluded_grouping_physical_columns=len(excluded),
        expression_fits_or_clinical_outcomes_read=False, all_real_arrays_RAM_only=True, per_cell_rows_or_gene_lists_exported=False,
        interpretation='Reference qualification only. RNA capture classes are operational, not formal lineage/functional or receptor truth; legacy cells are not assumed interchangeable with frozen immune-atlas rows.')
    (out / 'private_reference_coverage_qualification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__': main()
