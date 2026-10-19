"""Independent tuple/Counter membership and accepted-T-cache capture replay."""
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    source = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_AUDIT_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    terminal = json.loads((source / 'private_terminal.json').read_text())
    for name, proof in terminal['files'].items():
        assert sha(source / name) == proof['SHA256']
    d = pd.DataFrame(json.loads((source / 'private_all_source_scope_major_panel_descriptions.json').read_text()))
    j = pd.DataFrame(json.loads((source / 'private_all_source_scope_major_two_arm_overlaps.json').read_text()))
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    p = source / 'private_frozen_plan.json'
    plan = json.loads(p.read_text())
    cm_path = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    assert sha(cm_path) == plan['CM_material_metadata_SHA256']
    cm = pd.read_csv(cm_path, sep='\t', keep_default_na=False)
    cm_keys = {tuple(v.split('::', 1)) for v in cm.sampleID}
    pre_keys = {tuple(v.split('::', 1)) for v in cm.loc[cm.timepoint.eq('pre'), 'sampleID']}
    assert len(cm_keys) == len(cm) and pre_keys <= cm_keys
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    meta = {k: np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset', 'library_id', 'major']}
    counters = [Counter(), Counter(), Counter()]
    for src, lib, major in zip(meta['dataset'], meta['library_id'], meta['major']):
        key = (src, lib)
        counters[0][src, major] += 1
        if key in cm_keys: counters[1][src, major] += 1
        if key in pre_keys: counters[2][src, major] += 1
    scopes = plan['scopes']
    full_count_positions = 0
    for x in d.to_dict('records'):
        assert x['original_cells'] == counters[scopes.index(x['scope'])][x['source'], x['major']]
        full_count_positions += 1
    old = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    old_manifest = json.loads((old / 'private_terminal.json').read_text())
    for name in ['private_original_T_unit_scope.tsv', 'private_source_gene_measurement_contract.json']:
        assert sha(old / name) == old_manifest['files'][name]['sha256']
    assert sha(old / 'private_fixed34_gene_T_counts_and_depth.npz') == plan['original_T_cache_SHA256']
    a = dict(np.load(old / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False))
    u = pd.read_csv(old / 'private_original_T_unit_scope.tsv', sep='\t', keep_default_na=False)
    unit = a['original_T_unit_index']; raw = a['raw_gene_counts']; genes = a['genes'].tolist()
    fixed = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed) == plan['fixed_panel_contract_SHA256']
    panels = json.loads(fixed.read_text())['panels']
    t_sources = u.dataset.to_numpy()[unit]
    t_keys = [(s, lib) for s, lib in zip(u.dataset.to_numpy()[unit], u.library_id.to_numpy()[unit])]
    masks = [np.ones(len(unit), bool), np.fromiter((x in cm_keys for x in t_keys), bool, len(unit)), np.fromiter((x in pre_keys for x in t_keys), bool, len(unit))]
    maps = {x['dataset']: x for x in json.loads((old / 'private_source_gene_measurement_contract.json').read_text())}
    lookup = d.set_index(['source', 'scope', 'major', 'panel'])
    joint_lookup = j.set_index(['source', 'scope', 'major', 'panelA', 'panelB', 'threshold'])
    checked = 0
    for si, scope_name in enumerate(scopes):
        for src in sorted(d.source.unique()):
            mask = masks[si] & (t_sources == src)
            for threshold in [1, 2]:
                positive = {}
                for panel, gg in panels.items():
                    record = lookup.loc[(src, scope_name, 'T', panel)]
                    if not maps[src]['complete_panels'][panel]:
                        assert pd.isna(record['captured_cells_threshold' + str(threshold)])
                        continue
                    # Explicit gene loop, unlike producer's column-wise matrix count.
                    detections = np.zeros(int(mask.sum()), np.int8)
                    for gene in gg: detections += raw[mask, genes.index(gene)] >= threshold
                    positive[panel] = detections >= 2
                    assert int(positive[panel].sum()) == record['captured_cells_threshold' + str(threshold)]
                    checked += 1
                for left, right in plan['fixed_panel_pairs']:
                    record = joint_lookup.loc[(src, scope_name, 'T', left, right, threshold)]
                    if left in positive and right in positive:
                        assert int(np.logical_and(positive[left], positive[right]).sum()) == record.both_captured_cells
                        checked += 1
                    else: assert pd.isna(record.both_captured_cells)
    out.mkdir(mode=0o700)
    result = dict(status='ACCEPTED_ALL_SCOPE_METADATA_TUPLE_COUNTER_AND_ALL_ORIGINAL_T_PANEL_CAPTURE_REPLAY', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(source / 'private_terminal.json'), all_metadata_scope_panel_positions_replayed=full_count_positions, original_T_capture_and_pair_positions_replayed=checked, original_project_read_only=True, real_metadata_and_RNA_arrays_RAM_only=True, independent_full_nonT_raw_replay=False, individual_rows_or_gene_arrays_exported=False)
    (out / 'private_independent_acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
