"""Freeze all existing T/myeloid nodes and resources before external measurement.

No node pooling, gene expression, score, mapping or clinical effect is computed.
"""
import hashlib
import json
import os
import re
import shutil
import subprocess
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
    preflight = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2'
    proof = json.loads((preflight / 'private_independent_audit.json').read_text())
    assert proof['status'].startswith('PASS_ALL_ATLAS_LIBRARY_METADATA_')
    atlas_run = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    atlas = atlas_run / 'outputs/immune_validation_authority_v1.zarr'
    z = zarr.open_group(str(atlas), mode='r')
    columns = ['dataset', 'library_id', 'node_id', 'major']
    meta = pd.DataFrame({c: np.asarray(read_elem(z['obs'][c])).astype(str) for c in columns})
    text = read_elem(z['uns']['subcluster_catalog_json'])
    if isinstance(text, np.ndarray):
        text = text.item()
    catalog = pd.DataFrame(json.loads(text))
    selected_nodes = catalog[catalog.major.isin(['T', 'Myeloid'])]
    assert len(selected_nodes) == 29
    use = meta.node_id.isin(selected_nodes.node_id)
    keys = ['dataset', 'library_id', 'node_id', 'major']
    units = meta.loc[use].groupby(keys, sort=True).size().reset_index(name='n_cells')
    units['unit_index'] = np.arange(len(units))
    units = units.merge(selected_nodes[['node_id', 'subcluster', 'subtype_candidate', 'subtype_evidence_status', 'identity_promoted']], on='node_id', how='left', validate='many_to_one')
    assert units.node_id.nunique() == 29 and units.unit_index.tolist() == list(range(len(units)))
    assert units.n_cells.sum() == int(use.sum())
    index = pd.MultiIndex.from_frame(units[keys])
    rows = np.full(len(meta), -1, np.int32)
    rows[use.to_numpy()] = index.get_indexer(pd.MultiIndex.from_frame(meta.loc[use, keys]))
    assert (rows[use] >= 0).all() and (rows[~use] == -1).all()
    assert np.array_equal(np.bincount(rows[rows >= 0], minlength=len(units)), units.n_cells)
    out = root / 'existing_immune_fine_native_measurement_scope_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    units.to_csv(out / 'private_all_frozen_T_myeloid_node_units.tsv', sep='\t', index=False)
    np.save(out / 'private_original_cell_to_fine_node_unit.npy', rows)
    quota = subprocess.run(['quota', '-g', '-v', '-w', 'USER002'], capture_output=True, text=True, timeout=20)
    assert quota.returncode == 0
    match = re.search(r'^\s*/\S+\s+(\d+)\*?\s+\d+\s+(\d+)', quota.stdout, re.M)
    assert match, 'Exact KiB group quota not resolved'
    remaining = (int(match[2]) - int(match[1])) * 1024
    mem = {line.split(':')[0]: line.split(':')[1].strip() for line in Path('/proc/meminfo').read_text().splitlines()}
    available_ram = int(mem['MemAvailable'].split()[0]) * 1024
    native_features = int(z['layers']['counts'].attrs['shape'][1])
    dense_workspace_bytes = len(units) * native_features * 12
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), external_input_audit_sha256=sha(preflight / 'private_independent_audit.json'), atlas_release_sha256=sha(atlas_run / 'manifests/immune_authority_release_public.json'),
        scope='ALL29 frozen T/myeloid nodes, all source libraries, all node units including <20 cells. Node numbers and original identities unchanged; no pool substituted for scarce fine units.',
        prospective_native_processing='Read original counts and source masks, full native per-cell RNA denominator. Preserve raw gene availability versus true zero. Full native mean-log gene workspace may be held in RAM; avoid persisting another all-gene dense matrix.',
        prospective_outputs='Complete fine-unit metadata/quality, raw frozen programmes, all pre-fixed regulator target scores, limited native gene QA panel and coverage/identity flags. All four original molecular families retained. No effects or endpoint mapping in this stage.',
        interpretation='Candidate CD8/Treg groups are not newly promoted identities. Exact NSCLC endpoint correspondence requires a separate, effect-blind review; broader lineage extension cannot be called exact replication.',
        cells=len(meta), selected_original_cells=int(use.sum()), nodes=len(selected_nodes), native_units=len(units), native_features=native_features,
        source_CSR_stored_entries=int(z['layers']['counts']['data'].shape[0]), dense_workspace_bytes=dense_workspace_bytes, available_ram_bytes=available_ram, exact_group_quota_remaining_bytes=remaining, filesystem_free_bytes=shutil.disk_usage(out).free,
        quota_reserve_decimal_bytes=500000000000, future_output_cap_bytes=256*1024**2, quota_extra_headroom_bytes=1000000000,
        resource_gate=available_ram > dense_workspace_bytes*2.5 + 4*1024**3 and remaining > 500000000000 + 256*1024**2 + 1000000000,
        no_expression_scores_models_or_download=True)
    (out / 'plan_and_resource_preflight.json').write_text(json.dumps(plan, indent=2) + '\n')
    result = dict(status='FROZEN_ALL_ORIGINAL_FINE_NODE_UNITS_AND_RESOURCE_GATE_RECORDED_NO_EXPRESSION_MEASURED', finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), nodes=len(selected_nodes), native_units=len(units), selected_original_cells=int(use.sum()), resource_gate=plan['resource_gate'], files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_labels_or_scores_or_models_changed=True)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
