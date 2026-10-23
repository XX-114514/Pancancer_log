"""Read the frozen private pilot in RAM with full native RNA denominators."""
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from annotation_representation_pilot_v1 import ROUTES, feature_columns, normalized_log_counts
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import gather_csr_rows


def load_pilot(pilot_root, manifest, checkpoint=None, block_size=2048, feature_universe_binding=None):
    pilot = Path(pilot_root)
    spec = json.loads((pilot / 'config/pilot_spec_v1.json').read_text())
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    plan = json.loads((pilot / 'config/model_features_and_heldouts_v1.json').read_text())
    reader = FullContextAtlas(manifest)
    genes = plan['common_model_feature_candidates']
    heldout = sorted({g for panels in plan['heldout_RNA_readouts'].values() for panel in panels.values() for g in panel})
    if feature_universe_binding is not None:
        from common_measured_feature_contract_v1 import restrict_candidates
        genes = restrict_candidates(genes, heldout, feature_universe_binding)
    common_readouts = {
        name: [g for g in panel if all(v[g] for v in plan['source_measurement_coverage'].values())]
        for name, panel in plan['heldout_RNA_readouts']['state'].items()}
    folds = {row['dataset']: row['fold_group'] for row in spec['sources']}
    matrices, readouts, totals, sources, routes, fine = [], [], [], [], [], []
    for entry in inventory:
        ds = entry['dataset']
        path = pilot / entry['selector_file']
        if sha256_file(path) != entry['selector_sha256']:raise ValueError('Pilot selector changed')
        positions = np.load(path, allow_pickle=False)
        lookup = reader.feature_lookup(ds)
        columns = feature_columns(lookup, genes, heldout)
        measured = np.array([len(lookup.get(g, set())) == 1 for g in heldout])
        available = np.flatnonzero(measured)
        readout_columns = np.array([next(iter(lookup[heldout[j]])) for j in available], dtype=np.int64)
        z = reader.store(ds)
        key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        xs, ys, ts = [], [], []
        for low in range(0, len(positions), block_size):
            local = positions[low:low + block_size]
            raw = gather_csr_rows(z[key], reader.rows(ds)[local])
            if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                    not np.equal(raw.data, np.floor(raw.data)).all()):raise ValueError('Native count domain changed')
            total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
            x, ok = normalized_log_counts(raw[:, columns], total)
            y = np.full((len(local), len(heldout)), np.nan, dtype=np.float32)
            y[:, available] = raw[:, readout_columns].toarray()
            scale = np.zeros(len(total), dtype=np.float32);scale[ok] = 10000 / total[ok]
            y = np.log1p(y * scale[:, None])
            xs.append(x);ys.append(y);ts.append(total)
        labels = reader.frozen_labels(ds)
        matrices.append(sp.vstack(xs, format='csr'));readouts.append(np.vstack(ys));totals.append(np.concatenate(ts))
        sources.extend([ds] * len(positions))
        routes.extend([ROUTES.get(str(v), 'unresolved_context') for v in labels['major'][positions]])
        fine.extend(labels['fine'][positions].astype(str))
        if checkpoint:checkpoint('native_holdout_input_loaded', dataset=ds, n_cells=len(positions))
    x = sp.vstack(matrices, format='csr');readout = np.vstack(readouts);source = np.asarray(sources)
    state = {}
    for name, panel in common_readouts.items():
        if len(panel) >= 3:
            state[name] = readout[:, [heldout.index(g) for g in panel]].mean(axis=1)
        else:state[name] = np.full(len(source), np.nan)
    gamma_genes = ['TRDC', 'TRGC1', 'TRGC2', 'CD3D']
    a = {g: readout[:, heldout.index(g)] for g in gamma_genes}
    gamma_measured = np.all(np.isfinite(np.column_stack(list(a.values()))), axis=1)
    gamma = np.full(len(source), np.nan)
    gamma[gamma_measured] = ((a['TRDC'] > 0) & ((a['TRGC1'] > 0) | (a['TRGC2'] > 0)) & (a['CD3D'] > 0))[gamma_measured]
    return {'x': x, 'genes': genes, 'native_total': np.concatenate(totals), 'source': source,
            'route': np.asarray(routes), 'fine': np.asarray(fine),
            'fold': np.asarray([folds[s] for s in source]), 'state': state,
            'state_genes': common_readouts, 'gamma_RNA_candidate': gamma,
            'readout_genes': heldout, 'native_readout_logCP10K': readout}
