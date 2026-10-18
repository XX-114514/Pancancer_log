"""Finite conditional spatial mixing of fixed functional RNA capture patterns.

Not exact CM reconstruction, direct signalling, patient inference or new labels.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from qualify_GSE299886_fixed_panels_negative_controls_v1 import verify
from repair_GSE299886_original_coordinate_sidecars_v1 import dense_rows


BASE_SEED = 2026101801
N_PERM = 2000
SCALES = (2, 5, 10)
MINIMUM = dict(T_cells=100, arm_cells=30, edges=100, movable_fraction=0.5)


def adjustment(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def cross_count(state, edges):
    a, b = state[edges[:, 0]], state[edges[:, 1]]
    return int((((a & 1) > 0) & ((b & 2) > 0)).sum() + (((b & 1) > 0) & ((a & 2) > 0)).sum())


def strata_for(xy, total, negative, unit):
    grid = np.floor((xy - xy.min(axis=0)) / (50 * unit)).astype(np.int64)
    _, block = np.unique(grid, axis=0, return_inverse=True)
    depth = np.zeros(len(xy), np.int16)
    for b in np.unique(block):
        ix = np.flatnonzero(block == b)
        cuts = np.unique(np.quantile(total[ix], [.25, .5, .75]))
        depth[ix] = np.searchsorted(cuts, total[ix], side='left')
    _, group = np.unique(np.column_stack([block, depth, negative > 0]), axis=0, return_inverse=True)
    return group, grid, depth


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 64 * 1024**2)
    coordinate = root / 'GSE299886_original_coordinate_sidecars_20261018_v1'
    qc = root / 'GSE299886_fixed_panels_negative_controls_20261018_v1'
    verify(coordinate)
    verify(qc)
    reports = json.loads((coordinate / 'private_all_source_coordinate_reconstruction_receipt.json').read_text())
    qc_records = json.loads((qc / 'private_all3_negative_probe_geometry_QC.json').read_text())
    schemas = json.loads((root / 'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json').read_text())
    state_root = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    ap = state_root / 'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    axes = json.loads(ap.read_text())
    fp = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    frozen = json.loads(fp.read_text())
    out = root / 'CM02_MERFISH_conditional_spatial_mixing_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        budget=budget, coordinate_terminal_sha256=sha(coordinate / 'private_terminal.json'),
        QC_terminal_sha256=sha(qc / 'private_terminal.json'), fixed_panels_sha256=sha(fp),
        scope='All3 original fixed-panel-eligible assays; all18 slots = 3 assays x 3 radii x 2 RNA thresholds. No ligand/receptor selection or original CM/clinical refit.',
        lineage='At least2 distinct original pan_T genes with >=1 count; no author cell-state label or inferred patient identity.',
        arms='At least2 distinct genes from each complete unchanged regulatory_associated/cytotoxic panel. Main gene threshold>=1; stricter sensitivity>=2. Double-positive states retained.',
        spatial='T-to-T undirected pairs within [2,5,10] times median nearest-other-ALL-cell distance. Primary5. Source coordinates scale-invariant; no micrometre claim, self-edge or graph crossing a larger tissue gap.',
        null='Joint R/C four-state labels permuted within source-local square50-unit blocks x within-block RNA-total quartiles x negative-probe-present/absent. T positions, graph and joint within-stratum arm counts/overlap fixed. Zero negative probes cannot establish negative background measurement.',
        statistic='Directed R->C edge opportunities, cross_count=sum(R_i*C_j+R_j*C_i) for undirected edges; observed/permutation mean. Double-positive endpoints may contribute twice. Local composition is held approximately by blocks, not removed as proven biological confounding.',
        inference='Two-sided doubled smaller inclusive Monte Carlo tail with +1 over2001,2000 permutations. Full18-slot BH/BY, unavailable internal P=1 but NA shown. Empirical permutation reference quantiles are null intervals, NOT effect CI or patient uncertainty.',
        seed=BASE_SEED, permutations=N_PERM, minimum=MINIMUM,
        limits='Conditional within-assay exploratory organization only. FOV/core/patients and physical units incompletely reviewed, paired arrays may overlap people, not3 independent cohorts. Author labels, TMA selection, RNA sparsity, capture background and residual tissue structure can affect mixing. No direct binding, validated Treg state, ICI mechanism, clinical prediction or exact learned CM claim.')
    (out / 'plan_before_any_spatial_effects.json').write_text(json.dumps(plan, indent=2) + '\n')
    results = []
    for number, report in enumerate(reports):
        sid = report['sample_id']
        schema = next(x for x in schemas if x['sample_id'] == sid)
        axis = next(x for x in axes if x['sample_id'] == sid)
        path = Path(schema['analysis_file'])
        assert sha(path) == report['matrix_SHA256'] == axis['input_sha256']
        mapping = axis['symbol_to_feature_index']
        required = ['pan_T', 'regulatory_associated', 'cytotoxic']
        assert all(g in mapping for k in required for g in frozen['panels'][k])
        z = np.load(coordinate / report['output'], allow_pickle=False)
        quality = np.load(qc / ('private_source_%02d_negative_QC.npz' % number), allow_pickle=False)
        n = len(z['xy_source'])
        detection = {(k, threshold): np.zeros(n, np.int16) for k in required for threshold in [1, 2]}
        with h5py.File(path, 'r') as h:
            for start in range(0, n, 2048):
                stop = min(n, start + 2048)
                v = dense_rows(h['X'], start, stop, 500)
                for k, threshold in detection:
                    detection[k, threshold][start:stop] = (v[:, [mapping[g] for g in frozen['panels'][k]]] >= threshold).sum(axis=1)
        select = detection['pan_T', 1] >= 2
        index = np.flatnonzero(select)
        xy = z['xy_source'][select]
        unit = float(np.median(cKDTree(z['xy_source']).query(z['xy_source'], k=2, workers=1)[0][:, 1]))
        assert unit > 0 and np.isfinite(unit)
        group, grid, depth = strata_for(xy, quality['biological_total'][select], quality['negative_total'][select], unit)
        groups = [np.flatnonzero(group == g) for g in np.unique(group)]
        graphs = {scale: cKDTree(xy).query_pairs(unit * scale, output_type='ndarray').astype(np.int32) for scale in SCALES}
        states = {threshold: (detection['regulatory_associated', threshold][select] >= 2).astype(np.int8) + 2 * (detection['cytotoxic', threshold][select] >= 2).astype(np.int8) for threshold in [1, 2]}
        arrays = dict(existing_T_row_index=index, xy=xy, permutation_stratum=group, spatial_grid=grid, RNA_depth_quartile=depth, unit=np.array(unit))
        arrays.update({'edges_scale%d' % k: v for k, v in graphs.items()})
        arrays.update({'state_threshold%d' % k: v for k, v in states.items()})
        for threshold, state in states.items():
            move = [ix for ix in groups if len(np.unique(state[ix])) > 1]
            movable = sum(len(ix) for ix in move) / len(state) if len(state) else 0
            nr, nc = int(((state & 1) > 0).sum()), int(((state & 2) > 0).sum())
            seed = BASE_SEED + number * 101 + threshold * 31
            reasons = []
            if len(state) < MINIMUM['T_cells']: reasons.append('insufficient_T_capture')
            if min(nr, nc) < MINIMUM['arm_cells']: reasons.append('insufficient_arm_capture')
            if movable < MINIMUM['movable_fraction']: reasons.append('insufficient_joint_state_exchangeability')
            rng = np.random.default_rng(seed)
            null = np.full((N_PERM, len(SCALES)), -1, np.int64)
            active = [j for j, scale in enumerate(SCALES) if not reasons and len(graphs[scale]) >= MINIMUM['edges']]
            if active:
                for permutation in range(N_PERM):
                    shuffled = state.copy()
                    for ix in move:
                        shuffled[ix] = rng.permutation(state[ix])
                    for j in active:
                        null[permutation, j] = cross_count(shuffled, graphs[SCALES[j]])
            arrays['null_threshold%d' % threshold] = null
            for j, scale in enumerate(SCALES):
                local_reasons = reasons.copy()
                if len(graphs[scale]) < MINIMUM['edges']: local_reasons.append('insufficient_T_pair_edges')
                counts = null[:, j]
                valid = not local_reasons and (counts >= 0).all() and counts.std() > 0 and counts.mean() > 0
                if not local_reasons and not valid: local_reasons.append('degenerate_permutation_reference')
                observed = cross_count(state, graphs[scale])
                p = min(1., 2 * min((1 + (counts >= observed).sum()) / (N_PERM + 1), (1 + (counts <= observed).sum()) / (N_PERM + 1))) if valid else None
                results.append(dict(dataset_id='GSE299886', sample_id=sid, technical_input_number=number,
                    gene_count_threshold=threshold, radius_scale=scale, primary=(threshold == 1 and scale == 5),
                    T_cells=len(state), regulatory_capture_cells=nr, cytotoxic_capture_cells=nc,
                    double_positive_capture_cells=int((state == 3).sum()), graph_undirected_edges=len(graphs[scale]),
                    movable_fraction=movable, permutation_strata=len(groups), cross_count_observed=observed,
                    cross_count_null_mean=float(counts.mean()) if valid else None,
                    cross_count_null_q025=float(np.quantile(counts, .025)) if valid else None,
                    cross_count_null_q975=float(np.quantile(counts, .975)) if valid else None,
                    mixing_observed_over_null=float(observed / counts.mean()) if valid else None,
                    p_two_sided=p, available=bool(valid), reason=';'.join(local_reasons), seed=seed,
                    negative_background_uninformative=qc_records[number]['negative_transcripts'] == 0,
                    no_independent_patient_or_direct_communication_claim=True))
        np.savez_compressed(out / ('private_source_%02d_graph_state_and_null.npz' % number), **arrays)
        assert sha(path) == report['matrix_SHA256']
    assert len(results) == 18
    d = pd.DataFrame(results)
    pp = d.p_two_sided.fillna(1).to_numpy(float)
    bh = adjustment(pp)
    by = np.minimum(bh * np.sum(1 / np.arange(1, 19)), 1)
    d['BH_all18'] = np.where(d.available, bh, np.nan)
    d['BY_all18'] = np.where(d.available, by, np.nan)
    d.to_csv(out / 'private_all18_conditional_spatial_mixing.tsv', sep='\t', index=False)
    files = [p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files) < 64 * 1024**2
    terminal = dict(status='COMPLETE_ALL18_EXPLORATORY_WITHIN_ASSAY_SPATIAL_SLOTS_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), available=int(d.available.sum()),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}, independent_statistical_audit=False,
        biological_mechanism_or_patient_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
