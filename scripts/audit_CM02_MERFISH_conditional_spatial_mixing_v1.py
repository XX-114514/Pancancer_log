"""Independent raw sparse-RNA, geometry and full-random-stream acceptance.

No producer statistic, graph, stratification or FDR function is imported.
"""
import hashlib
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def count(state, edges):
    r = (state % 2).astype(np.int64)
    c = (state // 2).astype(np.int64)
    return int(np.dot(r[edges[:, 0]], c[edges[:, 1]]) + np.dot(r[edges[:, 1]], c[edges[:, 0]]))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 8 * 1024**2)
    source = root / 'CM02_MERFISH_conditional_spatial_mixing_20261018_v1'
    terminal = json.loads((source / 'private_terminal.json').read_text())
    assert sha(source / 'private_terminal.json') == '3fca0f4393ed3e153579802d3f17b42739cf0120f73af243a64f004261c92659'
    for n, v in terminal['files'].items():
        assert sha(source / n) == v['sha256']
    plan = json.loads((source / 'plan_before_any_spatial_effects.json').read_text())
    assert sha(root / 'CRC_fine_pair_DIALOGUE_code_20261014_v1/explore_CM02_MERFISH_conditional_spatial_mixing_v1.py') == plan['source_sha256']
    out = root / 'CM02_MERFISH_conditional_spatial_mixing_audit_20261018_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    (out / 'plan_before_complete_statistical_acceptance.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        producer_terminal_sha256=sha(source / 'private_terminal.json'),
        scope='Re-read every original input sparse RNA; independently derive captured T/arm states, graph all edges via balls rather than query_pairs, strata via pandas, all2000 random reference counts, two-sided tails and statsmodels full18 BH/BY. No data refit or new candidate.',
        synthetic='Every6 unique four-cell label allocation, positive bipartite and negative within-state graphs, all-degenerate null and arm-count insufficient-case gates.'), indent=2) + '\n')
    coordinate = root / 'GSE299886_original_coordinate_sidecars_20261018_v1'
    qc = root / 'GSE299886_fixed_panels_negative_controls_20261018_v1'
    assert sha(coordinate / 'private_terminal.json') == plan['coordinate_terminal_sha256']
    assert sha(qc / 'private_terminal.json') == plan['QC_terminal_sha256']
    for folder in [coordinate, qc]:
        for n, v in json.loads((folder / 'private_terminal.json').read_text())['files'].items():
            assert sha(folder / n) == v['sha256']
    fp = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fp) == plan['fixed_panels_sha256']
    panels = json.loads(fp.read_text())['panels']
    reports = json.loads((coordinate / 'private_all_source_coordinate_reconstruction_receipt.json').read_text())
    schemas = json.loads((root / 'CM02_existing_targeted_spatial_schema_20261018_v2/private_all_eligible_targeted_input_schema.json').read_text())
    state_root = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    ap = state_root / 'state_spatial_proxy_calibration_inputs_20261016_v1/fixed_all14_panel_gene_axes.json'
    qc_plan = json.loads((qc / 'plan_before_all_source_capture_and_negative_probe_QC.json').read_text())
    assert sha(ap) == qc_plan['axes_sha256']
    axes = json.loads(ap.read_text())
    d = pd.read_csv(source / 'private_all18_conditional_spatial_mixing.tsv', sep='\t', keep_default_na=False)
    pvalues = []
    checks = []
    for number, report in enumerate(reports):
        row_schema = next(x for x in schemas if x['sample_id'] == report['sample_id'])
        path = Path(row_schema['analysis_file'])
        assert sha(path) == report['matrix_SHA256']
        mapping = next(x for x in axes if x['sample_id'] == report['sample_id'])['symbol_to_feature_index']
        data = np.load(source / ('private_source_%02d_graph_state_and_null.npz' % number), allow_pickle=False)
        original = np.load(coordinate / report['output'], allow_pickle=False)
        quality = np.load(qc / ('private_source_%02d_negative_QC.npz' % number), allow_pickle=False)
        with h5py.File(path, 'r') as h:
            from anndata._io.specs import read_elem
            var_index = h['var'].attrs['_index']
            declared = [v.decode() if isinstance(v, bytes) else str(v) for v in read_elem(h['var'][var_index])]
            requested = [g for k in ['pan_T', 'regulatory_associated', 'cytotoxic'] for g in panels[k]]
            if not set(requested) <= set(declared):
                declared = [v.decode() if isinstance(v, bytes) else str(v) for v in read_elem(h['var']['gene_symbol'])]
            assert all(declared[mapping[g]] == g for g in requested)
            g = h['X']
            matrix = csr_matrix((g['data'][...], g['indices'][...], g['indptr'][...]), shape=tuple(g.attrs['shape']))
        raw_panels = {k: matrix[:, [mapping[gene] for gene in panels[k]]].toarray() for k in ['pan_T', 'regulatory_associated', 'cytotoxic']}
        index = np.flatnonzero((raw_panels['pan_T'] > 0).sum(axis=1) >= 2)
        assert np.array_equal(index, data['existing_T_row_index'])
        xy = original['xy_source'][index]
        assert np.array_equal(xy, data['xy'])
        all_tree = cKDTree(original['xy_source'])
        unit = float(np.median(all_tree.query(original['xy_source'], k=2, workers=1)[0][:, 1]))
        assert unit == float(data['unit'])
        box = np.floor((xy - xy.min(axis=0)) / (50 * unit)).astype(np.int64)
        assert np.array_equal(box, data['spatial_grid'])
        frame = pd.DataFrame(dict(x=box[:, 0], y=box[:, 1], total=quality['biological_total'][index],
                                  negative=quality['negative_total'][index] > 0))
        frame['depth'] = 0
        for _, loc in frame.groupby(['x', 'y'], sort=True).groups.items():
            cuts = np.unique(frame.loc[loc, 'total'].quantile([.25, .5, .75]).to_numpy())
            frame.loc[loc, 'depth'] = [int(np.count_nonzero(cuts < value)) for value in frame.loc[loc, 'total']]
        assert np.array_equal(frame.depth, data['RNA_depth_quartile'])
        # Sorted x/y/depth/negative tuples are equivalent to the producer's sorted block indices.
        tuples = list(zip(frame.x, frame.y, frame.depth, frame.negative))
        labels = {key: j for j, key in enumerate(sorted(set(tuples)))}
        strata = np.array([labels[key] for key in tuples])
        assert np.array_equal(strata, data['permutation_stratum'])
        groups = [np.flatnonzero(strata == j) for j in sorted(set(strata))]
        graph = {}
        tree = cKDTree(xy)
        for scale in [2, 5, 10]:
            pairs = [(i, j) for i, neighbors in enumerate(tree.query_ball_point(xy, unit * scale)) for j in neighbors if j > i]
            graph[scale] = np.asarray(sorted(pairs), np.int64).reshape(-1, 2)
            recorded = data['edges_scale%d' % scale]
            assert np.array_equal(graph[scale], recorded[np.lexsort((recorded[:, 1], recorded[:, 0]))])
        for threshold in [1, 2]:
            r = (raw_panels['regulatory_associated'][index] >= threshold).sum(axis=1) >= 2
            c = (raw_panels['cytotoxic'][index] >= threshold).sum(axis=1) >= 2
            state = r.astype(np.int8) + 2 * c.astype(np.int8)
            assert np.array_equal(state, data['state_threshold%d' % threshold])
            moving = [g for g in groups if len(set(state[g])) > 1]
            movable = sum(map(len, moving)) / len(state)
            gate = len(state) >= 100 and min(int(r.sum()), int(c.sum())) >= 30 and movable >= .5
            expected = data['null_threshold%d' % threshold]
            active = [j for j, scale in enumerate([2, 5, 10]) if gate and len(graph[scale]) >= 100]
            rng = np.random.default_rng(2026101801 + number * 101 + threshold * 31)
            for permutation in range(2000):
                if not active: break
                permuted = state.copy()
                for g in moving:
                    permuted[g] = state[g][rng.permutation(len(g))]
                for j in active:
                    assert count(permuted, graph[[2, 5, 10][j]]) == expected[permutation, j]
            for j, scale in enumerate([2, 5, 10]):
                row = d[(d.technical_input_number == number) & (d.gene_count_threshold == threshold) & (d.radius_scale == scale)]
                assert len(row) == 1
                row = row.iloc[0]
                observed = count(state, graph[scale])
                assert observed == int(row.cross_count_observed)
                counts = expected[:, j]
                valid = j in active and np.std(counts) > 0 and np.mean(counts) > 0
                assert valid == (str(row.available).lower() == 'true')
                if valid:
                    p = min(1., 2 * min((int(np.count_nonzero(counts <= observed)) + 1) / 2001,
                                        (int(np.count_nonzero(counts >= observed)) + 1) / 2001))
                    assert abs(p - float(row.p_two_sided)) < 1e-13
                    assert abs(observed / np.mean(counts) - float(row.mixing_observed_over_null)) < 1e-12
                    assert np.isclose(np.mean(counts), float(row.cross_count_null_mean))
                    assert np.isclose(np.quantile(counts, .025), float(row.cross_count_null_q025))
                    assert np.isclose(np.quantile(counts, .975), float(row.cross_count_null_q975))
                else:
                    p = 1.
                    assert row.p_two_sided == '' and row.BH_all18 == '' and row.BY_all18 == ''
                pvalues.append((int(row.name), p))
                checks.append(dict(technical_input_number=number, threshold=threshold, radius=scale,
                    raw_RNA_state_recomputed=True, all_graph_edges_rebuilt=True,
                    all_random_reference_counts_replayed=bool(j in active), available=bool(valid)))
        assert sha(path) == report['matrix_SHA256']
        del matrix, raw_panels
    assert len(pvalues) == len(d) == 18
    pp = np.array([p for _, p in sorted(pvalues)])
    for method, col in [('fdr_bh', 'BH_all18'), ('fdr_by', 'BY_all18')]:
        values = multipletests(pp, method=method)[1]
        for i, row in d.iterrows():
            if str(row.available).lower() == 'true': assert abs(values[i] - float(row[col])) < 1e-12
    # Exact small graph orientation, overlap and self-exclusion checks.
    small = np.array([[0, 2], [0, 3], [1, 2], [1, 3]])
    allocations = sorted(set(itertools.permutations([1, 1, 2, 2])))
    exact = [count(np.array(a), small) for a in allocations]
    assert sorted(exact) == [2, 2, 2, 2, 4, 4]
    assert count(np.array([3, 3]), np.array([[0, 1]])) == 2
    labels = np.array([1] * 50 + [2] * 50)
    positive = np.array([(i, j) for i in range(50) for j in range(50, 100)])
    negative = np.array([(i, j) for i in range(100) for j in range(i + 1, 100) if (i < 50) == (j < 50)])
    rng = np.random.default_rng(20261018)
    pos_null, neg_null = [], []
    for _ in range(2000):
        shuffled = labels[rng.permutation(100)]
        pos_null.append(count(shuffled, positive)); neg_null.append(count(shuffled, negative))
    assert count(labels, positive) > max(pos_null) and count(labels, negative) < min(neg_null)
    assert not (100 >= 100 and min(29, 30) >= 30 and .9 >= .5)
    degenerate = np.array([3] * 100)
    assert count(degenerate, positive) == 2 * len(positive)
    (out / 'private_complete_independent_acceptance.json').write_text(json.dumps(dict(
        status='ALL18_RAW_RNA_GEOMETRY_STRATA_FULL_RANDOM_STREAM_AND_BH_BY_REPLAY_ACCEPTED',
        checks=checks, synthetic_positive_negative_and_unavailable=True, exact_small_graph_enumeration=exact,
        note='Independent implementations, same deposited data and fixed random generator; not independent scientific replication or a proof of exchangeability.'), indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    t = dict(status='ACCEPTED_ALL18_CONDITIONAL_SPATIAL_NUMERICS_NOT_MECHANISM',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        producer_terminal_sha256=sha(source / 'private_terminal.json'),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
        independent_patient_or_biological_mechanism_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(t, indent=2) + '\n')
    print(json.dumps(dict(status=t['status'], completed_UTC=t['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
