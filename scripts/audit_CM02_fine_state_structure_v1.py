"""Replay all fine-state descriptions with Counter counts and SVD residuals.

Shares the pinned metadata decoder and already audited cell join materializer;
does not use the producer's count matrices, effect, or resampling functions.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def effect(y, c):
    if len(y) < 20:
        return np.full((7, 6), np.nan)
    d = c.copy()
    for j in [1, 2]:
        d[:, j] = pd.Series(d[:, j]).rank(method='average').to_numpy()
    d[:, 1:3] = (d[:, 1:3] - d[:, 1:3].mean(axis=0)) / np.maximum(d[:, 1:3].std(axis=0), 1e-15)
    u, s, _ = np.linalg.svd(d, full_matrices=False)
    rank = int(np.sum(s > max(d.shape) * np.finfo(float).eps * s[0]))
    if rank != 5 or len(y) - rank - 1 < 10 or s[0] / s[-1] > 1e6:
        return np.full((7, 6), np.nan)
    a = pd.DataFrame(y).rank(axis=0, method='average').to_numpy()
    e = a - u @ (u.T @ a)
    v = np.sqrt(np.sum(e * e, axis=0))
    result = e[:, :7].T @ e[:, 7:] / np.maximum(v[:7, None] * v[None, 7:], 1e-15)
    result[(v[:7, None] < 1e-10) | (v[None, 7:] < 1e-10)] = np.nan
    return result


def equal(a, b):
    assert (a is None) == (not np.isfinite(b))
    if a is not None:
        assert abs(a - b) < 1e-9


def main():
    parent = Path(os.environ['PANCANCER_FINE_STRUCTURE_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    t = json.loads((parent / 'private_terminal.json').read_text())
    for name, v in t['files'].items():
        assert sha(parent / name) == v['SHA256']
    source = parent.parent / 'describe_CM02_fine_state_structure_v1.py'
    assert sha(source) == t['source_SHA256']
    spec = importlib.util.spec_from_file_location('accepted_reference_inputs', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    labels, tab, cur, contexts, unique, nodes, catalog = module.load_inputs()
    indexed = {(r['context'], r['phase'], r['CD4_author_state'], r['CD8_author_state']): r for r in json.loads((parent / 'private_all168_fine_state_associations.json').read_text())}
    assert len(indexed) == 168
    keys = set()
    metadata = tab.to_dict('records')
    total = Counter(cur.library_id)
    current_t = Counter(cur.loc[cur.major.eq('T'), 'library_id'])
    for ci, context in enumerate(['source_fine_reference', 'CM02_naming_anchors_only']):
        counted = Counter()
        for r in unique.itertuples():
            if context == 'source_fine_reference' or (r.reference == 'CD4_reference' and r.node_id == nodes[0]) or (r.reference == 'CD8_reference' and r.node_id == nodes[1]):
                counted[r.library_id, r.label] += 1
        for pi, phase in enumerate(['pre', 'post']):
            g = [r for r in metadata if r['timepoint_registry'] == phase]
            patient_counts = Counter(r['patient_registry'] for r in g)
            ready = []
            for r in g:
                n4 = sum(counted[r['library_id'], k] for k in labels.CD4)
                n8 = sum(counted[r['library_id'], k] for k in labels.CD8)
                if n4 >= 20 and n8 >= 20 and patient_counts[r['patient_registry']] == 1 and r['neo_type'] in ['D', 'D+T'] and r['HPV/p16'] in ['positive', 'negative']:
                    ready.append((r, n4, n8))
            ready.sort(key=lambda z: z[0]['patient_registry'])
            n = len(ready)
            x = np.array([[counted[r['library_id'], k] for k in labels.CD4 + labels.CD8] for r, _, _ in ready], float).reshape((-1, 13))
            n4, n8 = x[:, :7].sum(axis=1), x[:, 7:].sum(axis=1)
            y = x / np.column_stack([np.repeat(n4[:, None], 7, axis=1), np.repeat(n8[:, None], 6, axis=1)])
            tx = np.array([r['neo_type'] == 'D+T' for r, _, _ in ready], int)
            hpv = np.array([r['HPV/p16'] == 'positive' for r, _, _ in ready], int)
            c = np.column_stack([np.ones(n), np.log1p(n4), np.log1p(n8), tx, hpv])
            point = effect(y, c)
            denominators = {'all_unique_reference_T': n4 + n8, 'all_current_T': np.array([current_t[r['library_id']] for r, _, _ in ready], float), 'all_current_source_cells': np.array([total[r['library_id']] for r, _, _ in ready], float)}
            sensitivity = {k: effect(x / v[:, None], c) for k, v in denominators.items()}
            loo = np.asarray([effect(np.delete(y, k, axis=0), np.delete(c, k, axis=0)) for k in range(n)]).reshape((-1, 7, 6))
            rng = np.random.default_rng(2026102103 + 10 * ci + pi)
            strata = 2 * tx + hpv
            blocks = [np.where(strata == k)[0] for k in sorted(set(strata))]
            boot = []
            if np.isfinite(point).any():
                for _ in range(2000):
                    ix = np.concatenate([rng.choice(block, size=len(block), replace=True) for block in blocks])
                    boot.append(effect(y[ix], c[ix]))
            boot = np.asarray(boot).reshape((-1, 7, 6))
            for i, a in enumerate(labels.CD4):
                for j, b in enumerate(labels.CD8):
                    key = (context, phase, a, b)
                    keys.add(key)
                    r = indexed[key]
                    assert r['eligible_author_patients'] == n
                    for arm in r['eligible_by_arm']:
                        assert arm['patients'] == sum(z['neo_type'] == arm['arm'] for z, _, _ in ready)
                    equal(r['partial_Spearman'], point[i, j])
                    bv = boot[:, i, j]; bv = bv[np.isfinite(bv)]
                    lv = loo[:, i, j]; lv = lv[np.isfinite(lv)]
                    assert r['valid_bootstraps'] == len(bv) and r['valid_leave_one'] == len(lv)
                    ci_expected = np.quantile(bv, [.025, .975]) if len(bv) >= 1800 else [np.nan, np.nan]
                    for a0, b0 in zip(r['pointwise_bootstrap_CI'], ci_expected):
                        equal(a0, b0)
                    equal(r['leave_one_min'], lv.min() if len(lv) else np.nan)
                    equal(r['leave_one_max'], lv.max() if len(lv) else np.nan)
                    if np.isfinite(point[i, j]):
                        assert r['leave_one_same_direction'] == int(np.sum(np.sign(lv) == np.sign(point[i, j])))
                    else:
                        assert r['leave_one_same_direction'] is None
                    for k, v in sensitivity.items():
                        equal(r['sensitivity_partial_Spearman'][k], v[i, j])
                    expected_status = 'DESCRIPTIVE_REFERENCE_STATE_ASSOCIATION' if np.isfinite(point[i, j]) else ('INSUFFICIENT_PATIENTS' if n < 20 else 'UNIDENTIFIABLE_OR_INVARIANT')
                    assert r['status'] == expected_status
    assert keys == set(indexed)
    result = dict(status='ACCEPTED_ALL168_FINE_STATE_DESCRIPTIONS_FULL_BOOTSTRAP_LOO_SENSITIVITY', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(parent / 'private_terminal.json'), all_counts_and_patient_selection_reconstructed_with_Counter=True, all_executed_points_2000_streams_LOO_and_three_sensitivities_replayed_with_Pandas_rank_SVD=True, all_insufficient_opportunities_preserved=True, shares_pinned_reference_join_materializer_not_independent_RDS_decoder=True, no_clinical_spatial_or_population_truth_claim=True, no_individual_rows_exported=True)
    out.mkdir(mode=0o700)
    (out / 'private_independent_acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], acceptance_SHA256=sha(out / 'private_independent_acceptance.json'))))


if __name__ == '__main__':
    main()
