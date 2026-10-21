"""All reference-state combinations in source and naming-anchor contexts.

Outcome-free descriptive analysis; count references are not population truth,
CM02 naming anchors are not the full module, and neither context is replication.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs():
    parent = Path(os.environ['PANCANCER_PATIENT_SCOPE_OUTPUT'])
    t = json.loads((parent / 'private_terminal.json').read_text())
    a = Path(os.environ['PANCANCER_PATIENT_SCOPE_AUDIT'])
    audit = json.loads(a.read_text())
    assert audit['parent_terminal_SHA256'] == sha(parent / 'private_terminal.json')
    assert audit['all_patient_scope_fractions_replayed_with_Fraction']
    src = parent.parent / 'qualify_CM02_fine_state_patient_scope_v1.py'
    assert sha(src) == t['source_SHA256']
    spec = importlib.util.spec_from_file_location('qualified_patient_scope', src)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    tab, cur, unique, shared, matrices, table = module.load_inputs()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    catalog_path = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1/tables/subcluster_catalog_refined_private.tsv'
    catalog = pd.read_csv(catalog_path, sep='\t', usecols=['node_id', 'subcluster'], dtype=str)
    labels = json.loads(os.environ['PANCANCER_ANCHOR_LABELS'])
    assert len(labels) == 2
    nodes = []
    for label in labels:
        found = catalog.loc[catalog.subcluster.eq(label), 'node_id'].unique()
        assert len(found) == 1
        nodes.append(found[0])
    u = unique[((unique.reference == 'CD4_reference') & (unique.node_id == nodes[0])) | ((unique.reference == 'CD8_reference') & (unique.node_id == nodes[1]))]
    anchored = u.groupby(['library_id', 'label']).size().unstack(fill_value=0).reindex(index=tab.library_id, columns=module.CD4 + module.CD8, fill_value=0).astype(int)
    contexts = {'source_fine_reference': matrices['unique_current_reference'], 'CM02_naming_anchors_only': anchored}
    return module, tab, cur, contexts, unique, nodes, catalog_path


def effect(y, c, qr_check=True):
    empty = np.full((7, 6), np.nan)
    if len(y) < 20:
        return empty
    d = c.copy()
    for j in [1, 2]:
        d[:, j] = rankdata(d[:, j])
    d[:, 1:3] -= d[:, 1:3].mean(axis=0)
    d[:, 1:3] /= np.maximum(d[:, 1:3].std(axis=0), 1e-15)
    if np.linalg.matrix_rank(d) != 5 or len(y) - 6 < 10 or np.linalg.cond(d) > 1e6:
        return empty
    a = rankdata(y, axis=0)
    e = a - d @ np.linalg.lstsq(d, a, rcond=1e-12)[0]
    n = np.linalg.norm(e, axis=0)
    r = e[:, :7].T @ e[:, 7:] / np.maximum(n[:7, None] * n[None, 7:], 1e-15)
    r[(n[:7, None] < 1e-10) | (n[None, 7:] < 1e-10)] = np.nan
    if qr_check:
        q, _ = np.linalg.qr(d, mode='reduced')
        z = a - q @ (q.T @ a)
        v = np.linalg.norm(z, axis=0)
        alternative = z[:, :7].T @ z[:, 7:] / np.maximum(v[:7, None] * v[None, 7:], 1e-15)
        good = np.isfinite(r)
        np.testing.assert_allclose(r[good], alternative[good], atol=1e-10, rtol=1e-10)
    return r


def materialize_context(module, tab, cur, counts, phase):
    g = tab[tab.timepoint_registry.eq(phase)].copy()
    x = counts.loc[g.library_id].to_numpy(float)
    n4, n8 = x[:, :7].sum(axis=1), x[:, 7:].sum(axis=1)
    keep = (n4 >= 20) & (n8 >= 20) & ~g.patient_registry.duplicated(keep=False).to_numpy() & g.neo_type.isin(['D', 'D+T']).to_numpy() & g['HPV/p16'].isin(['positive', 'negative']).to_numpy()
    g = g.loc[keep].sort_values('patient_registry')
    x = counts.loc[g.library_id].to_numpy(float)
    n4, n8 = x[:, :7].sum(axis=1), x[:, 7:].sum(axis=1)
    y = np.column_stack([x[:, :7] / n4[:, None], x[:, 7:] / n8[:, None]])
    tx, hpv = g.neo_type.eq('D+T').to_numpy(int), g['HPV/p16'].eq('positive').to_numpy(int)
    c = np.column_stack([np.ones(len(g)), np.log1p(n4), np.log1p(n8), tx, hpv])
    denominators = {'all_unique_reference_T': n4 + n8, 'all_current_T': cur[cur.major.eq('T')].groupby('library_id').size().reindex(g.library_id, fill_value=0).to_numpy(float), 'all_current_source_cells': cur.groupby('library_id').size().reindex(g.library_id, fill_value=0).to_numpy(float)}
    assert all((d > 0).all() for d in denominators.values())
    return g, x, y, c, 2 * tx + hpv, denominators


def main():
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    save('plan_before_fine_state_associations.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), contexts=['source_fine_reference', 'CM02_naming_anchors_only'], all_positions=168, primary_phase='pre', secondary_phase='post_same_cohort', all7_CD4_by6_CD8_labels_preserved=True, main_denominators='within respective seven CD4/six CD8 unique author-reference state pools', minimum_parent_cells=20, minimum_patients=20, covariates=['log1p_CD4_reference_parent_count', 'log1p_CD8_reference_parent_count', 'treatment_DplusT', 'HPV_positive'], design='full five columns, condition below 1e6, at least10 residual df; never drop covariates', bootstrap='2000 patient draws within treatment x HPV; rerank every draw; 1800 finite required for pointwise CI', sensitivity_denominators=['all_unique_reference_T', 'all_current_T', 'all_current_source_cells'], seed=2026102103, no_P_FDR_or_confirmatory_pair_promotion=True, no_new_clinical_outcome_values_read=True, source_context_not_exact_module_validation=True, anchors_not_all_CM02_members=True, no_cell_or_patient_arrays_saved=True))
    module, tab, cur, contexts, unique, nodes, catalog = load_inputs()
    save('private_input_reference_proof.json', dict(catalog_SHA256=sha(catalog), anchor_labels=json.loads(os.environ['PANCANCER_ANCHOR_LABELS']), exact_anchor_node_ids=nodes, source_reference='GSE286827 processed author fine state subset', patient_scope_audit_SHA256=sha(Path(os.environ['PANCANCER_PATIENT_SCOPE_AUDIT']))))
    rows, proofs = [], []
    for ci, (context, counts) in enumerate(contexts.items()):
        for pi, phase in enumerate(['pre', 'post']):
            g, x, y, c, strata, denominators = materialize_context(module, tab, cur, counts, phase)
            n = len(g)
            main = effect(y, c)
            sensitivity = {k: effect(x / d[:, None], c) for k, d in denominators.items()}
            loo = np.array([effect(np.delete(y, i, 0), np.delete(c, i, 0)) for i in range(n)])
            rng = np.random.default_rng(2026102103 + 10 * ci + pi)
            blocks = [np.flatnonzero(strata == v) for v in np.unique(strata)]
            boot = []
            if np.isfinite(main).any():
                for _ in range(2000):
                    draw = np.concatenate([rng.choice(ix, len(ix), replace=True) for ix in blocks])
                    boot.append(effect(y[draw], c[draw]))
            boot = np.asarray(boot).reshape((-1, 7, 6))
            for i, a in enumerate(module.CD4):
                for j, b in enumerate(module.CD8):
                    point = float(main[i, j]) if np.isfinite(main[i, j]) else None
                    bv = boot[:, i, j]; bv = bv[np.isfinite(bv)]
                    lv = loo[:, i, j] if n else np.array([]); lv = lv[np.isfinite(lv)]
                    rows.append(dict(source='GSE286827', context=context, phase=phase, CD4_author_state=a, CD8_author_state=b, eligible_author_patients=n, eligible_by_arm=[dict(arm=arm, patients=int(g.neo_type.eq(arm).sum())) for arm in ['D', 'D+T']], partial_Spearman=point, pointwise_bootstrap_CI=np.quantile(bv, [.025, .975]).tolist() if len(bv) >= 1800 else [None, None], valid_bootstraps=len(bv), valid_leave_one=len(lv), leave_one_min=float(lv.min()) if len(lv) else None, leave_one_max=float(lv.max()) if len(lv) else None, leave_one_same_direction=int(np.sum(np.sign(lv) == np.sign(point))) if point is not None else None, sensitivity_partial_Spearman={k: float(v[i, j]) if np.isfinite(v[i, j]) else None for k, v in sensitivity.items()}, status='DESCRIPTIVE_REFERENCE_STATE_ASSOCIATION' if point is not None else ('INSUFFICIENT_PATIENTS' if n < 20 else 'UNIDENTIFIABLE_OR_INVARIANT'), original_labels_unchanged=True, no_P_or_multiplicity_controlled_claim=True))
            proofs.append(dict(context=context, phase=phase, eligible_patients=n, all_executed_point_bootstrap_LOO_sensitivity_QR_checked=True, no_effect_computed_below_patient_gate=n < 20))
    assert len(rows) == 168
    save('private_all168_fine_state_associations.json', rows)
    save('private_all4_context_numeric_proofs.json', proofs)
    terminal = dict(status='COMPLETE_ALL168_FINE_REFERENCE_AND_ANCHOR_STATE_DESCRIPTIONS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), no_cell_or_patient_arrays_saved=True, no_new_clinical_or_spatial_validation=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()})
    save('private_terminal.json', terminal)
    print(json.dumps(dict(status=terminal['status'], terminal_SHA256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
