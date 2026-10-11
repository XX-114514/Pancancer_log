"""Independent cached-native endpoint, graph and all45 conditional-null replay.

No producer import. Native-cache values were independently audited in parent v3.
Uses query_ball_point distance checks and CSR multiplication, versus query_pairs
and edge-wise sums in production. Random draws share the predeclared seed contract.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix


def digest(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''): h.update(b)
    return h.hexdigest()


def stamp(): return datetime.now(timezone.utc).isoformat()


def independent_graph(xy, area):
    r = np.sqrt(area / np.pi); n = len(r)
    if n == 0: return [csr_matrix((0, 0), dtype=np.int64)] * 3
    tree = cKDTree(xy); domains = tree.query_ball_point(xy, 6 * max(r))
    rows, cols, ratios = [], [], []
    for i, domain in enumerate(domains):
        j = np.asarray(domain, int); j = j[j != i]
        dd = np.sqrt(((xy[j] - xy[i]) ** 2).sum(axis=1)) / (r[i] + r[j])
        good = dd <= 3
        rows.extend([i] * int(good.sum())); cols.extend(j[good]); ratios.extend(dd[good])
    rows, cols, ratios = np.asarray(rows, int), np.asarray(cols, int), np.asarray(ratios)
    return [csr_matrix((np.ones(int((ratios <= k).sum()), np.int64), (rows[ratios <= k], cols[ratios <= k])), shape=(n, n)) for k in [1, 2, 3]]


def fixture_check():
    xy = np.array([[0., 0.], [1., 0.], [5., 0.], [1., 0.]])
    area = np.full(4, np.pi)
    for k, a in zip([1, 2, 3], independent_graph(xy, area)):
        brute = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=2) <= k * 2
        np.fill_diagonal(brute, False)
        assert np.array_equal(a.toarray(), brute)
        assert not a.diagonal().any() and (a != a.T).nnz == 0
    # Overlapping endpoint at cell0 cannot create a self-contact; coincident
    # distinct cells are allowed and must be kept rather than silently dropped.
    l = np.array([1, 0, 0, 0]); r = np.array([1, 0, 0, 0])
    assert l @ independent_graph(xy, area)[1] @ r == 0


def masks(z, row, pp, mode):
    gene_index = dict(zip(z['genes'], range(len(z['genes']))))
    def ix(genes): return [gene_index[g] for g in genes]
    lig = row.ligand.split('_'); rec = row.receptor.split('_')
    nm = [row.sender_family.rsplit('_', 1)[1], row.receiver_family.rsplit('_', 1)[1]]
    missing = sorted(set(g for g in lig + rec + nm if not z['present'][gene_index[g]]))
    zero = np.zeros(len(z['qc']), bool)
    if missing: return (zero,) * 4 + ('MISSING_REQUIRED_GENE:' + ';'.join(missing),)
    y = z['counts'] > np.maximum(0, z['background_q95'][:, None])
    if mode == 'author_count_ge2': y = y & (z['counts'] >= 2)
    sk = 'T' if row.sender_family.startswith('T_') else 'B'; rk = 'T' if row.receiver_family.startswith('T_') else 'B'
    if mode.startswith('author'):
        if row.sample_id.startswith('GSM8746298'): return (zero,) * 4 + ('AUTHOR_NOT_RETRIEVED',)
        if row.sample_id.startswith(('GSM9046088', 'GSM9046089')) and 'B' in [sk, rk]:
            return (zero,) * 4 + ('AUTHOR_B_CATEGORY_UNAVAILABLE',)
        lineage = {k: np.array([s in pp['author_' + k + '_labels'] for s in z['author_labels']]) for k in ['T', 'B']}
    else:
        banned = set(lig + rec + pp['family_name_genes'])
        panels = {k: [g for g in gs if g not in banned and z['present'][gene_index[g]]] for k, gs in pp['panels'].items()}
        if min(len(panels[sk]), len(panels[rk])) < 2: return (zero,) * 4 + ('LINEAGE_PANEL_UNMEASURABLE',)
        lineage = {k: np.count_nonzero(y[:, ix(gs)], axis=1) >= 2 for k, gs in panels.items()}
        if mode == 'RNA_pure':
            original = dict(lineage)
            for k in [sk, rk]: lineage[k] = original[k] & ~np.column_stack([v for j, v in original.items() if j != k]).any(axis=1)
    s = lineage[sk] & z['qc'] & y[:, gene_index[nm[0]]]
    r = lineage[rk] & z['qc'] & y[:, gene_index[nm[1]]]
    l = s & y[:, ix(lig)].all(axis=1); q = r & y[:, ix(rec)].all(axis=1)
    return s, r, l, q, 'MEASURED'


def independent_permutations(code, groups, adjacency, seed):
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    strata = [np.where(groups == value)[0] for value in sorted(set(groups))]
    draws = np.zeros((2000, len(code)), np.uint8)
    for b in range(2000):
        for members in strata:
            order = rng.permutation(members)
            draws[b, members] = np.take(code, order)
    # Every draw must conserve the full joint two-bit histogram within each bin,
    # including overlapping family and same-cell positive expression.
    for members in strata:
        for val in range(4): assert np.all((draws[:, members] == val).sum(axis=1) == (code[members] == val).sum())
    left = (draws & 1).astype(np.int64); right = (draws >> 1).astype(np.int64)
    scores = np.array([np.sum(left * (a @ right.T).T, axis=1) for a in adjacency])
    move = sum(len(x) for x in strata if len(set(code[x])) > 1)
    return scores, move


def main():
    start = stamp(); fixture_check()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); p = root / 'immune_same_name_family_source_review_v1'
    old = p / 'spatial_contact_native_preflight_v3'; out = p / 'spatial_contact_conditional_neighborhood_v1'
    assert not (out / 'independent_audit.json').exists()
    inv = json.loads((out / 'inventory.json').read_text()); plan = json.loads((out / 'predeclared_plan.json').read_text())
    assert digest(out / 'predeclared_plan.json') == inv['plan_sha256']
    assert digest(old / 'inventory.json') == plan['parent']['inventory_sha256']
    assert digest(old / 'independent_audit.json') == plan['parent']['audit_sha256']
    for name, h in inv['outputs'].items(): assert digest(out / name) == h
    pp = json.loads((old / 'predeclared_plan.json').read_text())
    previous = pd.read_csv(old / 'all_contact_record_testability.tsv', sep='\t')
    result = pd.read_csv(out / 'all45_interface_results.tsv', sep='\t')
    radii = pd.read_csv(out / 'all135_radius_results.tsv', sep='\t')
    sens = pd.read_csv(out / 'all135_endpoint_sensitivity.tsv', sep='\t')
    mapping = pd.read_csv(out / 'all65_record_mapping.tsv', sep='\t')
    axes = previous[['sender_family', 'receiver_family', 'ligand', 'receptor']].drop_duplicates().sort_values(['sender_family', 'receiver_family', 'ligand', 'receptor']).reset_index(drop=True)
    assert len(axes) == 9 and len(result) == 45 and len(mapping) == 65
    assert mapping[['sample_id', 'model_key']].equals(previous[['sample_id', 'model_key']])
    parent_inv = json.loads((old / 'inventory.json').read_text()); checks = []; cells = 0; replayed = 0
    for fi, inp in enumerate(sorted(parent_inv['inputs'], key=lambda v: v['sample_id'])):
        sid = inp['sample_id']; path = old / 'private_native_cell_caches' / (sid + '.npz')
        assert digest(path) == inv['private_input_cache_sha256'][sid] == parent_inv['outputs'][str(path.relative_to(old))]
        z = dict(np.load(path, allow_pickle=False)); cells += len(z['qc']); fovs = sorted(set(z['fov']), key=int)
        domain = {}
        for fv in fovs:
            idx = np.where((z['fov'] == fv) & z['qc'])[0]
            a = independent_graph(z['coords'][idx] / .12028, z['area'][idx])
            for graph in a: assert not graph.diagonal().any() and (graph != graph.T).nnz == 0
            # Direct full distance matrix for one FOV/file independently proves
            # the area-normalized bound beyond the small synthetic fixture.
            if fv == fovs[0] and len(idx) <= 5000:
                x = z['coords'][idx] / .12028; r = np.sqrt(z['area'][idx] / np.pi)
                dd = np.sqrt(((x[:, None] - x[None, :]) ** 2).sum(axis=2)) / (r[:, None] + r[None, :])
                for k, graph in zip([1, 2, 3], a):
                    direct = dd <= k; np.fill_diagonal(direct, False); assert np.array_equal(direct, graph.toarray())
            domain[fv] = idx, a
        for ui, axis in axes.iterrows():
            row = axis.copy(); row['sample_id'] = sid; iid = 'LR%03d' % (ui + 1)
            target = result[(result.sample_id == sid) & (result.interface_id == iid)].iloc[0]
            s, r, l, q, why = masks(z, row, pp, 'author_background')
            if why == 'MEASURED' and (l.sum() == 0 or q.sum() == 0): why = 'NO_POSITIVE_ENDPOINT'
            obs, null, movable = [], [], 0
            for vi, fv in enumerate(fovs):
                idx, a = domain[fv]; keep = np.where(s[idx] | r[idx])[0]
                graphs = [g[keep][:, keep] for g in a]
                labels = l[idx[keep]].astype(np.uint8) + (q[idx[keep]].astype(np.uint8) << 1)
                family = s[idx[keep]].astype(int) + 2 * r[idx[keep]].astype(int)
                counts = z['native_total'][idx]; density = np.asarray(a[1].sum(axis=1)).ravel()
                depth = counts[keep] > np.median(counts) if len(idx) else np.empty(0, bool)
                crowded = density[keep] > np.median(density) if len(idx) else np.empty(0, bool)
                groups = 4 * family + 2 * depth + crowded
                observed = np.array([l[idx].astype(np.int64) @ g @ q[idx].astype(np.int64) for g in a])
                simulated, move = independent_permutations(labels, groups, graphs, [20261011, fi, ui, vi])
                obs.append(observed); null.append(simulated); movable += move; replayed += 2000
            obs, null = np.asarray(obs), np.asarray(null)
            cached = np.load(out / 'private_null_cache' / (sid + '__' + iid + '.npz'), allow_pickle=False)
            assert np.array_equal(obs, cached['observed']) and np.array_equal(null, cached['null'])
            assert str(cached['reason']) == why and np.array_equal(cached['fovs'], np.asarray(fovs))
            total_obs, total_null = obs.sum(axis=0), null.sum(axis=0)
            measured = why in ['MEASURED', 'NO_POSITIVE_ENDPOINT']
            if measured:
                for name, val in [('sender_positive', l.sum()), ('receiver_positive', q.sum()), ('overlap_same_cell', (l & q).sum()), ('sender_family_cells', s.sum()), ('receiver_family_cells', r.sum()), ('movable_cells', movable)]: assert target[name] == val
            for k, radius in enumerate([1, 2, 3]):
                rr = radii[(radii.sample_id == sid) & (radii.interface_id == iid) & (radii.radius_multiplier == radius)].iloc[0]
                if not measured:
                    assert pd.isna(rr.observed_pairs) and pd.isna(rr.p_enrichment); continue
                assert rr.observed_pairs == total_obs[k] and np.isclose(rr.null_mean, total_null[k].mean(), atol=1e-10)
                sd = total_null[k].std(ddof=1)
                if why == 'MEASURED' and sd > 0:
                    assert rr.status == 'ESTIMABLE_CONDITIONAL'
                    expect = (1 + (total_null[k] >= total_obs[k]).sum()) / 2001
                    assert np.isclose(rr.p_enrichment, expect, rtol=0, atol=1e-14)
                else: assert pd.isna(rr.p_enrichment)
            for mode in ['RNA_background', 'RNA_pure', 'author_count_ge2']:
                ss, rs, ls, rq, reason = masks(z, row, pp, mode)
                rr = sens[(sens.sample_id == sid) & (sens.interface_id == iid) & (sens['mode'] == mode)].iloc[0]
                assert rr.status == reason
                if reason == 'MEASURED':
                    counted = sum(ls[ix].astype(np.int64) @ aa[1] @ rq[ix].astype(np.int64) for ix, aa in domain.values())
                    assert rr.observed_pairs == counted and rr.sender_positive == ls.sum() and rr.receiver_positive == rq.sum()
                else: assert pd.isna(rr.observed_pairs)
            checks.append(dict(interface_id=iid, sample_id=sid, status='PASS_ALL_COUNTS_GEOMETRY_AND2000_CONDITIONAL_DRAWS', endpoint_status=why, observed_primary=int(total_obs[1]) if measured else None))
            print(json.dumps(dict(audited_interface=iid, sample_id=sid, at_UTC=stamp())), flush=True)
    assert cells == 254458 and len(checks) == 45
    pvals = result.p_enrichment.fillna(1).to_numpy(); ranks = np.argsort(pvals)
    corrected = np.minimum.accumulate((45 * pvals[ranks] / np.arange(1, 46))[::-1])[::-1].clip(0, 1)
    expected = np.empty(45); expected[ranks] = corrected
    assert np.allclose(expected[result.p_enrichment.notna()], result.loc[result.p_enrichment.notna(), 'BH_full45'])
    assert result.loc[result.p_enrichment.isna(), 'BH_full45'].isna().all()
    pd.DataFrame(checks).to_csv(out / 'independent_all45_checks.tsv', sep='\t', index=False, na_rep='NA')
    audit = dict(schema='conditional_spatial_neighborhood_audit_v1', status='PASS_ALL45_NATIVE_ENDPOINT_GEOMETRY_JOINT_NULL_AND_FULL_FAMILY_BH', started_UTC=start, finished_UTC=stamp(), cells=cells, tests=45, mapped_records=65, radius_checks=135, sensitivity_checks=135, FOV_permutation_draws_replayed=replayed, producer_inventory_sha256=digest(out / 'inventory.json'), auditor_sha256=digest(__file__), fixture='brute force distances, overlap self-edge, coincident distinct cells', joint_histograms_all2000_draws_verified=True, geometric_self_edges=0, between_FOV_edges=0, patient_mapping_validated=False, independence='No producer import; separately constructed query-ball graph, CSR bilinear counts, independent endpoint parsing. Parent v3 native cache accepted; identical null specification and RNG seed deliberately replayed. This is computational verification, not biological replication.', scientific_release=False, axis_validation_complete=False)
    (out / 'independent_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit), flush=True)


if __name__ == '__main__': main()
