"""Frozen conditional spatial RNA neighbourhood tests; individual cells stay private.

PANCANCER_PREFLIGHT_ROOT=<run> python scripts/spatial_conditional_neighborhood_v1.py
New output directory must contain predeclared_plan.json. No patient inference.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda: f.read(8388608), b''): h.update(x)
    return h.hexdigest()


def now(): return datetime.now(timezone.utc).isoformat()


def geometry(coords, area):
    radius = np.sqrt(area / np.pi)
    pairs = cKDTree(coords).query_pairs(6 * radius.max(), output_type='ndarray')
    if not len(pairs): return np.empty(0, int), np.empty(0, int), np.empty(0), np.zeros(len(area), int)
    d = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
    ratio = d / (radius[pairs[:, 0]] + radius[pairs[:, 1]])
    good = ratio <= 3
    pairs, ratio = pairs[good], ratio[good]
    ii = np.r_[pairs[:, 0], pairs[:, 1]]; jj = np.r_[pairs[:, 1], pairs[:, 0]]
    dd = np.r_[ratio, ratio]
    degree = np.bincount(ii[dd <= 2], minlength=len(area))
    assert (ii != jj).all()
    return ii, jj, dd, degree


def endpoints(z, row, parent_plan, mode):
    genes = list(z['genes']); lookup = {g: i for i, g in enumerate(genes)}
    def cols(s): return [lookup[g] for g in s]
    required = row['ligand'].split('_') + row['receptor'].split('_')
    names = [row['sender_family'].split('_')[-1], row['receiver_family'].split('_')[-1]]
    missing = [g for g in required + names if not z['present'][lookup[g]]]
    n = len(z['qc']); empty = np.zeros(n, bool)
    if missing: return empty, empty, empty, empty, 'MISSING_REQUIRED_GENE:' + ';'.join(sorted(set(missing)))
    pos = (z['counts'] > 0) & (z['counts'] > z['background_q95'][:, None])
    if mode == 'author_count_ge2': pos &= z['counts'] >= 2
    sline = 'T' if row['sender_family'].startswith('T_') else 'B'
    rline = 'T' if row['receiver_family'].startswith('T_') else 'B'
    if mode.startswith('author'):
        sid = row['sample_id']
        if sid.startswith('GSM8746298'): return empty, empty, empty, empty, 'AUTHOR_NOT_RETRIEVED'
        if sid.startswith(('GSM9046088', 'GSM9046089')) and 'B' in [sline, rline]:
            return empty, empty, empty, empty, 'AUTHOR_B_CATEGORY_UNAVAILABLE'
        lin = {k: np.isin(z['author_labels'], parent_plan['author_' + k + '_labels']) for k in ['T', 'B']}
        purity = {'T': np.ones(n, bool), 'B': np.ones(n, bool)}
    else:
        forbidden = set(required + parent_plan['family_name_genes'])
        panels = {k: [g for g in pp if g not in forbidden and z['present'][lookup[g]]] for k, pp in parent_plan['panels'].items()}
        if len(panels[sline]) < 2 or len(panels[rline]) < 2:
            return empty, empty, empty, empty, 'LINEAGE_PANEL_UNMEASURABLE'
        lin = {k: pos[:, cols(pp)].sum(axis=1) >= 2 for k, pp in panels.items()}
        purity = {k: ~np.logical_or.reduce([lin[j] for j in lin if j != k]) if mode == 'RNA_pure' else np.ones(n, bool) for k in ['T', 'B']}
    ss = z['qc'] & lin[sline] & purity[sline] & pos[:, lookup[names[0]]]
    rr = z['qc'] & lin[rline] & purity[rline] & pos[:, lookup[names[1]]]
    ll = ss & pos[:, cols(row['ligand'].split('_'))].all(axis=1)
    rc = rr & pos[:, cols(row['receptor'].split('_'))].all(axis=1)
    return ss, rr, ll, rc, 'MEASURED'


def conditional_draws(codes, groups, ii, jj, distance, seed, permutations):
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    draws = np.empty((permutations, len(codes)), np.uint8)
    strata = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    for b in range(permutations):
        for ix in strata: draws[b, ix] = codes[rng.permutation(ix)]
    score = np.zeros((3, permutations), np.int64)
    for k, rad in enumerate([1, 2, 3]):
        x, y = ii[distance <= rad], jj[distance <= rad]
        for begin in range(0, permutations, 100):
            dd = draws[begin:begin + 100]
            score[k, begin:begin + len(dd)] = (((dd[:, x] & 1) > 0) & ((dd[:, y] & 2) > 0)).sum(axis=1)
    movable = sum(len(ix) for ix in strata if len(np.unique(codes[ix])) > 1)
    return score, movable, len(strata)


def stats(obs, draws, reason):
    m = float(draws.mean()); sd = float(draws.std(ddof=1))
    status = reason if reason != 'MEASURED' else ('ZERO_NULL_VARIANCE' if sd == 0 else 'ESTIMABLE_CONDITIONAL')
    valid = status == 'ESTIMABLE_CONDITIONAL'
    up = (1 + int((draws >= obs).sum())) / (len(draws) + 1) if valid else None
    lo = (1 + int((draws <= obs).sum())) / (len(draws) + 1) if valid else None
    measured = reason in ['MEASURED', 'NO_POSITIVE_ENDPOINT']
    return dict(status=status, observed_pairs=int(obs) if measured else None, null_mean=m if measured else None,
                null_sd=sd if measured else None, null_q025=float(np.quantile(draws, .025)) if measured else None,
                null_q975=float(np.quantile(draws, .975)) if measured else None,
                excess_pairs=float(obs - m) if measured else None,
                enrichment_ratio=float(obs / m) if measured and m > 0 else None,
                p_enrichment=up, p_depletion=lo, p_two_sided=min(1, 2 * min(up, lo)) if valid else None)


def main():
    start = now(); root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    p = root / 'immune_same_name_family_source_review_v1'; old = p / 'spatial_contact_native_preflight_v3'
    out = p / 'spatial_contact_conditional_neighborhood_v1'
    assert out.is_dir() and not (out / 'inventory.json').exists()
    plan = json.loads((out / 'predeclared_plan.json').read_text())
    assert sha(old / 'inventory.json') == plan['parent']['inventory_sha256']
    assert sha(old / 'independent_audit.json') == plan['parent']['audit_sha256']
    assert sha(p / 'axis_validation_feasibility_v1/automatic_next_validation_records.tsv') == plan['parent']['candidate_sha256']
    inv = json.loads((old / 'inventory.json').read_text()); pp = json.loads((old / 'predeclared_plan.json').read_text())
    assert inv['cells'] == 254458 and len(inv['inputs']) == 5
    assert plan['null']['permutations'] == 2000 and plan['null']['fixed_seed'] == 20261011
    assert sha(old / 'all_contact_record_testability.tsv') == inv['outputs']['all_contact_record_testability.tsv']
    records = pd.read_csv(old / 'all_contact_record_testability.tsv', sep='\t')
    keys = ['sender_family', 'receiver_family', 'ligand', 'receptor']
    interfaces = records[keys].drop_duplicates().sort_values(keys).reset_index(drop=True)
    assert len(records) == 65 and len(interfaces) == 9
    interfaces['interface_id'] = ['LR%03d' % (i + 1) for i in range(9)]
    mapping = records[['model_key', 'sample_id', 'dataset_id', 'discovery_cancer', 'spatial_cancer', 'program'] + keys].merge(interfaces, on=keys, validate='many_to_one')
    mapping.to_csv(out / 'all65_record_mapping.tsv', sep='\t', index=False)
    private = out / 'private_null_cache'; private.mkdir(exist_ok=False)
    all_rows, radii_rows, sensitivity, fov_rows, source_hashes = [], [], [], [], {}
    for fi, inp in enumerate(sorted(inv['inputs'], key=lambda x: x['sample_id'])):
        sid = inp['sample_id']; path = old / 'private_native_cell_caches' / (sid + '.npz')
        rel = str(path.relative_to(old)); assert sha(path) == inv['outputs'][rel]
        source_hashes[sid] = sha(path)
        z = dict(np.load(path, allow_pickle=False)); fovs = sorted(set(z['fov']), key=int)
        graph = {}
        for vi, fv in enumerate(fovs):
            ix = np.flatnonzero((z['fov'] == fv) & z['qc'])
            if not len(ix): graph[fv] = (ix, np.empty(0, int), np.empty(0, int), np.empty(0), np.empty(0, int)); continue
            ii, jj, dist, dens = geometry(z['coords'][ix] / .12028, z['area'][ix])
            graph[fv] = ix, ii, jj, dist, dens
        for ui, interface in interfaces.iterrows():
            row = dict(interface); row['sample_id'] = sid
            ss, rr, ll, rc, reason = endpoints(z, row, pp, 'author_background')
            if reason == 'MEASURED' and (not ll.any() or not rc.any()): reason = 'NO_POSITIVE_ENDPOINT'
            total_obs = np.zeros(3, np.int64); total_null = np.zeros((3, 2000), np.int64)
            byfov_obs, byfov_null, moving, family_n = [], [], 0, int((ss | rr).sum())
            for vi, fv in enumerate(fovs):
                ix, ii, jj, dist, dens = graph[fv]
                observed = np.array([int((ll[ix[ii[dist <= rad]]] & rc[ix[jj[dist <= rad]]]).sum()) for rad in [1, 2, 3]])
                total_obs += observed
                union = ss[ix] | rr[ix]
                loc = np.flatnonzero(union); lookup = np.full(len(ix), -1, int); lookup[loc] = np.arange(len(loc))
                keep = union[ii] & union[jj]; ni, nj, nd = lookup[ii[keep]], lookup[jj[keep]], dist[keep]
                code = ll[ix[loc]].astype(np.uint8) + 2 * rc[ix[loc]].astype(np.uint8)
                membership = ss[ix[loc]].astype(np.int64) + 2 * rr[ix[loc]].astype(np.int64)
                depth = (z['native_total'][ix[loc]] > np.median(z['native_total'][ix])).astype(int) if len(ix) else np.empty(0, int)
                density = (dens[loc] > np.median(dens)).astype(int) if len(ix) else np.empty(0, int)
                groups = membership * 4 + depth * 2 + density
                draws, movable, strata = conditional_draws(code, groups, ni, nj, nd, [20261011, fi, ui, vi], 2000)
                total_null += draws; moving += movable; byfov_obs.append(observed); byfov_null.append(draws)
                fov_rows.append(dict(sample_id=sid, interface_id=row['interface_id'], fov=fv, QC_cells=len(ix), family_union=len(loc), sender=int(ll[ix].sum()), receiver=int(rc[ix].sum()), overlap=int((ll[ix] & rc[ix]).sum()), strata=strata, movable_cells=movable, observed_primary=int(observed[1]), null_mean=float(draws[1].mean()), primary_variance=float(draws[1].var())))
            byfov_obs = np.asarray(byfov_obs); byfov_null = np.asarray(byfov_null)
            np.savez_compressed(private / (sid + '__' + row['interface_id'] + '.npz'), observed=byfov_obs, null=byfov_null, fovs=np.asarray(fovs), reason=np.asarray(reason))
            base = dict(interface_id=row['interface_id'], sample_id=sid, dataset_id=records.loc[records.sample_id == sid, 'dataset_id'].iloc[0], spatial_cancer=records.loc[records.sample_id == sid, 'spatial_cancer'].iloc[0], **{k: row[k] for k in keys}, sender_family_cells=int(ss.sum()) if reason in ['MEASURED', 'NO_POSITIVE_ENDPOINT'] else None, receiver_family_cells=int(rr.sum()) if reason in ['MEASURED', 'NO_POSITIVE_ENDPOINT'] else None, sender_positive=int(ll.sum()) if reason in ['MEASURED', 'NO_POSITIVE_ENDPOINT'] else None, receiver_positive=int(rc.sum()) if reason in ['MEASURED', 'NO_POSITIVE_ENDPOINT'] else None, overlap_same_cell=int((ll & rc).sum()) if reason in ['MEASURED', 'NO_POSITIVE_ENDPOINT'] else None, family_union=family_n, movable_cells=moving, FOVs=len(fovs), paired_positive_FOVs=sum(x['sender'] > 0 and x['receiver'] > 0 for x in fov_rows[-len(fovs):]), patient_level_inference=False)
            for k, rad in enumerate([1, 2, 3]): radii_rows.append(dict(base, radius_multiplier=rad, **stats(total_obs[k], total_null[k], reason)))
            s = stats(total_obs[1], total_null[1], reason)
            if s['status'] == 'ESTIMABLE_CONDITIONAL':
                leave = (total_obs[1] - byfov_obs[:, 1]) - (total_null[1].mean() - byfov_null[:, 1].mean(axis=1))
                s['leave_one_FOV_excess_min'] = float(leave.min()); s['leave_one_FOV_excess_max'] = float(leave.max())
                s['max_FOV_pair_fraction'] = float(byfov_obs[:, 1].max() / total_obs[1]) if total_obs[1] > 0 else None
            else: s.update(leave_one_FOV_excess_min=None, leave_one_FOV_excess_max=None, max_FOV_pair_fraction=None)
            all_rows.append(dict(base, **s))
            for mode in ['RNA_background', 'RNA_pure', 'author_count_ge2']:
                a, b, ls, rs, why = endpoints(z, row, pp, mode)
                obs = sum(int((ls[ix[ii[dist <= 2]]] & rs[ix[jj[dist <= 2]]]).sum()) for ix, ii, jj, dist, dens in graph.values())
                sensitivity.append(dict(interface_id=row['interface_id'], sample_id=sid, mode=mode, status=why, sender_positive=int(ls.sum()) if why == 'MEASURED' else None, receiver_positive=int(rs.sum()) if why == 'MEASURED' else None, overlap_same_cell=int((ls & rs).sum()) if why == 'MEASURED' else None, observed_pairs=obs if why == 'MEASURED' else None, pvalue_not_computed=True))
            print(json.dumps(dict(completed_interface=row['interface_id'], sample_id=sid, status=s['status'], observed=s['observed_pairs'], p=s['p_enrichment'], at_UTC=now())), flush=True)
    result = pd.DataFrame(all_rows); assert len(result) == 45
    raw = result.p_enrichment.fillna(1).to_numpy(); order = np.argsort(raw)
    adj = np.minimum.accumulate((raw[order] * 45 / np.arange(1, 46))[::-1])[::-1].clip(0, 1)
    q = np.empty(45); q[order] = adj; result['BH_full45'] = np.where(result.p_enrichment.notna(), q, np.nan)
    result.to_csv(out / 'all45_interface_results.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(radii_rows).to_csv(out / 'all135_radius_results.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(sensitivity).to_csv(out / 'all135_endpoint_sensitivity.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(fov_rows).to_csv(out / 'private_FOV_diagnostics.tsv', sep='\t', index=False, na_rep='NA')
    outputs = {str(q.relative_to(out)): sha(q) for q in out.rglob('*') if q.is_file()}
    inventory = dict(schema='conditional_spatial_neighborhood_inventory_v1', status='PRODUCER_COMPLETE_PENDING_INDEPENDENT_AUDIT', started_UTC=start, finished_UTC=now(), producer_sha256=sha(__file__), plan_sha256=sha(out / 'predeclared_plan.json'), private_input_cache_sha256=source_hashes, outputs=outputs, tests=45, mapped_records=65, sensitivity_rows=135, radius_rows=135, permutations=2000, patient_level_inference=False, axis_validation_complete=False, scientific_release=False)
    (out / 'inventory.json').write_text(json.dumps(inventory, indent=2) + '\n')
    print(json.dumps(dict(status=inventory['status'], tests=45, significant_primary=int((result.BH_full45 < .05).sum()), estimable=int(result.status.eq('ESTIMABLE_CONDITIONAL').sum()), finished_UTC=inventory['finished_UTC'])), flush=True)


if __name__ == '__main__': main()
