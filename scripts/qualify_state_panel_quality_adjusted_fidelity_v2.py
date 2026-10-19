"""Outcome-blind technical transfer triage; original full-signature gates stay closed."""
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem
from scipy.stats import rankdata


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def weighted_center(a, w):
    return a - np.sum(w[:, None] * a, axis=0)


def correlations(a, b, w):
    den = np.sqrt(np.sum(w[:, None]*a*a, axis=0)*np.sum(w[:, None]*b*b, axis=0))
    return np.divide(np.sum(w[:, None]*a*b, axis=0), den,
                     out=np.full(a.shape[1], np.nan), where=den > 1e-14)


def fit_scores(p, o, quality, blocks):
    """Equal block weight; each cell within a block receives equal share."""
    _, inv, counts = np.unique(blocks, return_inverse=True, return_counts=True)
    w = 1.0 / (len(counts)*counts[inv])
    q = rankdata(quality, axis=0, method='average')
    q = weighted_center(q, w)
    scale = np.sqrt(np.sum(w[:, None]*q*q, axis=0))
    if np.any(scale <= 1e-12):
        return None, 'CONSTANT_REQUIRED_QUALITY_COVARIATE'
    q /= scale
    X = np.column_stack([np.ones(len(w)), q])
    sw = np.sqrt(w)
    Q, R = np.linalg.qr(sw[:, None]*X, mode='reduced')
    if np.linalg.matrix_rank(R) < 3 or np.linalg.cond(R) >= 1e8:
        return None, 'REQUIRED_QUALITY_DESIGN_UNIDENTIFIABLE'
    a = rankdata(p, axis=0, method='average')
    b = rankdata(o, axis=0, method='average')
    ac = weighted_center(a, w); bc = weighted_center(b, w)
    raw = correlations(ac, bc, w)
    ar = (sw[:, None]*a - Q@(Q.T@(sw[:, None]*a)))/sw[:, None]
    br = (sw[:, None]*b - Q@(Q.T@(sw[:, None]*b)))/sw[:, None]
    adjusted = correlations(ar, br, w)
    # Algebraically independent Schur-complement check of every cell-level value.
    G = q.T@(w[:, None]*q)
    ap = q.T@(w[:, None]*ac); bp = q.T@(w[:, None]*bc)
    cov = np.sum(w[:, None]*ac*bc, axis=0)-np.sum(ap*np.linalg.solve(G, bp), axis=0)
    va = np.sum(w[:, None]*ac*ac, axis=0)-np.sum(ap*np.linalg.solve(G, ap), axis=0)
    vb = np.sum(w[:, None]*bc*bc, axis=0)-np.sum(bp*np.linalg.solve(G, bp), axis=0)
    den = np.sqrt(np.maximum(va, 0)*np.maximum(vb, 0))
    check = np.divide(cov, den, out=np.full_like(cov, np.nan), where=den > 1e-14)
    # Near-exactly explained scores are unidentifiable under either calculation.
    original_va = np.sum(w[:, None]*ac*ac, axis=0)
    original_vb = np.sum(w[:, None]*bc*bc, axis=0)
    keep = (original_va > 1e-12)&(original_vb > 1e-12)&(va > original_va*1e-10)&(vb > original_vb*1e-10)
    adjusted[~keep] = np.nan; check[~keep] = np.nan
    np.testing.assert_allclose(adjusted, check, rtol=1e-7, atol=2e-8, equal_nan=True)
    ba = np.stack([ar[inv == j].mean(0) for j in range(len(counts))])
    bb = np.stack([br[inv == j].mean(0) for j in range(len(counts))])
    block_vary = (np.var(ba, axis=0) > np.sum(w[:, None]*ar*ar, axis=0)*1e-10)&(np.var(bb, axis=0) > np.sum(w[:, None]*br*br, axis=0)*1e-10)
    ba = rankdata(ba, axis=0, method='average'); bb = rankdata(bb, axis=0, method='average')
    bw = np.full(len(counts), 1/len(counts))
    ba = weighted_center(ba, bw); bb = weighted_center(bb, bw)
    block = correlations(ba, bb, bw); block[(~keep)|(~block_vary)] = np.nan
    dep = [correlations(ac, np.repeat(q[:, j:j+1], ac.shape[1], axis=1), w) for j in range(2)]
    err = np.abs(adjusted-check)
    return dict(raw=raw, adjusted=adjusted, block=block, RNA=dep[0], detected=dep[1],
                max_schur_error=float(np.nanmax(err)) if np.isfinite(err).any() else 0), 'ESTIMABLE_TECHNICAL_DESCRIPTION'


def main():
    V = Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    planp = Path(os.environ['PANCANCER_QUALITY_PROXY_PLAN']); plan = json.loads(planp.read_text())
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    assert shutil.disk_usage('/dev/shm').free > 512*1024**2
    assert sha(__file__) == plan['producer_sha256']
    dirs = {k: V/v['directory'] for k, v in plan['inputs'].items()}
    manifests = {}
    for k, p in dirs.items():
        target = p/plan['inputs'][k]['manifest']
        assert sha(target) == plan['inputs'][k]['sha256']
        manifests[k] = json.loads(target.read_text())
    assert manifests['audit']['calibration_inventory_sha256'] == plan['inputs']['scores']['sha256']
    dictionary = json.loads((dirs['dictionary']/plan['inputs']['dictionary']['manifest']).read_text())['programmes']
    recs = {r['programme_id']:r for r in dictionary}; assert len(recs) == 814
    prior = pd.read_csv(dirs['scores']/'public_all_source_programme_axis_summary.tsv', sep='\t')
    assert sha(dirs['scores']/'public_all_source_programme_axis_summary.tsv') == manifests['scores']['public_outputs']['public_all_source_programme_axis_summary.tsv']
    lookup = prior.set_index(['gse_id','major','programme_id','axis_id'])
    assert lookup.index.is_unique
    axes = sorted(prior.axis_id.unique()); assert len(axes) == 4
    score_files = {r['relative_file']:r for r in manifests['scores']['private_outputs'] if r['relative_file'].endswith('.npz')}
    assert len(score_files) == 2836
    out.mkdir(mode=0o700)
    (out/'plan_before_quality_adjusted_scores.json').write_text(json.dumps(plan, indent=2)+'\n')
    rows = []; processed = 0; checked = 0; maxerr = 0.; metadata_cells = 0
    for ordinal, info in enumerate(manifests['native']['private_native_files']):
        p = dirs['native']/info['relative_file']; assert sha(p) == info['sha256']
        with h5py.File(p, 'r') as h:
            obs = pd.DataFrame({k:read_elem(h['obs'][k]) for k in ['core_position','unit_index','full_RNA_total','full_detected_genes','node_id','stratum','effective_resampling_block']})
        n = len(obs); assert n == info['cells']; metadata_cells += n
        groups = pd.factorize(pd.MultiIndex.from_frame(obs[['node_id','stratum']]), sort=True)[0]
        target_groups = [np.flatnonzero(groups == k) for k in np.unique(groups)]
        for axis in axes:
            rel = 'private_cell_scores/'+Path(info['relative_file']).stem+'_'+axis[:10]+'.npz'
            si = score_files[rel]; p = dirs['scores']/rel; assert sha(p) == si['sha256']
            with np.load(p, allow_pickle=False) as z:
                np.testing.assert_array_equal(z['global_rows'], obs.core_position)
                np.testing.assert_array_equal(z['unit_index'], obs.unit_index)
                ids = z['programme_ids'].tolist(); proxy = z['proxy'].astype(float); outside = z['outside_reference'].astype(float)
            assert len(ids) == info['programmes'] and proxy.shape == outside.shape == (n, len(ids))
            for stratum, ix in enumerate(target_groups):
                blocks = obs.effective_resampling_block.iloc[ix].astype(str).to_numpy()
                nb = len(np.unique(blocks)); common = dict(gse_id=info['gse_id'], major=info['major'], axis_id=axis,
                                                          native_file_ordinal=ordinal, identity_stratum_ordinal=stratum,
                                                          cells=len(ix), analysis_blocks=nb)
                take = np.isfinite(proxy[ix]).all(0)&np.isfinite(outside[ix]).all(0)
                status = 'FEWER_THAN50CELLS_OR10BLOCKS' if len(ix) < 50 or nb < 10 else 'READY'
                result = None
                if status == 'READY' and take.any():
                    result, status = fit_scores(proxy[ix][:, take], outside[ix][:, take],
                                               obs[['full_RNA_total','full_detected_genes']].iloc[ix].to_numpy(float), blocks)
                idx = {j:k for k, j in enumerate(np.flatnonzero(take))}
                for j, pid in enumerate(ids):
                    old = lookup.loc[(info['gse_id'],info['major'],pid,axis)]
                    row = dict(common, programme_id=pid, partial_signature_genes=int(old.partial_signature_genes),
                               status=status if take[j] else 'PROXY_OR_OUTSIDE_REFERENCE_UNMEASURABLE',
                               raw_block_weighted_rho=np.nan, adjusted_cell_rho=np.nan, adjusted_block_rho=np.nan,
                               proxy_RNA_rho=np.nan, proxy_detected_rho=np.nan)
                    if result is not None and take[j]:
                        k = idx[j]
                        for field, src in [('raw_block_weighted_rho','raw'),('adjusted_cell_rho','adjusted'),('adjusted_block_rho','block'),('proxy_RNA_rho','RNA'),('proxy_detected_rho','detected')]:
                            row[field] = float(result[src][k])
                        if not np.isfinite(row['adjusted_cell_rho']) or not np.isfinite(row['adjusted_block_rho']):
                            row['status'] = 'CONSTANT_OR_QUALITY_EXPLAINED_SCORE'
                        checked += 1
                    rows.append(row)
                if result is not None: maxerr = max(maxerr, result['max_schur_error'])
            processed += 1
        (out/'progress.json').write_text(json.dumps(dict(stage='EXISTING_METADATA_AND_SCORES_NO_NEW_RNA_OR_SPACE_EFFECT',native_files=ordinal+1,score_files=processed,metadata_cells=metadata_cells,descriptor_rows=len(rows),observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert metadata_cells == 91447 and processed == 2836
    table = pd.DataFrame(rows); assert len(table) == 725004
    table.to_csv(out/'private_all_identity_quality_adjusted_fidelity.tsv.gz', sep='\t', index=False, compression='gzip')
    metrics = ['raw_block_weighted_rho','adjusted_cell_rho','adjusted_block_rho','proxy_RNA_rho','proxy_detected_rho']
    source = []
    keys = ['gse_id','major','programme_id','axis_id']
    for key, g in table.groupby(keys, sort=True):
        use = g[g.status.eq('ESTIMABLE_TECHNICAL_DESCRIPTION')]
        row = dict(zip(keys, key)); row.update(identity_strata=len(g),estimable_identity_strata=len(use),
                                              partial_signature_genes=int(g.partial_signature_genes.iloc[0]))
        for field in metrics: row[field] = float(use[field].median()) if len(use) else np.nan
        source.append(row)
    source = pd.DataFrame(source); assert len(source) == 127636
    source.to_csv(out/'private_all_source_major_programme_axis_fidelity.tsv.gz',sep='\t',index=False,compression='gzip')
    triage = []
    for pid, r in sorted(recs.items()):
        for major in r['allowed_majors']:
            for axis in axes:
                all_g = source[source.programme_id.eq(pid)&source.major.eq(major)&source.axis_id.eq(axis)]
                g = all_g[all_g.estimable_identity_strata.gt(0)&all_g.partial_signature_genes.ge(5)].copy()
                assert g.gse_id.is_unique
                cell = g.adjusted_cell_rho.to_numpy(float); block = g.adjusted_block_rho.to_numpy(float)
                ns = len(g); cm = float(np.median(cell)) if ns else np.nan; bm = float(np.median(block)) if ns else np.nan
                concord = float(((cell >= .2)&(block >= .2)).mean()) if ns else np.nan
                loo = [min(np.median(np.delete(cell,k)), np.median(np.delete(block,k))) for k in range(ns)] if ns > 1 else []
                loomin = float(min(loo)) if loo else np.nan
                stable = bool(r['algorithm_stable_candidate'])
                passed = bool(stable and ns >= 3 and cm >= .3 and bm >= .3 and concord >= .75 and loomin >= .3)
                triage.append(dict(programme_id=pid,major=major,axis_id=axis,available_sources=len(all_g),qualifiable_sources=ns,
                                   algorithm_stable_candidate=stable,source_median_adjusted_cell_rho=cm,
                                   source_median_adjusted_block_rho=bm,source_fraction_both_ge_point2=concord,
                                   leave_one_source_min_of_two_medians=loomin,exploratory_proxy_triage_pass=passed,
                                   full_signature_spatial_gate_changed=False,formal_state_or_mechanism_validation=False))
    triage = pd.DataFrame(triage)
    triage.to_csv(out/'private_complete_programme_major_axis_triage.tsv',sep='\t',index=False)
    coverage = []
    for (major,axis),g in triage.groupby(['major','axis_id'],sort=True):
        coverage.append(dict(major=major,axis_id=axis,programme_positions=len(g),positions_ge3_qualifiable_sources=int(g.qualifiable_sources.ge(3).sum()),
                             provisional_pass_positions=int(g.exploratory_proxy_triage_pass.sum()),no_new_spatial_effect_or_clinical_inference=True))
    pd.DataFrame(coverage).to_csv(out/'private_major_axis_coverage_counts.tsv',sep='\t',index=False)
    terminal = dict(status='COMPLETE_ALL_EXISTING_CELLS_SCORES_QUALITY_ADJUSTED_DISJOINT_REFERENCE_TECHNICAL_TRIAGE_PENDING_AUDIT',
                    completed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),
                    metadata_cells=metadata_cells,score_files=processed,identity_descriptions=len(table),source_descriptions=len(source),
                    triage_positions=len(triage),schur_positions_checked=checked,max_schur_error=maxerr,
                    no_new_expression_or_spatial_effects_or_outcomes_read=True,original_full_signature_gates_changed=False,
                    no_formal_state_or_clinical_or_interaction_validation=True,scientific_outputs_public=False,
                    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(d['bytes'] for d in terminal['files'].values()) < 64*1024**2
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps({k:terminal[k] for k in ['status','metadata_cells','score_files','identity_descriptions','source_descriptions','triage_positions','schur_positions_checked','max_schur_error']}))


if __name__ == '__main__':
    main()
