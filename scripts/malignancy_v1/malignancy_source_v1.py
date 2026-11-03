# Uniform malignancy calling for one source of integration v2 (all 72 sources rerun with the same
# pipeline, 2026-11-02). Two passes over the units (samples) of the source:
#  pass 1, per unit: quick clustering of raw counts (normalize 1e4, log1p, 2000 HVG, 30 PCs, 15-NN,
#    Leiden 1.0) and cluster lineage by marker means: immune (PTPRC), stromal = max of fibroblast
#    (COL1A1/COL1A2/DCN), endothelial (PECAM1/VWF) and pericyte (RGS5) scores, else candidate. For DLBCL
#    sources, B-cell clusters (MS4A1/CD79A/CD19) are candidates, not reference. scATOMIC on all cells of
#    the unit; immune/stromal clusters that scATOMIC calls >=50% tumour become candidates
#    ("candidate_from_<lineage>") so they are scored rather than used as CNV baseline. Cached per unit
#    in WORK so an interrupted source resumes without rerunning scATOMIC.
#  pass 2, per unit: inferCNV. Reference pool = the unit's remaining immune/stromal cells; when the unit
#    has fewer than MIN_OWN (200) of them, the pool is padded to POOL_TARGET (400) with immune/stromal
#    cells drawn from the other units of the same source (same technology and processing; flagged
#    "ref_borrowed" in the unit record). Up to 1/4 of the pool (max 500) is held out and scored as
#    "control"; the rest (max 1500) is the reference; observations = the unit's candidate cells.
#  CNV metrics per cell: score = mean squared deviation of the denoised inferCNV profile from 1;
#    corr = correlation with the mean profile of the top 5% scoring observation cells (Tirosh/Neftel).
#  CNV call: thresholds = 99th percentiles of the control cells' score and corr (needs >= 20 controls);
#    a candidate cluster is positive when its median score or median corr exceeds its threshold
#    (cluster-level call used).
#  consensus: immune/stromal lineage -> non_malignant; candidate with scATOMIC tumour and CNV positive ->
#    malignant; scATOMIC normal and CNV negative -> non_malignant; disagreement or CNV not evaluable ->
#    uncertain (kept in the TME and flagged, per user decision).
# Output: per-source TSV of per-cell evidence + per-unit status JSON lines.
# Usage: malignancy_source_v1.py SHARD.h5ad WORK_DIR OUT.tsv RSCRIPT WORKER GENE_POS NCORES SEED
import sys, os, json, subprocess, shutil
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io, anndata as ad, scanpy as sc

shard, work, out, rscript, worker, gene_pos, ncores, seed = sys.argv[1:9]
ncores, seed = int(ncores), int(seed)
IMM, BCELL = ['PTPRC'], ['MS4A1', 'CD79A', 'CD19']
MIN_OWN, POOL_TARGET = 200, 400
rng = np.random.default_rng(seed)
a = ad.read_h5ad(shard)
counts = a.layers['counts'].tocsr() if 'counts' in a.layers else a.X.tocsr()
main = str(a.obs['analysis_main'].iloc[0])
dlbcl = main == 'Diffuse large B-cell lymphoma'
groups = list(a.obs.groupby('unit_id', observed=True).indices.items())


def run_r(d, mode, cores=None, timeout=None):
    # own process group so a timeout also kills forked workers (scATOMIC/MAGIC hung >20 h on 4 units)
    p = subprocess.Popen([rscript, worker, d, gene_pos, str(cores or ncores), str(seed), mode], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True, start_new_session=True)
    try:
        o, e = p.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, 9); o, e = p.communicate()
        e = f'timeout after {timeout} s\n' + (e or '')
    return subprocess.CompletedProcess(p.args, p.returncode, o, e)


def run_scatomic(d):
    # first try with all cores (2 h), then single-core (4 h); a unit that still fails keeps scATOMIC "failed"
    r = run_r(d, 'scatomic', timeout=7200)
    if not os.path.exists(os.path.join(d, 'scatomic.tsv')):
        r = run_r(d, 'scatomic', cores=1, timeout=14400)
    return r


def write_unit(d, idx, cid, lin, grp):
    sub = counts[idx]
    keep_g = np.asarray((sub > 0).sum(0)).ravel() > 0
    scipy.io.mmwrite(os.path.join(d, 'counts.mtx'), sub[:, keep_g].T.tocoo())
    open(os.path.join(d, 'genes.txt'), 'w').write('\n'.join(a.var_names[keep_g].astype(str)) + '\n')
    pd.DataFrame({'cell': cid, 'lineage': lin, 'infercnv_group': grp}).to_csv(os.path.join(d, 'cells.tsv'), sep='\t', index=False)


# pass 1
p1 = {}
for k, (u, idx) in enumerate(groups):
    idx = np.asarray(idx); cid = a.obs_names[idx].astype(str)
    d = os.path.join(work, f'u{k:04d}'); os.makedirs(d, exist_ok=True)
    cache = os.path.join(d, 'pass1.tsv.gz')
    if os.path.exists(cache):
        p1[u] = pd.read_csv(cache, sep='\t', dtype={'quick_cluster': str}).set_index('cell'); continue
    sub = counts[idx]
    keep_g = np.asarray((sub > 0).sum(0)).ravel() > 0
    b = ad.AnnData(X=sub[:, keep_g].astype(np.float32), var=pd.DataFrame(index=a.var_names[keep_g].astype(str)))
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    sc.pp.highly_variable_genes(b, n_top_genes=min(2000, b.n_vars), flavor='seurat')
    sc.pp.pca(b, n_comps=min(30, b.n_obs - 1), mask_var='highly_variable', random_state=seed)
    sc.pp.neighbors(b, n_neighbors=15, random_state=seed)
    sc.tl.leiden(b, resolution=1.0, random_state=seed, flavor='igraph', n_iterations=2, directed=False)
    cl = b.obs['leiden'].astype(str).values

    def mscore(gs):
        gs = [g for g in gs if g in b.var_names]
        return pd.Series(np.asarray(b[:, gs].X.mean(1)).ravel() if gs else np.zeros(b.n_obs), index=range(b.n_obs)).groupby(cl).mean()
    sim, sb = mscore(IMM), mscore(BCELL)
    sst = pd.concat([mscore(['COL1A1', 'COL1A2', 'DCN']), mscore(['PECAM1', 'VWF']), mscore(['RGS5'])], axis=1).max(1)  # fibroblast / endothelial / pericyte
    lin_cl = {}
    for c in sim.index:
        if dlbcl and sb[c] > 1.0:
            lin_cl[c] = 'candidate_B'
        elif sim[c] >= 1.0 and sim[c] >= sst[c]:
            lin_cl[c] = 'immune'
        elif sst[c] >= 1.0:
            lin_cl[c] = 'stromal'
        else:
            lin_cl[c] = 'candidate'
    lin = np.array([lin_cl[c] for c in cl])
    write_unit(d, idx, cid, lin, 'none')
    r1 = run_scatomic(d)
    df = pd.DataFrame({'cell': cid, 'quick_cluster': cl, 'lineage': lin})
    tumour = np.zeros(len(idx), bool)
    if os.path.exists(os.path.join(d, 'scatomic.tsv')):
        scat = pd.read_csv(os.path.join(d, 'scatomic.tsv'), sep='\t').set_index('cell').reindex(cid)
        pred = scat['scATOMIC_pred'].astype(str)
        pc_ = scat['pan_cancer_cluster'].astype(str) if 'pan_cancer_cluster' in scat else pd.Series('', index=scat.index)
        tumour = (pc_.eq('Cancer') | pred.str.contains('Cancer Cell|Tumor|Tumour|Neoplastic', case=False)).values
        df['scATOMIC_pred'] = pred.values; df['scATOMIC_pan_cancer_cluster'] = pc_.values
        df['scATOMIC_tumour'] = np.where(tumour, 'tumour', np.where(pred.isin(['nan', '', 'unclassified_any_cell', 'unclassified']).values, 'unclassified', 'normal'))
    else:
        df['scATOMIC_pred'] = ''; df['scATOMIC_pan_cancer_cluster'] = ''; df['scATOMIC_tumour'] = 'not_run'
    # immune/stromal clusters that scATOMIC calls mostly tumour are not trusted as CNV reference
    tum_frac = pd.Series(tumour).groupby(cl).mean()
    df['lineage'] = [('candidate_from_' + l) if (l in ('immune', 'stromal') and tum_frac[c] >= 0.5) else l for l, c in zip(lin, cl)]
    if 'scatomic' not in (json.load(open(os.path.join(d, 'status.json'))) if os.path.exists(os.path.join(d, 'status.json')) else {}):
        json.dump({'scatomic': 'failed: worker error ' + r1.stderr[-300:]}, open(os.path.join(d, 'status.json'), 'w'))
    df.to_csv(cache, sep='\t', index=False, compression='gzip')
    p1[u] = df.set_index('cell')
    print(json.dumps({'pass1': u, 'cells': len(idx)}), flush=True)

# pass 2
pos = pd.Series(np.arange(a.n_obs), index=a.obs_names.astype(str))
allref = pd.concat([pd.DataFrame({'cell': p.index[p['lineage'].isin(['immune', 'stromal'])], 'unit': u}) for u, p in p1.items()], ignore_index=True)
rows, units = [], []
for k, (u, idx) in enumerate(groups):
    idx = np.asarray(idx); cid = a.obs_names[idx].astype(str)
    d = os.path.join(work, f'u{k:04d}')
    p = p1[u].reindex(cid)
    lin = p['lineage'].astype(str).values
    own = np.flatnonzero(np.isin(lin, ['immune', 'stromal']))
    rng.shuffle(own)
    borrowed = np.array([], dtype=str)
    if len(own) < MIN_OWN:
        other = allref.loc[allref.unit != u, 'cell'].values
        borrowed = rng.choice(other, min(POOL_TARGET - len(own), len(other)), replace=False) if len(other) else borrowed
    pool_cells = np.concatenate([cid[own], borrowed])
    nctrl = min(500, len(pool_cells) // 4)
    ctrl_cells, ref_cells = pool_cells[:nctrl], pool_cells[nctrl:][:1500]
    obs_mask = np.char.startswith(lin.astype(str), 'candidate')
    grp = pd.Series('none', index=cid, dtype=object)
    grp[obs_mask] = 'observation'
    grp[[c for c in ref_cells if c in grp.index]] = 'reference'; grp[[c for c in ctrl_cells if c in grp.index]] = 'control'
    bor = [c for c in borrowed]
    bgrp = pd.Series('reference', index=bor, dtype=object)
    bgrp[[c for c in bor if c in set(ctrl_cells)]] = 'control'
    bgrp = bgrp[[c in set(ref_cells) or c in set(ctrl_cells) for c in bor]]
    all_cells = np.concatenate([cid, bgrp.index.values])
    all_grp = np.concatenate([grp.values, bgrp.values])
    all_lin = np.concatenate([lin, np.array(['borrowed_reference'] * len(bgrp))])
    write_unit(d, pos[all_cells].values, all_cells, all_lin, all_grp)
    st_path = os.path.join(d, 'status.json')
    if os.path.exists(os.path.join(d, 'cnv_scores.tsv')):
        os.remove(os.path.join(d, 'cnv_scores.tsv'))
    st0 = json.load(open(st_path)) if os.path.exists(st_path) else {}
    st0.pop('infercnv', None); json.dump(st0, open(st_path, 'w'))
    r2 = run_r(d, 'infercnv', timeout=21600)
    st = json.load(open(st_path)) if os.path.exists(st_path) else {}
    if 'infercnv' not in st:
        st['infercnv'] = 'failed: worker error ' + r2.stderr[-300:]
    df = pd.DataFrame({'cell_id': cid, 'unit_id': u, 'quick_cluster': p['quick_cluster'].values, 'lineage': lin,
                       'infercnv_group': grp.values, 'scATOMIC_pred': p['scATOMIC_pred'].values,
                       'scATOMIC_pan_cancer_cluster': p['scATOMIC_pan_cancer_cluster'].values, 'scATOMIC_tumour': p['scATOMIC_tumour'].values})
    thr = np.nan
    if os.path.exists(os.path.join(d, 'cnv_scores.tsv')):
        cs = pd.read_csv(os.path.join(d, 'cnv_scores.tsv'), sep='\t').set_index('cell')
        df['cnv_score'] = cs['cnv_score'].reindex(cid).values; df['cnv_corr'] = cs['cnv_corr'].reindex(cid).values
        ctrl_scores = cs.reindex(ctrl_cells).dropna()
        if len(ctrl_scores) >= 20:
            thr = float(np.quantile(ctrl_scores['cnv_score'], 0.99)); thr_c = float(np.quantile(ctrl_scores['cnv_corr'], 0.99))
            df['cnv_corr_threshold'] = thr_c
            df['cnv_cell_positive'] = (df['cnv_score'] > thr) | (df['cnv_corr'] > thr_c)
            obs = df.infercnv_group == 'observation'
            clmed = df[obs].groupby('quick_cluster')[['cnv_score', 'cnv_corr']].median()
            pos_cl = (clmed['cnv_score'] > thr) | (clmed['cnv_corr'] > thr_c)
            df['cnv_call'] = np.where(obs, np.where(df['quick_cluster'].map(pos_cl).fillna(False).astype(bool), 'cnv_positive', 'cnv_negative'), 'reference_lineage')
        else:
            df['cnv_cell_positive'] = np.nan; df['cnv_call'] = np.where(df.infercnv_group == 'observation', 'not_evaluable', 'reference_lineage')
    else:
        df['cnv_score'] = np.nan; df['cnv_corr'] = np.nan; df['cnv_cell_positive'] = np.nan
        df['cnv_call'] = np.where(df.infercnv_group == 'observation', 'not_evaluable', 'reference_lineage')
    df['cnv_threshold'] = thr
    cand = df['lineage'].str.startswith('candidate')
    tum, cnv = df['scATOMIC_tumour'], df['cnv_call']
    df['malignancy_call'] = np.select(
        [~cand, cand & (tum == 'tumour') & (cnv == 'cnv_positive'), cand & (tum.isin(['normal'])) & (cnv == 'cnv_negative')],
        ['non_malignant', 'malignant', 'non_malignant'], 'uncertain')
    df['malignancy_rule'] = np.select(
        [~cand, df['malignancy_call'] == 'malignant', df['malignancy_call'] == 'non_malignant', cnv == 'not_evaluable'],
        ['immune_or_stromal_lineage', 'scATOMIC_tumour_and_CNV_positive', 'scATOMIC_normal_and_CNV_negative', 'CNV_not_evaluable'], 'methods_disagree')
    rows.append(df)
    rec = {'unit': u, 'cells': int(len(idx)), 'candidates': int(cand.sum()), 'own_reference_pool': int(len(own)), 'ref_borrowed': int(len(bgrp)),
           'reference': int(len(ref_cells)), 'control': int(len(ctrl_cells)), 'cnv_threshold': thr,
           **{kk: v for kk, v in st.items() if not kk.endswith('_seconds')},
           'seconds': round(sum(v for kk, v in st.items() if kk.endswith('_seconds')), 1),
           'calls': df['malignancy_call'].value_counts().to_dict()}
    units.append(rec); print(json.dumps(rec, default=str), flush=True)
pd.concat(rows).to_csv(out, sep='\t', index=False, compression='gzip')
with open(out + '.units.jsonl', 'w') as f:
    for r in units:
        f.write(json.dumps(r, default=str) + '\n')
shutil.rmtree(work, ignore_errors=True)
