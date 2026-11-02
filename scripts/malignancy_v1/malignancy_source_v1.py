# Uniform malignancy calling for one source of integration v2 (all 72 sources rerun with the same
# pipeline, 2026-11-02). Per unit (sample):
#  1 quick clustering of raw counts (normalize 1e4, log1p, 2000 HVG, 30 PCs, 15-NN, Leiden 1.0) and
#    cluster lineage by marker means: immune (PTPRC), stromal = max of fibroblast (COL1A1/COL1A2/DCN),
#    endothelial (PECAM1/VWF) and pericyte (RGS5) scores,
#    else candidate. For DLBCL sources, B-cell clusters (MS4A1/CD79A/CD19) are candidates, not reference.
#  2 R worker: scATOMIC on all cells first; immune/stromal clusters that scATOMIC calls >=50% tumour
#    are moved to candidates; then inferCNV with reference = up to 1500 remaining immune/stromal cells,
#    a held-out up to 500 of them scored as "control", observations = candidate cells.
#  3 CNV call: positive if a cell's score > 99th percentile of the control scores; candidate clusters
#    are positive when their median score exceeds that threshold (cluster-level call used).
#  4 consensus: immune/stromal lineage -> non_malignant; candidate with scATOMIC tumour and CNV
#    positive -> malignant; both negative -> non_malignant; disagreement or CNV not evaluable ->
#    uncertain (kept in the TME and flagged, per user decision).
# Output: per-source TSV of per-cell evidence + per-unit status JSON lines.
# Usage: malignancy_source_v1.py SHARD.h5ad WORK_DIR OUT.tsv RSCRIPT WORKER GENE_POS NCORES SEED
import sys, os, json, subprocess, shutil
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io, anndata as ad, scanpy as sc

shard, work, out, rscript, worker, gene_pos, ncores, seed = sys.argv[1:9]
ncores, seed = int(ncores), int(seed)
IMM, BCELL = ['PTPRC'], ['MS4A1', 'CD79A', 'CD19']
rng = np.random.default_rng(seed)
a = ad.read_h5ad(shard)
counts = a.layers['counts'].tocsr() if 'counts' in a.layers else a.X.tocsr()
main = str(a.obs['analysis_main'].iloc[0])
dlbcl = main == 'Diffuse large B-cell lymphoma'
rows, units = [], []
for u, idx in a.obs.groupby('unit_id', observed=True).indices.items():
    idx = np.asarray(idx)
    sub = counts[idx]
    keep_g = np.asarray((sub > 0).sum(0)).ravel() > 0
    sub = sub[:, keep_g]; genes = a.var_names[keep_g].astype(str)
    b = ad.AnnData(X=sub.astype(np.float32), var=pd.DataFrame(index=genes))
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
    d = os.path.join(work, f'u{len(units):04d}'); os.makedirs(d, exist_ok=True)
    scipy.io.mmwrite(os.path.join(d, 'counts.mtx'), sub.T.tocoo())
    open(os.path.join(d, 'genes.txt'), 'w').write('\n'.join(genes) + '\n')
    cid = a.obs_names[idx].astype(str)
    run_r = lambda mode: subprocess.run([rscript, worker, d, gene_pos, str(ncores), str(seed), mode], capture_output=True, text=True)
    pd.DataFrame({'cell': cid, 'lineage': lin, 'infercnv_group': 'none'}).to_csv(os.path.join(d, 'cells.tsv'), sep='\t', index=False)
    r1 = run_r('scatomic')
    tumour = np.zeros(len(idx), bool); scat = None
    if os.path.exists(os.path.join(d, 'scatomic.tsv')):
        scat = pd.read_csv(os.path.join(d, 'scatomic.tsv'), sep='\t').set_index('cell').reindex(cid)
        p_ = scat['scATOMIC_pred'].astype(str)
        pc_ = scat['pan_cancer_cluster'].astype(str) if 'pan_cancer_cluster' in scat else pd.Series('', index=scat.index)
        tumour = (pc_.eq('Cancer') | p_.str.contains('Cancer Cell|Tumor|Tumour|Neoplastic', case=False)).values
    # immune/stromal clusters that scATOMIC calls mostly tumour are not trusted as CNV reference: they
    # become candidates (observations) so their CNV is scored rather than used as baseline
    tum_frac = pd.Series(tumour).groupby(cl).mean()
    lin = np.array([('candidate_from_' + l) if (l in ('immune', 'stromal') and tum_frac[c] >= 0.5) else l for l, c in zip(lin, cl)])
    refpool = np.flatnonzero(np.isin(lin, ['immune', 'stromal']))
    rng.shuffle(refpool)
    ctrl, ref = refpool[:min(500, len(refpool) // 4)], refpool[min(500, len(refpool) // 4):][:1500]
    grp = np.array(['none'] * len(idx), dtype=object)
    grp[np.char.startswith(lin.astype(str), 'candidate')] = 'observation'
    grp[ref] = 'reference'; grp[ctrl] = 'control'
    pd.DataFrame({'cell': cid, 'lineage': lin, 'infercnv_group': grp}).to_csv(os.path.join(d, 'cells.tsv'), sep='\t', index=False)
    r2 = run_r('infercnv')
    st = json.load(open(os.path.join(d, 'status.json'))) if os.path.exists(os.path.join(d, 'status.json')) else {'worker_error': (r1.stderr + r2.stderr)[-300:]}
    df = pd.DataFrame({'cell_id': cid, 'unit_id': u, 'quick_cluster': cl, 'lineage': lin, 'infercnv_group': grp})
    # scATOMIC
    if scat is not None:
        pred = scat['scATOMIC_pred'].astype(str).values
        df['scATOMIC_pred'] = pred
        df['scATOMIC_pan_cancer_cluster'] = scat['pan_cancer_cluster'].astype(str).values if 'pan_cancer_cluster' in scat else ''
        df['scATOMIC_tumour'] = np.where(tumour, 'tumour', np.where(pd.Series(pred).isin(['nan', '', 'unclassified_any_cell', 'unclassified']), 'unclassified', 'normal'))
    else:
        df['scATOMIC_pred'] = ''; df['scATOMIC_pan_cancer_cluster'] = ''; df['scATOMIC_tumour'] = 'not_run'
    # inferCNV
    thr = np.nan
    if os.path.exists(os.path.join(d, 'cnv_scores.tsv')):
        cs = pd.read_csv(os.path.join(d, 'cnv_scores.tsv'), sep='\t').set_index('cell')['cnv_score']
        df['cnv_score'] = cs.reindex(cid).values
        ctrl_scores = df.loc[df.infercnv_group == 'control', 'cnv_score'].dropna()
        if len(ctrl_scores) >= 20:
            thr = float(np.quantile(ctrl_scores, 0.99))
            df['cnv_cell_positive'] = df['cnv_score'] > thr
            obs = df.infercnv_group == 'observation'
            clmed = df[obs].groupby('quick_cluster')['cnv_score'].median()
            df['cnv_call'] = np.where(obs, np.where(df['quick_cluster'].map(clmed) > thr, 'cnv_positive', 'cnv_negative'), 'reference_lineage')
        else:
            df['cnv_cell_positive'] = np.nan; df['cnv_call'] = np.where(df.infercnv_group == 'observation', 'not_evaluable', 'reference_lineage')
    else:
        df['cnv_score'] = np.nan; df['cnv_cell_positive'] = np.nan
        df['cnv_call'] = np.where(df.infercnv_group == 'observation', 'not_evaluable', 'reference_lineage')
    df['cnv_threshold'] = thr
    # consensus
    cand = df['lineage'].str.startswith('candidate')
    tum, cnv = df['scATOMIC_tumour'], df['cnv_call']
    df['malignancy_call'] = np.select(
        [~cand, cand & (tum == 'tumour') & (cnv == 'cnv_positive'), cand & (tum.isin(['normal'])) & (cnv == 'cnv_negative')],
        ['non_malignant', 'malignant', 'non_malignant'], 'uncertain')
    df['malignancy_rule'] = np.select(
        [~cand, df['malignancy_call'] == 'malignant', df['malignancy_call'] == 'non_malignant', cnv == 'not_evaluable'],
        ['immune_or_stromal_lineage', 'scATOMIC_tumour_and_CNV_positive', 'scATOMIC_normal_and_CNV_negative', 'CNV_not_evaluable'], 'methods_disagree')
    rows.append(df)
    rec = {'unit': u, 'cells': int(len(idx)), 'candidates': int(cand.sum()), 'reference': int(len(ref)), 'control': int(len(ctrl)),
           'cnv_threshold': thr, **{k: v for k, v in st.items() if not k.endswith('_seconds')},
           'seconds': round(sum(v for k, v in st.items() if k.endswith('_seconds')), 1),
           'calls': df['malignancy_call'].value_counts().to_dict()}
    units.append(rec); print(json.dumps(rec, default=str), flush=True)
    shutil.rmtree(d, ignore_errors=True)
pd.concat(rows).to_csv(out, sep='\t', index=False, compression='gzip')
with open(out + '.units.jsonl', 'w') as f:
    for r in units:
        f.write(json.dumps(r, default=str) + '\n')
