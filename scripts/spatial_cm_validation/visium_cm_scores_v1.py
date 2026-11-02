#!/usr/bin/env python3
"""Per-section CM activity, colocalization and aggregation scores on cell2location output.

Usage: visium_cm_scores_v1.py --c2l DIR --fit FIT_RUN_DIR --ct75 NODE_CATALOGUE.tsv --output DIR
"""
import argparse, json, sys, zlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, binomtest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cm_spatial_scores_v1 import (load_cm, cm_activity, scale_activity, spearman_matrix, knn_weights,
                                  bivariate_moran_matrix, coloc_score, random_member_null, edge_values)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--c2l', required=True)
    ap.add_argument('--fit', required=True)
    ap.add_argument('--ct75', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--abundance', default='q05', choices=['q05', 'mean'])
    ap.add_argument('--n-null', type=int, default=1000)
    a = ap.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    W, members, edges = load_cm(a.fit)
    cat = pd.read_csv(a.ct75, sep='\t'); pool = cat.node.tolist(); strata = cat.set_index('node').parent.to_dict()
    secs = sorted((Path(a.c2l) / 'sections').glob('GSM*')) or [Path(a.c2l)]
    spots = pd.concat([pd.read_csv(d / 'visium_spots_v1.tsv.gz', sep='\t', index_col=0) for d in secs])
    ab = pd.concat([pd.read_csv(d / f'visium_{a.abundance}_abundance_v1.tsv.gz', sep='\t', index_col=0) for d in secs])
    assert ab.index.equals(spots.index) and ab.notna().all().all()
    names = ab.columns.tolist()
    act = scale_activity(cm_activity(ab, W, members))
    act_all = scale_activity(cm_activity(ab, W, {cm: [n for n in W.index] for cm in W.columns}))
    act.join(act_all, rsuffix='_all75').join(spots[['sample']]).to_csv(out / 'visium_CM_activity_v1.tsv.gz', sep='\t')
    rows, erows = [], []
    for s, idx in spots.groupby('sample').groups.items():
        X = ab.loc[idx].to_numpy(); coords = spots.loc[idx, ['pxl_col', 'pxl_row']].to_numpy()
        C = spearman_matrix(X); M = bivariate_moran_matrix(X, knn_weights(coords, 6))
        for cm, mem in members.items():
            for stat, mat in [('colocalization_Spearman', C), ('aggregation_bivariate_Moran', M)]:
                score, r = coloc_score(mat, names, mem)
                sd = zlib.crc32(f"{s}|{cm}|{stat}".encode())
                null = random_member_null(mat, names, mem, pool, a.n_null, seed=sd)
                lnull = random_member_null(mat, names, mem, pool, a.n_null, seed=sd + 1, strata=strata)
                rows.append(dict(sample=s, n_spots=len(idx), cm=cm, statistic=stat, n_members=len(mem), score=score,
                                 within_CM_median=r, null_median=float(np.median(null)) if len(null) else np.nan,
                                 null_percentile=float((null < score).mean()) if len(null) else np.nan,
                                 lineage_null_median=float(np.median(lnull)) if len(lnull) else np.nan,
                                 lineage_null_percentile=float((lnull < score).mean()) if len(lnull) else np.nan,
                                 mean_activity=float(act.loc[idx, cm].mean())))
            for e, v in edge_values(C, names, edges[cm]).items():
                erows.append(dict(sample=s, cm=cm, edge=e, spearman=v, moran=edge_values(M, names, [tuple(e.split('~'))])[e]))
    res = pd.DataFrame(rows); res.to_csv(out / 'visium_CM_section_scores_v1.tsv', sep='\t', index=False)
    er = pd.DataFrame(erows); er.to_csv(out / 'visium_CM_edge_section_values_v1.tsv', sep='\t', index=False)
    summ = []
    for (cm, stat), g in res.groupby(['cm', 'statistic']):
        sc_ = g.score.dropna(); pc = g.null_percentile.dropna(); lp = g.lineage_null_percentile.dropna()
        summ.append(dict(cm=cm, statistic=stat, n_sections=len(sc_), median_score=float(sc_.median()),
                         n_sections_score_gt_0_5=int((sc_ > .5).sum()),
                         wilcoxon_vs_0_5_p=float(wilcoxon(sc_ - .5).pvalue) if len(sc_) >= 5 else np.nan,
                         median_null_percentile=float(pc.median()),
                         wilcoxon_null_percentile_vs_0_5_p=float(wilcoxon(pc - .5).pvalue) if len(pc) >= 5 else np.nan,
                         median_lineage_null_percentile=float(lp.median()),
                         wilcoxon_lineage_null_percentile_vs_0_5_p=float(wilcoxon(lp - .5).pvalue) if len(lp) >= 5 else np.nan))
    esum = []
    for (cm, e), g in er.groupby(['cm', 'edge']):
        v = g.spearman.dropna()
        esum.append(dict(cm=cm, edge=e, n_sections=len(v), median_spearman=float(v.median()), n_positive=int((v > 0).sum()),
                         sign_test_p=float(binomtest(int((v > 0).sum()), len(v)).pvalue) if len(v) else np.nan,
                         median_moran=float(g.moran.median())))
    pd.DataFrame(summ).to_csv(out / 'visium_CM_score_summary_v1.tsv', sep='\t', index=False)
    pd.DataFrame(esum).to_csv(out / 'visium_CM_edge_summary_v1.tsv', sep='\t', index=False)
    print(pd.DataFrame(summ).to_string()); print(pd.DataFrame(esum).to_string())


if __name__ == '__main__':
    main()
