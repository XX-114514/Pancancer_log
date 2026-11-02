#!/usr/bin/env python3
"""Post hoc control: deconvolution-induced correlation between subsets with similar reference signatures.

cell2location tends to split a shared signal between subsets whose reference
signatures are alike, so their spot abundances correlate whatever the tissue
organisation. Per section, each pair's Spearman correlation is replaced by its
residual from the median correlation of all pairs within +/- --bin of the same
signature similarity (Pearson of log1p CP10K signatures). The CoVarNet
colocalization score and both nulls are recomputed on the residual matrix, and
each CM edge gets its residual and the percentile among signature-matched pairs.

Usage: visium_signature_control_v1.py --c2l DIR --fit FIT_RUN --ct75 CATALOGUE --output DIR
"""
import argparse, sys, zlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, binomtest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cm_spatial_scores_v1 import load_cm, spearman_matrix, coloc_score, random_member_null


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--c2l', required=True); ap.add_argument('--fit', required=True)
    ap.add_argument('--ct75', required=True); ap.add_argument('--output', required=True)
    ap.add_argument('--bin', type=float, default=0.01); ap.add_argument('--n-null', type=int, default=1000)
    a = ap.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    W, members, edges = load_cm(a.fit)
    cat = pd.read_csv(a.ct75, sep='\t'); pool = cat.node.tolist(); strata = cat.set_index('node').parent.to_dict()
    S = pd.read_csv(Path(a.c2l) / 'reference_signatures_v1.tsv', sep='\t', index_col=0)
    names = list(S.columns); L = np.log1p(S / S.sum() * 1e4); sig = np.corrcoef(L.T.to_numpy())
    iu = np.triu_indices(len(names), 1); sv = sig[iu]
    order = np.argsort(sv); svs = sv[order]
    lo = np.searchsorted(svs, sv - a.bin); hi = np.searchsorted(svs, sv + a.bin, side='right')
    rows, erows = [], []
    for d in sorted((Path(a.c2l) / 'sections').glob('GSM*')):
        ab = pd.read_csv(d / 'visium_q05_abundance_v1.tsv.gz', sep='\t', index_col=0)[names]
        s = pd.read_csv(d / 'visium_spots_v1.tsv.gz', sep='\t', index_col=0)['sample'].iloc[0]
        C = spearman_matrix(ab.to_numpy()); cv = C[iu]; cvs = cv[order]
        expect = np.array([np.median(cvs[l:h]) for l, h in zip(lo, hi)])
        pct = np.array([(cvs[l:h] < c).mean() for l, h, c in zip(lo, hi, cv)])
        Rm = np.zeros_like(C); Rm[iu] = cv - expect; Rm = Rm + Rm.T
        for cm, mem in members.items():
            score, r = coloc_score(Rm, names, mem)
            sd = zlib.crc32(f'{s}|{cm}|resid'.encode())
            null = random_member_null(Rm, names, mem, pool, a.n_null, seed=sd)
            lnull = random_member_null(Rm, names, mem, pool, a.n_null, seed=sd + 1, strata=strata)
            rows.append(dict(sample=s, cm=cm, statistic='signature_adjusted_colocalization', score=score, within_CM_median_residual=r,
                             null_median=float(np.median(null)), null_percentile=float((null < score).mean()),
                             lineage_null_median=float(np.median(lnull)), lineage_null_percentile=float((lnull < score).mean())))
            for e in edges[cm]:
                i, j = sorted((names.index(e[0]), names.index(e[1])))
                k = np.flatnonzero((iu[0] == i) & (iu[1] == j))[0]
                erows.append(dict(sample=s, cm=cm, edge=f'{e[0]}~{e[1]}', signature_r=sv[k], spearman=cv[k],
                                  signature_matched_expected=expect[k], residual=cv[k] - expect[k], percentile_among_signature_matched=pct[k]))
    res = pd.DataFrame(rows); res.to_csv(out / 'visium_signature_adjusted_section_scores_v1.tsv', sep='\t', index=False)
    er = pd.DataFrame(erows); er.to_csv(out / 'visium_signature_adjusted_edge_values_v1.tsv', sep='\t', index=False)
    summ = []
    for cm, g in res.groupby('cm'):
        summ.append(dict(cm=cm, n_sections=len(g), median_score=g.score.median(), n_sections_score_gt_0_5=int((g.score > .5).sum()),
                         median_null_percentile=g.null_percentile.median(), wilcoxon_null_p=wilcoxon(g.null_percentile - .5).pvalue,
                         median_lineage_null_percentile=g.lineage_null_percentile.median(), wilcoxon_lineage_null_p=wilcoxon(g.lineage_null_percentile - .5).pvalue))
    esum = []
    for (cm, e), g in er.groupby(['cm', 'edge']):
        esum.append(dict(cm=cm, edge=e, signature_r=g.signature_r.iloc[0], median_spearman=g.spearman.median(),
                         median_signature_matched_expected=g.signature_matched_expected.median(), median_residual=g.residual.median(),
                         n_sections_residual_gt_0=int((g.residual > 0).sum()), n_sections=len(g),
                         sign_test_p=binomtest(int((g.residual > 0).sum()), len(g)).pvalue,
                         median_percentile_among_signature_matched=g.percentile_among_signature_matched.median()))
    pd.DataFrame(summ).to_csv(out / 'visium_signature_adjusted_summary_v1.tsv', sep='\t', index=False)
    pd.DataFrame(esum).to_csv(out / 'visium_signature_adjusted_edge_summary_v1.tsv', sep='\t', index=False)
    print(pd.DataFrame(summ).round(4).to_string()); print(pd.DataFrame(esum).round(4).to_string())


if __name__ == '__main__':
    main()
