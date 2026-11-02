#!/usr/bin/env python3
"""Aggregate CT75 spatial validation outputs into publishable tables and figures (no cell or spot rows).

Usage: report_spatial_cm_v1.py --run RUN_DIR --fit FIT_RUN --output DIR
"""
import argparse, glob, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

MEMBERS = {'CM01': ['N_44', 'N_55', 'N_62', 'N_114'], 'CM02': ['N_176', 'N_243', 'N_234', 'N_238', 'N_44']}


def unit_of(section):
    if section.startswith('S07_'):
        return 'S07_breast'
    return section


def summarise(scores, unit_col):
    rows = []
    for (cm, stat), g in scores.groupby(['cm', 'statistic']):
        u = g.groupby(unit_col)[['score', 'null_percentile', 'null_median', 'lineage_null_percentile', 'lineage_null_median']].mean().dropna(subset=['score'])
        rows.append(dict(cm=cm, statistic=stat, n_units=len(u), median_score=u.score.median(),
                         n_units_score_gt_0_5=int((u.score > .5).sum()), median_null_score=u.null_median.median(),
                         median_null_percentile=u.null_percentile.median(), n_units_null_percentile_gt_0_5=int((u.null_percentile > .5).sum()),
                         wilcoxon_score_vs_0_5_p=wilcoxon(u.score - .5).pvalue if len(u) >= 5 else np.nan,
                         wilcoxon_null_percentile_vs_0_5_p=wilcoxon(u.null_percentile - .5).pvalue if len(u) >= 5 else np.nan,
                         median_lineage_null_score=u.lineage_null_median.median(), median_lineage_null_percentile=u.lineage_null_percentile.median(),
                         n_units_lineage_null_percentile_gt_0_5=int((u.lineage_null_percentile > .5).sum()),
                         wilcoxon_lineage_null_percentile_vs_0_5_p=wilcoxon(u.lineage_null_percentile.dropna() - .5).pvalue if u.lineage_null_percentile.notna().sum() >= 5 else np.nan))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True); ap.add_argument('--fit', required=True); ap.add_argument('--output', required=True)
    a = ap.parse_args()
    R = Path(a.run); out = Path(a.output); out.mkdir(parents=True, exist_ok=True); F = out / 'figures'; F.mkdir(exist_ok=True)
    cat = pd.read_csv(R.parent / '20261101_TME_CT75_UMAP_CM_review_v1/outputs/CT75_node_catalogue_v1.tsv', sep='\t').set_index('node')
    lab = lambda n: f"{n} | {cat.loc[n, 'identity'].replace('_reference_related', '')} | C{cat.loc[n, 'numeric_cluster']:02d}"
    mem = list(dict.fromkeys(MEMBERS['CM01'] + MEMBERS['CM02']))
    summary = {}

    # 1. annotation transfer validation
    rows, f1 = [], {}
    names = {}
    for f in sorted(glob.glob(str(R / 'annotation/panel_*_heldout_summary.json'))):
        s = json.load(open(f)); t = pd.read_csv(f.replace('_summary.json', '_node_F1.tsv'), sep='\t').set_index('node')
        src = s['first_section'].split('__')[0].replace('_GSM9054471', '')
        src = {'S10_10x_IO_ovarian': 'S10-S12 Xenium IO (480)'}.get(src, src)
        names[s['panel']] = src
        rows.append(dict(panel=src, n_panel_genes=s['n_panel_genes'], parent_accuracy=s['parent_accuracy'],
                         parent_cluster_accuracy=s['parent_cluster_accuracy'], node_accuracy=s['node_accuracy_all'],
                         n_CT75_nodes_F1_ge_0_5=int((t.F1 >= .5).sum()))); f1[src] = t.F1
    fg = R / 'diagnostics/fullgene_heldout_summary.json'
    if fg.exists():
        s = json.load(open(fg)); t = pd.read_csv(R / 'diagnostics/fullgene_heldout_node_F1.tsv', sep='\t').set_index('node')
        rows.append(dict(panel='All 12,783 genes (upper bound)', n_panel_genes=s['n_genes'], parent_accuracy=s['parent_accuracy'],
                         parent_cluster_accuracy=s['parent_cluster_accuracy'], node_accuracy=s['node_accuracy_all'],
                         n_CT75_nodes_F1_ge_0_5=int((t.F1 >= .5).sum()))); f1['All 12,783 genes (upper bound)'] = t.F1
    val = pd.DataFrame(rows); val.to_csv(out / 'annotation_transfer_validation.tsv', sep='\t', index=False)
    F1 = pd.DataFrame(f1); F1.index.name = 'node'
    F1.loc[mem].rename(index=lab).to_csv(out / 'CM_member_heldout_F1.tsv', sep='\t')
    F1.to_csv(out / 'CT75_node_heldout_F1_all_panels.tsv', sep='\t')
    fig, ax = plt.subplots(figsize=(1.6 + 1.3 * F1.shape[1], 4.2))
    im = ax.imshow(F1.loc[mem].to_numpy(), cmap='viridis', vmin=0, vmax=1, aspect='auto')
    for i in range(len(mem)):
        for j in range(F1.shape[1]):
            v = F1.loc[mem].iat[i, j]; ax.text(j, i, f'{v:.2f}', ha='center', va='center', color='white' if v < .6 else 'black', fontsize=8)
    ax.set_yticks(range(len(mem)), [lab(n) for n in mem], fontsize=8); ax.set_xticks(range(F1.shape[1]), F1.columns, rotation=30, ha='right', fontsize=8)
    ax.set_title('Held-out atlas cells: CT75-scheme annotation F1 per CM member\n(pre-registered measurability threshold F1 >= 0.5)', fontsize=9)
    fig.colorbar(im, ax=ax, shrink=.7, label='F1'); fig.tight_layout(); fig.savefig(F / '01_CM_member_annotation_F1.png', dpi=170); plt.close(fig)
    summary['annotation'] = val.to_dict('records')

    # 2. Xenium neighbourhoods
    for tag in ['primary_F1ge0.5', 'exploratory_allmembers']:
        p = R / f'neighbourhood_{tag}/neighbourhood_CM_section_scores_v1.tsv'
        if not p.exists():
            continue
        sc = pd.read_csv(p, sep='\t'); sc['unit'] = [unit_of(s) if not str(pt).startswith('GSE300147') else pt for s, pt in zip(sc.section, sc.patient)]
        sc['source'] = sc.section.str.split('_').str[0]
        sc.groupby(['source', 'unit', 'cm', 'statistic'])[['n_members_used', 'score', 'null_median', 'null_percentile', 'lineage_null_median', 'lineage_null_percentile']].mean().reset_index().to_csv(out / f'xenium_{tag}_unit_scores.tsv', sep='\t', index=False)
        sm = summarise(sc, 'unit'); sm.to_csv(out / f'xenium_{tag}_summary.tsv', sep='\t', index=False)
        summary[f'xenium_{tag}'] = sm.to_dict('records')
        ed = pd.read_csv(R / f'neighbourhood_{tag}/neighbourhood_CM_edge_section_values_v1.tsv', sep='\t')
        ed['unit'] = [unit_of(s) if not str(pt).startswith('GSE300147') else pt for s, pt in zip(ed.section, ed.patient)]
        e = ed.groupby(['cm', 'edge', 'unit'])[['nhood_z', 'pseudospot_spearman']].mean().reset_index()
        es = e.groupby(['cm', 'edge']).agg(n_units=('nhood_z', lambda v: int(v.notna().sum())), median_nhood_z=('nhood_z', 'median'),
                                          n_units_z_gt_0=('nhood_z', lambda v: int((v > 0).sum())), median_pseudospot_spearman=('pseudospot_spearman', 'median')).reset_index()
        es.to_csv(out / f'xenium_{tag}_edge_summary.tsv', sep='\t', index=False)
        if tag.startswith('exploratory'):
            u = sc.groupby(['source', 'unit', 'cm', 'statistic'])[['score', 'null_median']].mean().reset_index()
            stats = sorted(u.statistic.unique()); fig, axs = plt.subplots(1, len(stats), figsize=(4.2 * len(stats), 3.8), sharey=True)
            for ax, st in zip(np.atleast_1d(axs), stats):
                g = u[u.statistic == st]
                for k, cm in enumerate(['CM01', 'CM02']):
                    h = g[g.cm == cm]; x = k + np.random.default_rng(k).uniform(-.12, .12, len(h))
                    for src, col in zip(sorted(u.source.unique()), plt.cm.tab10.colors):
                        m = h.source.to_numpy() == src; ax.scatter(x[m], h.score[m], s=22, color=col, label=src if k == 0 else None)
                    ax.scatter([k + .3] * len(h), h.null_median, s=10, marker='_', color='grey')
                ax.axhline(.5, ls='--', c='k', lw=.7); ax.set_xticks([0, 1], ['CM01', 'CM02']); ax.set_title(st.replace('_', ' '), fontsize=9)
            np.atleast_1d(axs)[0].set_ylabel('CoVarNet-style score (grey: random-member null median)'); np.atleast_1d(axs)[-1].legend(fontsize=7, frameon=False)
            fig.suptitle('Xenium cellular neighbourhoods - EXPLORATORY (all members, regardless of annotation F1)', fontsize=10); fig.tight_layout()
            fig.savefig(F / '03_xenium_exploratory_neighbourhood_scores.png', dpi=170); plt.close(fig)

    # 3. Visium
    vs = R / 'visium_scores/visium_CM_section_scores_v1.tsv'
    if vs.exists():
        sc = pd.read_csv(vs, sep='\t'); sm = summarise(sc, 'sample'); sm.to_csv(out / 'visium_summary.tsv', sep='\t', index=False)
        summary['visium'] = sm.to_dict('records')
        sc.drop(columns=['n_spots']).to_csv(out / 'visium_section_scores.tsv', sep='\t', index=False)
        es = pd.read_csv(R / 'visium_scores/visium_CM_edge_summary_v1.tsv', sep='\t'); es.to_csv(out / 'visium_edge_summary.tsv', sep='\t', index=False)
        summary['visium_edges'] = es.to_dict('records')
        stats = sorted(sc.statistic.unique()); fig, axs = plt.subplots(1, len(stats), figsize=(4.5 * len(stats), 3.8), sharey=True)
        for ax, st in zip(np.atleast_1d(axs), stats):
            g = sc[sc.statistic == st]
            for k, cm in enumerate(['CM01', 'CM02']):
                h = g[g.cm == cm]; x = k + np.random.default_rng(k).uniform(-.12, .12, len(h))
                ax.scatter(x, h.score, s=20, color='#3979B5'); ax.scatter([k + .3] * len(h), h.null_median, s=10, marker='_', color='grey')
            ax.axhline(.5, ls='--', c='k', lw=.7); ax.set_xticks([0, 1], ['CM01', 'CM02']); ax.set_title(st.replace('_', ' '), fontsize=9)
        np.atleast_1d(axs)[0].set_ylabel('Score per patient section (grey: null median)')
        fig.suptitle('GSE273952 ccRCC Visium, cell2location (CoVarNet method), 19 patients', fontsize=10); fig.tight_layout()
        fig.savefig(F / '02_visium_CM_scores.png', dpi=170); plt.close(fig)
        act = pd.read_csv(R / 'visium_scores/visium_CM_activity_v1.tsv.gz', sep='\t', index_col=0)
        spots = pd.concat([pd.read_csv(d / 'visium_spots_v1.tsv.gz', sep='\t', index_col=0) for d in sorted((R / 'visium_c2l/sections').glob('GSM*'))])
        show = sorted(spots['sample'].unique())[:4]
        fig, axs = plt.subplots(2, len(show), figsize=(3.3 * len(show), 6.4))
        for j, s in enumerate(show):
            ix = spots.index[spots['sample'] == s]
            for i, cm in enumerate(['CM01', 'CM02']):
                ax = axs[i, j]; ax.scatter(spots.loc[ix, 'pxl_col'], -spots.loc[ix, 'pxl_row'], c=act.loc[ix, cm].clip(0, 1), s=3, cmap='magma', vmin=0, vmax=1)
                ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f'{s.split("_", 1)[1]} {cm}', fontsize=8)
        fig.suptitle('CM activity (W-weighted member abundance; 99th percentile = 1)', fontsize=10); fig.tight_layout()
        fig.savefig(F / '04_visium_CM_activity_maps.png', dpi=150); plt.close(fig)
        sig = R / 'visium_c2l/reference_signatures_v1.tsv'
        if sig.exists():
            S = pd.read_csv(sig, sep='\t', index_col=0)[mem]; L = np.log1p(S / S.sum() * 1e4)
            C = np.corrcoef(L.T.to_numpy()); pd.DataFrame(C, index=[lab(n) for n in mem], columns=mem).to_csv(out / 'CM_member_reference_signature_correlation.tsv', sep='\t')
    json.dump(summary, open(out / 'summary.json', 'w'), indent=2, default=float)


if __name__ == '__main__':
    main()
