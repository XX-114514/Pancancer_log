# Annotation v2 step A3 helper (2026-11-06): PROPOSE a major label per cluster of every cancer type from the A2
# evidence, for manual review (the reviewed labels live in inventories/annotation_v2/major_labels_reviewed.tsv).
# Rules (methods/tme-annotation-v2-lodi-workflow-20261105.md):
#  lineage score = highest detection fraction among the lineage's identity genes (measured cells only);
#  lineage is supported when one of its genes is a positive DE marker (avg_log2FC>=.25, pct diff>=.25, adj p<.01)
#  1 low quality : median genes below the object's 10th percentile AND no supported lineage
#  2 doublet-like: two different supported lineages with score >= 0.4 AND cluster doublet evidence z >= 1
#                  (only one of the two conditions -> 'review_two_lineages')
#  3 tissue-specific population (Lodi): its module has >= 2 DE genes in the matching organ context
#  4 otherwise the supported lineage with the highest score; nothing supported -> 'unresolved'
# Usage: propose_major_labels_v1.py MAJOR_ROOT OUT.tsv
import sys, os, glob
import numpy as np, pandas as pd

root, outf = sys.argv[1:3]
LIN = {'T': ['CD3D', 'CD3E', 'CD3G', 'TRAC'], 'NK': ['KLRF1', 'NCR1', 'KLRD1'], 'B': ['MS4A1', 'CD19', 'CD22', 'BANK1'],
       'Plasma': ['MZB1', 'JCHAIN', 'DERL3', 'TNFRSF17'], 'Mac_Mono': ['CD68', 'CD14', 'C1QA', 'FCN1', 'CSF1R', 'LYZ'],
       'Neutrophil': ['CSF3R', 'FCGR3B', 'CXCR2'], 'pDC': ['LILRA4', 'CLEC4C'], 'cDC': ['CLEC9A', 'XCR1', 'CD1C', 'FCER1A', 'LAMP3'],
       'Mast': ['TPSAB1', 'CPA3', 'MS4A2'], 'EC': ['PECAM1', 'VWF', 'CLDN5', 'CDH5'], 'Lymphatic_EC': ['PROX1', 'CCL21', 'MMRN1'],
       'Fibroblast': ['COL1A1', 'DCN', 'LUM', 'PDGFRA'], 'Mural': ['RGS5', 'MYH11', 'NOTCH3', 'CSPG4'],
       'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT5'], 'Mesothelial': ['MSLN', 'UPK3B', 'CALB2'],
       'Melanocytic': ['PMEL', 'MLANA', 'TYR', 'DCT'], 'Glial_neural': ['GFAP', 'PLP1', 'OLIG2', 'MBP', 'SNAP25'],
       'Hepatocyte': ['ALB', 'APOA1', 'HP', 'TTR'], 'Osteo_chondro': ['COL2A1', 'SP7', 'IBSP', 'ACAN'],
       'Adipocyte': ['ADIPOQ', 'PLIN1'], 'Skeletal_muscle': ['ACTA1', 'CKM', 'MYH2'], 'Erythroid': ['HBB', 'HBA1', 'HBA2']}
TISSUE = {'Kupffer_like': 'Kupffer', 'Microglia_like': 'Microglia', 'Alveolar_macrophage_like': 'Alveolar_macrophage',
          'Alveolar_EC_like': 'Alveolar_EC', 'Sinusoidal_EC_like': 'Sinusoidal_EC'}
out = []
for d in sorted(glob.glob(f'{root}/*/clusters.tsv')):
    od = os.path.dirname(d); cancer = os.path.basename(od)
    if cancer.startswith('_'):
        continue
    S = pd.read_csv(d, sep='\t', index_col=0); S.index = S.index.astype(str)
    E = pd.read_csv(f'{od}/expression.tsv', sep='\t', dtype={'cluster': str}).pivot(index='cluster', columns='gene', values='pct_expressing_measured')
    DE = pd.read_csv(f'{od}/de.tsv', sep='\t', dtype={'cluster': str})
    ok = DE[(DE.avg_log2FC >= .25) & ((DE['pct.1'] - DE['pct.2']) >= .25) & (DE.p_val_adj < .01)]
    pos = ok.groupby('cluster')['gene'].agg(set)
    EV = pd.read_csv(f'{od}/evidence/all_panel_evidence.tsv', sep='\t', dtype={'cluster': str})
    g10 = S['genes_pct_of_object_p10']
    for c, r in S.iterrows():
        de = pos.get(c, set())
        sc_ = {k: float(E.loc[c, [x for x in v if x in E.columns]].max()) for k, v in LIN.items() if any(x in E.columns for x in v)}
        sup = {k: s for k, s in sc_.items() if set(LIN[k]) & de}
        top = sorted(sup.items(), key=lambda x: -x[1])
        tis = EV[(EV.cluster == c) & EV.panel.isin(TISSUE) & (EV.organ_context == 'context_present') & EV.positive_DEG_genes.notna()]
        tis = tis[tis.positive_DEG_genes.str.count(',') >= 1]
        flags = []
        if bool(g10[c]) and not top:
            lab = 'Low_quality'
        elif len(top) >= 2 and top[1][1] >= 0.4 and not ({top[0][0], top[1][0]} <= {'T', 'NK'} or {top[0][0], top[1][0]} <= {'B', 'Plasma'}
                                                         or {top[0][0], top[1][0]} <= {'Mac_Mono', 'cDC'} or {top[0][0], top[1][0]} <= {'EC', 'Lymphatic_EC'}
                                                         or {top[0][0], top[1][0]} <= {'Fibroblast', 'Mural'}):
            lab = 'Doublet_like' if (r.get('doublet_score_z', 0) or 0) >= 1 else 'review_two_lineages'
            flags.append(f'{top[0][0]}+{top[1][0]}')
        elif len(tis):
            lab = TISSUE[tis.sort_values('positive_DEG_genes', key=lambda s: s.str.count(','), ascending=False).panel.iloc[0]]
        elif top:
            lab = top[0][0]
        else:
            lab = 'unresolved'
        if r['uncertain_share'] >= 0.5:
            flags.append('uncertain_majority')
        if r['top_source_share'] >= 0.8 or r['sources'] < 3:
            flags.append('source_restricted')
        out.append({'cancer': cancer, 'cluster': c, 'cells': r['cells'], 'non_malignant': r['non_malignant'], 'uncertain_share': r['uncertain_share'],
                    'global_c13_share': r['global_c13_share'], 'proposed_major': lab, 'flags': ';'.join(flags),
                    'lineage_scores_supported': ';'.join(f'{k}:{v:.2f}' for k, v in top[:4]),
                    'lineage_scores_all_ge0.3': ';'.join(f'{k}:{v:.2f}' for k, v in sorted(sc_.items(), key=lambda x: -x[1]) if v >= 0.3),
                    'doublet_score_z': r.get('doublet_score_z'), 'median_genes': r['median_genes'], 'median_pct_mt': r['median_pct_mt'],
                    'context_modules_with_DEG': r.get('context_modules_with_DEG'), 'top_markers': r['top_markers']})
P = pd.DataFrame(out)
P.to_csv(outf, sep='\t', index=False)
print(P.groupby('proposed_major').agg(clusters=('cells', 'size'), cells=('cells', 'sum'), non_malignant=('non_malignant', 'sum')).sort_values('cells', ascending=False).to_string())
