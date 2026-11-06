# Annotation v2 step A3 (2026-11-06): reviewed major label per cluster of each cancer type.
# Replaces the first proposal rules after reading all ambiguous clusters: promiscuous genes removed from the
# competing-lineage sets; a lineage needs DE support to count as a SECOND lineage, while the dominant lineage may be
# supported by detection >= 0.5 alone (in T-rich objects T genes are not differential); pair resolvers and tissue
# population gene requirements below; MANUAL holds cluster-level overrides from reading the markers, each with a reason.
# major groups: T_NK, B_Plasma (Lodi pool these); Myeloid pool = Mac_Mono, Neutrophil, cDC, pDC, Mast (DC and mast are split
# further at higher resolution in the lineage step, as Lodi do); tissue-specific populations stay separate.
# Writes reviewed labels and, per label, how many non-malignant cells of global cluster 13 landed there.
# Usage: review_major_labels_v1.py MAJOR_ROOT OUT_TSV OUT_C13_TSV
import sys, os, glob
import numpy as np, pandas as pd, anndata as ad

root, outf, outc = sys.argv[1:4]
LIN = {'T_NK': ['CD3D', 'CD3E', 'CD3G', 'TRAC', 'KLRF1', 'NCR1', 'KLRD1', 'NKG7'], 'B_Plasma': ['MS4A1', 'CD19', 'BANK1', 'MZB1', 'TNFRSF17', 'DERL3'],
       'Mac_Mono': ['CD68', 'CD14', 'C1QA', 'FCN1', 'CSF1R'], 'Neutrophil': ['FCGR3B', 'CXCR2'], 'pDC': ['LILRA4', 'CLEC4C'],
       'cDC': ['CLEC9A', 'XCR1', 'CD1C', 'FCER1A', 'LAMP3', 'CD207'], 'Mast': ['TPSAB1', 'CPA3', 'MS4A2'], 'EC': ['PECAM1', 'VWF', 'CLDN5', 'CDH5'],
       'Lymphatic_EC': ['CCL21', 'MMRN1'], 'Fibroblast': ['COL1A1', 'DCN', 'LUM', 'PDGFRA'], 'Mural': ['RGS5', 'MYH11', 'NOTCH3', 'CSPG4'],
       'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'KRT5'], 'Mesothelial': ['MSLN', 'UPK3B', 'ITLN1'],
       'Melanocytic': ['PMEL', 'MLANA', 'TYR', 'DCT'], 'Glial_neural': ['GFAP', 'PLP1', 'OLIG2'], 'Hepatocyte': ['ALB', 'APOA1', 'HP', 'TTR'],
       'Osteo_chondro': ['COL2A1', 'SP7', 'IBSP', 'ACAN'], 'Adipocyte': ['ADIPOQ', 'PLIN1'], 'Skeletal_muscle': ['ACTA1', 'CKM', 'MYH2'], 'Erythroid': ['HBB', 'HBA1', 'HBA2']}
GROUP = {'T_NK': 'immune', 'B_Plasma': 'immune', 'Mac_Mono': 'immune', 'Neutrophil': 'immune', 'pDC': 'immune', 'cDC': 'immune', 'Mast': 'immune',
         'Kupffer': 'immune', 'Microglia': 'immune', 'Alveolar_macrophage': 'immune', 'Osteoclast_like': 'immune', 'Myeloid_mixed': 'immune',
         'EC': 'endothelial', 'Lymphatic_EC': 'endothelial', 'Sinusoidal_EC': 'endothelial', 'Alveolar_EC': 'endothelial',
         'Fibroblast': 'stromal', 'Mural': 'stromal', 'Osteo_chondro': 'stromal', 'Adipocyte': 'stromal', 'Skeletal_muscle': 'stromal',
         'Epithelial': 'parenchymal', 'Hepatocyte': 'parenchymal', 'Mesothelial': 'parenchymal', 'Melanocytic': 'parenchymal', 'Glial_neural': 'parenchymal',
         'Low_quality': 'excluded', 'Doublet_like': 'excluded', 'Erythroid': 'excluded', 'Cycling_mixed': 'split_at_higher_resolution', 'unresolved': 'review'}
CYC = {'MKI67', 'TOP2A', 'STMN1', 'UBE2C', 'CENPF', 'TYMS', 'NUSAP1'}
# (cancer, cluster): (label, reason) from reading markers of the clusters listed in the run log
MANUAL = {
    ('Renal_cell_carcinoma', '16'): ('Mac_Mono', 'cycling macrophages (C1Q, APOE + TUBA1B, STMN1); high doublet score explained by cycling, no second-lineage markers'),
    ('Renal_cell_carcinoma', '6'): ('T_NK', 'mostly CD8 T (TOX, CD8A, CCL5) with tumour RNA signal (KRT18, NNMT, NDUFA4L2) and 1,393 scATOMIC kidney cancer calls; re-check in T lineage QC round'),
    ('Pancreatic_neuroendocrine_tumor', '24'): ('Epithelial', 'neuroendocrine tumour-like (CHGB, CPE, APLP1), 99.7% uncertain; not a doublet'),
    ('Pancreatic_neuroendocrine_tumor', '16'): ('Epithelial', 'neuroendocrine genes (PCSK1N, CHGB, SCG5) at low depth, uncertain 0.54'),
    ('Pancreatic_neuroendocrine_tumor', '22'): ('Low_quality', 'acinar transcripts (PLA2G1B, CTRC, CEL) at 380 genes: ambient-like'),
    ('Pancreatic_neuroendocrine_tumor', '25'): ('cDC', 'Langerhans-like DC (CD207, CD1A, FCER1A)'),
    ('Endometrial_cancer', '20'): ('Epithelial', 'ciliated epithelium (CAPS, RSPH1, WFDC2); doublet score alone does not establish a doublet'),
    ('Renal_medullary_carcinoma', '20'): ('pDC', 'LILRA4, CLEC4C, IL3RA, IRF7'),
    ('Head_and_neck_squamous_cell_carcinoma', '23'): ('Hepatocyte', 'APOA1, TF, HAMP, CPS1, 100% uncertain, 61 cells; hepatocyte-like cells in a site-unusual sample, review source'),
    ('Head_and_neck_squamous_cell_carcinoma', '10'): ('Doublet_like', 'COL1A1/TAGLN + KRT18/KRT15, 49 cells'),
    ('Osteosarcoma', '22'): ('Osteoclast_like', 'CTSK, MMP9, ACP5, ATP6V1 (41 cells, doublet score high: osteoclasts are multinucleate; review in myeloid pool)'),
    ('Osteosarcoma', '15'): ('Epithelial', 'alveolar epithelium (NAPSA, SFTPD, SFTA2), 26 cells: lung metastasis context'),
    ('Osteosarcoma', '5'): ('Low_quality', '683 genes, no specific markers'),
    ('Soft_tissue_sarcoma', '17'): ('Osteoclast_like', 'MMP9, CTSK, ACP5, ATP6V1 family'),
    ('Diffuse_glioma', '1'): ('Glial_neural', 'glioma-lineage programme (CHI3L1, PTN, SCRG1, CLU), 97% uncertain'),
    ('Diffuse_glioma', '23'): ('Glial_neural', 'oligodendrocyte genes (MBP, CNP, ERMN, TF), 90% uncertain'),
    ('Diffuse_glioma', '24'): ('Glial_neural', 'glioma-lineage (CLU, CHI3L1, NRCAM, TNC), 99.7% uncertain; rule call Mesothelial was wrong'),
    ('Diffuse_glioma', '22'): ('Low_quality', 'HLA/C1Q/HBA2/MT3/GNLY mixed at 729 genes, 90% uncertain: ambient-like'),
    ('Diffuse_glioma', '6'): ('EC', 'CLDN5, VWF, ABCG2 with ordinary gene count: brain endothelium'),
    ('Cutaneous_melanoma', '20'): ('Melanocytic', 'neural-crest-like melanoma state (MPZ, NES, SCRG1), 83% uncertain; low PMEL does not exclude melanocytic origin (Tirosh 2016)'),
    ('Cutaneous_melanoma', '22'): ('Melanocytic', 'hypoxic programme (BNIP3, ENO1, VEGFA), 86% uncertain, no other lineage'),
    ('Cutaneous_melanoma', '19'): ('Epithelial', 'keratinocytes (KRT17, KRT10, SFN, DMKN)'),
    ('Uveal_melanoma', '2'): ('Melanocytic', 'hypoxic melanocytic (NDRG1, BNIP3; PMEL-lineage score 0.70), 99% uncertain'),
    ('Uveal_melanoma', '7'): ('Melanocytic', 'melanocytic score 0.91, 100% uncertain, MT 19%'),
    ('Uveal_melanoma', '13'): ('Melanocytic', 'melanocytic score 0.84, 100% uncertain'),
    ('Uveal_melanoma', '20'): ('Glial_neural', 'photoreceptor retina cells (NRL, ROM1, PDE6G, IMPG2)'),
    ('Basal_cell_carcinoma', '13'): ('Epithelial', 'CLDN1, TACSTD2, LYPD3; rule call Glial_neural was from a weak gene'),
    ('Mesothelioma', '10'): ('Mesothelial', 'ITLN1, CALB2, C1S, C1R, SLPI'), ('Mesothelioma', '12'): ('Mesothelial', 'ITLN1, CALB2, C1S'),
    ('Mesothelioma', '16'): ('Mesothelial', 'PRG4, ITLN1, TM4SF1, VCAM1'), ('Mesothelioma', '5'): ('Low_quality', '399 genes, no markers'),
    ('Ovarian_cancer', '11'): ('Fibroblast', 'TIMP1, BGN, CCDC80, C3 with KRT detection 0.51; possible mesothelial-like stroma, flag for review'),
    ('Papillary_thyroid_carcinoma', '26'): ('Doublet_like', 'RGS5 + PLVAP/FLT1, 6,162 genes'), ('Hepatocellular_carcinoma', '13'): ('Doublet_like', 'RGS5 + CLDN5/PLVAP, 4,414 genes'),
    ('Colorectal_cancer', '11'): ('Doublet_like', 'RGS5 + PLVAP/CLDN5, 3,577 genes'), ('Renal_cell_carcinoma', '23'): ('Doublet_like', 'RGS5 + PLVAP, 4,292 genes, doublet z 2.6'),
    ('Acral_melanoma', '6'): ('Doublet_like', 'RGS5 + PLVAP + COL4A1 with MKI67, 3,917 genes, doublet z 3.0'), ('Soft_tissue_sarcoma', '1'): ('Doublet_like', 'RGS5 + PLVAP, 4,735 genes, doublet z 2.2'),
    ('Diffuse_glioma', '0'): ('Glial_neural', 'OPC-like glioma (EGFR, PDGFRA, GPM6A, SOX4), 99.7% uncertain; PDGFRA made the rule call Fibroblast'),
    ('Diffuse_glioma', '25'): ('Glial_neural', 'glioma programme (CHI3L1, PTN, TTYH1), 100% uncertain; rule call Mural was wrong'),
    ('Basal_cell_carcinoma', '12'): ('Epithelial', 'CLDN4, LYPD3, AQP3, DMKN, 96% uncertain; S100A8/9 made the rule call Mac_Mono'),
    ('Basal_cell_carcinoma', '6'): ('Epithelial', 'DSC2, SDCBP2, SERPINB2, 99.6% uncertain; S100A9/IL1RN are not myeloid evidence here'),
    ('Basal_cell_carcinoma', '7'): ('Epithelial', 'CNFN, CSTA, DSC2, DMKN, LYPD3, 79% uncertain'),
    ('Acral_melanoma', '7'): ('Melanocytic', 'co-located amplicon genes (CHIC2+PDGFRA; TSPAN31+METTL1, CDK4 region; inferred tumour copy-number gain), 50% uncertain; not fibroblasts'),
    ('Cutaneous_melanoma', '23'): ('Doublet_like', 'melanocytic + pericyte genes (ABCC9, NDUFA4L2), 4,390 genes, 52% uncertain'),
    ('Intrahepatic_cholangiocarcinoma', '24'): ('Epithelial', 'KRT7, S100A14 are the only markers, 33% uncertain, c13 share 0.31; low depth'),
    ('Ovarian_cancer', '9'): ('T_NK', 'T detection 0.49, IFN markers MX1/MX2'),
    ('Ovarian_cancer', '2'): ('Mac_Mono', 'MHC-II, CCR7, GPR183 myeloid-like; DC split in the myeloid step'),
    ('Acral_melanoma', '18'): ('T_NK', 'T detection 0.48, 678 genes, no markers; re-check in T lineage QC round'),
    ('Cutaneous_melanoma', '26'): ('Doublet_like', 'RGS5 + CLDN5/PLVAP, 4,700 genes'), ('Pancreatic_neuroendocrine_tumor', '19'): ('Doublet_like', 'RGS5 + PLVAP/RAMP2, 4,504 genes, doublet z 2.9'),
}
rows, c13 = [], []
for d in sorted(glob.glob(f'{root}/*/clusters.tsv')):
    od = os.path.dirname(d); cancer = os.path.basename(od)
    if cancer.startswith('_'):
        continue
    S = pd.read_csv(d, sep='\t', index_col=0); S.index = S.index.astype(str)
    E = pd.read_csv(f'{od}/expression.tsv', sep='\t', dtype={'cluster': str}).pivot(index='cluster', columns='gene', values='pct_expressing_measured')
    DE = pd.read_csv(f'{od}/de.tsv', sep='\t', dtype={'cluster': str})
    pos = DE[(DE.avg_log2FC >= .25) & ((DE['pct.1'] - DE['pct.2']) >= .25) & (DE.p_val_adj < .01)].groupby('cluster')['gene'].agg(set)
    for c, r in S.iterrows():
        de = pos.get(c, set())
        sc_ = {k: float(E.loc[c, [x for x in v if x in E.columns]].max()) for k, v in LIN.items() if any(x in E.columns for x in v)}
        sup = sorted(((k, s) for k, s in sc_.items() if set(LIN[k]) & de), key=lambda x: -x[1])
        pres = sorted(((k, s) for k, s in sc_.items() if s >= 0.5), key=lambda x: -x[1])
        cyc = len(CYC & set(str(r['top_markers']).split(',')[:12])) >= 3
        low = bool(r['genes_pct_of_object_p10']) and not sup and not (pres and pres[0][1] >= 0.5)
        reason, lab = '', None
        if (cancer, c) in MANUAL:
            lab, reason = MANUAL[(cancer, c)]; basis = 'manual'
        else:
            basis = 'rule'
            if low:
                lab, reason = 'Low_quality', 'genes below object p10, no lineage support'
            elif sup or pres:
                first = sup[0] if sup else pres[0]
                second = next((x for x in sup if x[0] != first[0]), None)
                if not sup:
                    reason = f'{first[0]} detection {first[1]:.2f} without DE (lineage dominant in this object)'
                else:
                    reason = f'{first[0]} detection {first[1]:.2f} with DE support'
                lab = first[0]
                if second and second[1] >= 0.4:
                    pair = {first[0], second[0]}
                    if pair == {'B_Plasma', 'pDC'} and 'LILRA4' in de:
                        lab, reason = 'pDC', 'LILRA4 marker; MZB1 shared with pDC'
                    elif pair == {'Mac_Mono', 'cDC'}:
                        lab, reason = 'Mac_Mono', 'myeloid; DCs are split in the myeloid lineage step'
                    elif pair == {'EC', 'Lymphatic_EC'}:
                        lab, reason = 'Lymphatic_EC', 'CCL21/MMRN1 markers'
                    elif pair == {'Fibroblast', 'Mural'}:
                        pass
                    elif pair == {'EC', 'Mural'}:
                        if r['median_genes'] >= 1.5 * S['median_genes'].median() or (r['doublet_score_z'] or 0) >= 1:
                            lab, reason = 'Doublet_like', 'EC and mural markers with high complexity or doublet score'
                    elif 'Epithelial' in pair and pair & {'Hepatocyte', 'Mesothelial', 'Melanocytic', 'Mac_Mono'}:
                        other = (pair - {'Epithelial'}).pop()
                        lab = other if other == 'Hepatocyte' else 'Epithelial'
                        reason = f'{lab}: {other} genes shared with this epithelium or ambient'
                    elif (r['doublet_score_z'] or 0) >= 1 and not cyc:
                        lab, reason = 'Doublet_like', f'{first[0]}+{second[0]} markers, doublet z {r["doublet_score_z"]}'
                    else:
                        reason += f'; second lineage {second[0]} {second[1]:.2f} without doublet evidence'
            else:
                lab, reason = 'unresolved', 'no lineage support'
            # tissue-specific populations need the parent and their key genes as markers
            if lab == 'Mac_Mono':
                if cancer in ('Hepatocellular_carcinoma', 'Intrahepatic_cholangiocarcinoma') and 'VSIG4' in de and len(de & {'CD163', 'MAFB', 'LYVE1', 'CD209', 'PLTP', 'SELENOP', 'MERTK', 'FOLR2'}) >= 1:
                    lab, reason = 'Kupffer', 'VSIG4 + Kupffer module marker'
                elif cancer == 'NSCLC' and 'MARCO' in de and len(de & {'FABP4', 'PPARG', 'INHBA', 'ALDH2', 'FBP1'}) >= 1:
                    lab, reason = 'Alveolar_macrophage', 'MARCO + alveolar macrophage module marker'
                elif cancer == 'Diffuse_glioma' and len(de & {'P2RY12', 'TMEM119', 'P2RY13', 'SLC2A5', 'CX3CR1', 'GPR34', 'BHLHE41'}) >= 2:
                    lab, reason = 'Microglia', 'microglia markers (P2RY12/P2RY13/CX3CR1/GPR34/BHLHE41 at least 2)'
            if lab == 'EC' and cancer in ('Hepatocellular_carcinoma', 'Intrahepatic_cholangiocarcinoma') and len(de & {'CLEC4G', 'CLEC4M', 'DNASE1L3', 'FCN2', 'CLEC1B', 'CCL14'}) >= 1:
                lab, reason = 'Sinusoidal_EC', 'liver sinusoidal EC marker'
            if lab == 'EC' and cancer == 'NSCLC' and len(de & {'EDNRB', 'IL1RL1', 'HPGD'}) >= 2:
                lab, reason = 'Alveolar_EC', 'alveolar EC markers'
            if cyc and sup and lab not in ('Doublet_like',):
                top2 = sorted(sc_.items(), key=lambda x: -x[1])[:2]
                if top2[0][1] < 0.6 or (len(top2) > 1 and top2[1][1] > 0.75 * top2[0][1]):
                    lab, reason = 'Cycling_mixed', 'cycling programme with several lineages; split at higher resolution'
            elif cyc and not sup and (not pres or pres[0][1] < 0.6):
                lab, reason = 'Cycling_mixed', 'cycling programme with no dominant lineage; split at higher resolution'
        fl = []
        if r['uncertain_share'] >= 0.5: fl.append('uncertain_majority')
        if r['top_source_share'] >= 0.8 or r['sources'] < 3: fl.append('source_restricted')
        if cyc: fl.append('cycling')
        if r['stress_score_mean'] > S['stress_score_mean'].quantile(.9): fl.append('stress_high')
        conf = 'low' if lab in ('unresolved', 'Cycling_mixed') or (basis == 'rule' and not sup) else ('medium' if 'source_restricted' in fl or 'uncertain_majority' in fl else 'high')
        rows.append({'cancer': cancer, 'cluster': c, 'cells': r['cells'], 'non_malignant': r['non_malignant'], 'uncertain_share': r['uncertain_share'],
                     'major_label': lab, 'group': GROUP[lab], 'confidence': conf, 'basis': basis, 'reason': reason, 'flags': ';'.join(fl),
                     'median_genes': r['median_genes'], 'median_pct_mt': r['median_pct_mt'], 'doublet_z': r['doublet_score_z'], 'top_markers': str(r['top_markers'])[:150]})
    L = pd.DataFrame([x for x in rows if x['cancer'] == cancer]).set_index('cluster')['major_label']
    A = ad.read_h5ad(f'{od}/object.h5ad', backed='r'); o = A.obs
    cl = o['leiden_r1'].astype(str).values
    cell = pd.DataFrame({'cell_id': o.index.astype(str), 'cancer': cancer, 'cluster': cl, 'major_label': pd.Series(cl).map(L).values, 'split_from': ''})
    # mixed cycling clusters: each cell goes to the lineage with the highest mean identity-gene expression
    mix = np.flatnonzero(cell['major_label'].values == 'Cycling_mixed')
    if len(mix):
        B = A[mix].to_memory()
        names = [k for k in LIN if k not in ('Erythroid', 'Neutrophil', 'pDC', 'cDC', 'Lymphatic_EC', 'Adipocyte', 'Skeletal_muscle', 'Osteo_chondro')]
        M = np.column_stack([np.asarray(B[:, [g for g in LIN[k] if g in B.var_names]].X.mean(1)).ravel() for k in names])
        best = np.array(names, dtype=object)[M.argmax(1)]; best[M.max(1) == 0] = 'Cycling_unassigned'
        cell.loc[mix, 'major_label'] = best; cell.loc[mix, 'split_from'] = 'Cycling_mixed'
    cell.to_csv(f'{outf}.cells.{cancer}.tsv.gz', sep='\t', index=False, compression='gzip')
    x = o[(o['leiden_tme'].astype(str) == '13')]
    lab = cell.set_index('cell_id')['major_label'].reindex(x.index.astype(str)).values
    c13.append(pd.DataFrame({'cancer': cancer, 'major_label': lab, 'non_malignant': (x['tme_flag'].astype(str) == 'non_malignant').values}))
    print(cancer, flush=True)
P = pd.DataFrame(rows); P.to_csv(outf, sep='\t', index=False)
C = pd.concat(c13)
T = C.groupby('major_label').agg(c13_cells=('non_malignant', 'size'), c13_non_malignant=('non_malignant', 'sum')).sort_values('c13_non_malignant', ascending=False)
T.to_csv(outc, sep='\t')
pd.set_option('display.width', 200)
print(P.groupby(['group', 'major_label']).agg(clusters=('cells', 'size'), cells=('cells', 'sum'), non_malignant=('non_malignant', 'sum')).to_string())
print(P.groupby('confidence').cells.sum().to_dict(), P.groupby('basis').cells.size().to_dict())
print(T.to_string())
print(C.groupby('cancer').non_malignant.sum().sort_values(ascending=False).head(12).to_string())
