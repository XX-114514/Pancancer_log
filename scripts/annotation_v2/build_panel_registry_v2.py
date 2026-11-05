# Build the merged annotation panel registry v2 (2026-11-05).
#   base     : resources/annotation/v7_panel_repair_v1/registry.json (112 modules, unchanged)
#   + Lodi   : author-code panels from lambrechtslab/Pan-cancer_Lodietal (commit 35f054b): the 29-gene major-type
#              check (9 modules), 5 cancer-type-specific populations, 2 explicit state signatures, 3 named external
#              signatures (gene lists not in the repository), 21 published DEG profiles (marker_genes.txt)
#   + v1     : TME de novo annotation v1 modules that V7 lacks (organ epithelium, bone, adipose, skeletal muscle,
#              apCAF, PI16 fibroblast, CXCL9 macrophage), each with a module-level reference verified 2026-11-05
# Outputs: registry.json (same schema as V7 + 'source_layer', 'source_locator'), all_panels.tsv, references.tsv,
#          gene_reference_trace_v2_additions.tsv, lodi_subtype_vocabulary.tsv, ALL_PANELS.md
# Usage: build_panel_registry_v2.py V7_REGISTRY LODI_REPO OUT_DIR
import sys, os, json, csv, collections
import pandas as pd

v7p, lodi, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
R = json.load(open(v7p))
LODI_COMMIT = '35f054bb81e1b321718495e096171580b3bc130c'
LODI_URL = 'https://github.com/lambrechtslab/Pan-cancer_Lodietal'

# ---------- references added in v2 (module-level; verified 2026-11-05 via PubMed esummary or Crossref DOI) ----------
NEWREF = [
    ('Lodi2025paper', 'Decoding tumor heterogeneity: A spatially informed pan-cancer analysis of the tumor microenvironment', 'https://doi.org/10.1016/j.xcrm.2025.102416', 'PMID 41086810; PMC12629804. Methods read 2026-11-05: QC, per-cancer major types, pooled subtype clustering, regression, LQ/doublet cluster criteria, abundance thresholds'),
    ('Smillie2019', 'Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis', 'https://doi.org/10.1016/j.cell.2019.06.029', 'PMID 31348891; human colon epithelial atlas (colonocyte, goblet, stem, enteroendocrine). Normal/UC colon, not tumour'),
    ('ZhangP2019', 'Dissecting the Single-Cell Transcriptome Network Underlying Gastric Premalignant Lesions and Early Gastric Cancer', 'https://doi.org/10.1016/j.celrep.2019.04.052', 'PMID 31067475; gastric pit/mucous neck/chief/parietal lineages'),
    ('Baron2016', 'A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure', 'https://doi.org/10.1016/j.cels.2016.08.011', 'PMID 27667365; acinar, ductal, islet endocrine types'),
    ('Kumar2023', 'A spatially resolved single-cell genomic atlas of the adult human breast', 'https://doi.org/10.1038/s41586-023-06252-9', 'PMID 37380767; luminal hormone-sensing / secretory and basal-myoepithelial lineages'),
    ('YuZ2019', 'Single-Cell Transcriptomic Map of the Human and Mouse Bladders', 'https://doi.org/10.1681/ASN.2019040335', 'Crossref-verified DOI; urothelial basal/intermediate/umbrella (UPK) lineages'),
    ('Lake2023', 'An atlas of healthy and injured cell states and niches in the human kidney', 'https://doi.org/10.1038/s41586-023-05769-3', 'PMID 37468583; nephron segments (TAL UMOD/SLC12A1, collecting duct AQP2, podocyte NPHS1/2)'),
    ('Pu2021', 'Single-cell transcriptomic analysis of the tumor ecosystems underlying initiation and progression of papillary thyroid carcinoma', 'https://doi.org/10.1038/s41467-021-26343-3', 'PMID 34663816; thyroid follicular cells (TG/TPO). Also a source dataset in this project (GSE184362)'),
    ('HuZ2020', 'The Repertoire of Serous Ovarian Cancer Non-genetic Heterogeneity Revealed by Single-Cell Sequencing of Normal Fallopian Tube Epithelial Cells', 'https://doi.org/10.1016/j.ccell.2020.01.003', 'PMID 32049047; fallopian secretory (PAX8, OVGP1) vs ciliated epithelium'),
    ('WangW2020', 'Single-cell transcriptomic atlas of the human endometrium during the menstrual cycle', 'https://doi.org/10.1038/s41591-020-1040-z', 'PMID 32929266; endometrial epithelium (PAEP, MSX1, SOX17)'),
    ('Madissoon2019', 'scRNA-seq assessment of the human lung, spleen, and esophagus tissue stability after cold preservation', 'https://doi.org/10.1186/s13059-019-1906-x', 'PMID 31892341; oesophageal squamous basal/suprabasal/superficial (KRT4, KRT13, SPRR, CRNN)'),
    ('ChengJB2018', 'Transcriptional Programming of Normal and Inflamed Human Epidermis at Single-Cell Resolution', 'https://doi.org/10.1016/j.celrep.2018.09.006', 'PMID 30355494; epidermal keratinocyte differentiation (KRT1/KRT10, LOR, FLG)'),
    ('ZhouY2020', 'Single-cell RNA landscape of intratumoral heterogeneity and immunosuppressive microenvironment in advanced osteosarcoma', 'https://doi.org/10.1038/s41467-020-20059-6', 'PMID 33303760; osteoblastic/chondroblastic (RUNX2, SP7, COL2A1, ACAN, IBSP) and osteoclast (ACP5, CTSK, MMP9) in osteosarcoma. Also source GSE152048'),
    ('Emont2022', 'A single-cell atlas of human and mouse white adipose tissue', 'https://doi.org/10.1038/s41586-022-04518-2', 'PMID 35296864; adipocyte (ADIPOQ, PLIN1, LEP)'),
    ('DeMicheli2020', 'A reference single-cell transcriptomic atlas of human skeletal muscle tissue reveals bifurcated muscle stem cell populations', 'https://doi.org/10.1186/s13395-020-00236-3', 'PMID 32624006; myofibre / muscle stem cell markers'),
    ('Elyada2019', 'Cross-Species Single-Cell Analysis of Pancreatic Ductal Adenocarcinoma Reveals Antigen-Presenting Cancer-Associated Fibroblasts', 'https://doi.org/10.1158/2159-8290.CD-19-0094', 'PMID 31197017; apCAF (CD74, HLA-DR), iCAF, myCAF in PDAC'),
    ('Bill2023', 'CXCL9:SPP1 macrophage polarity identifies a network of cellular programs that control human cancers', 'https://doi.org/10.1126/science.ade2292', 'PMID 37535729; CXCL9 vs SPP1 macrophage polarity. Also source GSE234933'),
    ('MacParland2018', 'Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations', 'https://doi.org/10.1038/s41467-018-06318-7', 'Crossref-verified DOI; cited by Lodi for Kupffer cells and sinusoidal EC'),
    ('Cappuyns2023', 'PD-1- CD45RA+ effector-memory CD8 T cells and CXCL10+ macrophages are associated with response to atezolizumab plus bevacizumab in advanced hepatocellular carcinoma', 'https://doi.org/10.1038/s41467-023-43381-1', 'Crossref-verified DOI; cited by Lodi (as Cappuyins et al 2023) for Kupffer cells; identity of the cited paper inferred from author/journal/year'),
    ('Darmanis2017', 'Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma', 'https://doi.org/10.1016/j.celrep.2017.10.030', 'Crossref-verified DOI; cited by Lodi for microglia'),
    ('Masuda2019', 'Spatial and temporal heterogeneity of mouse and human microglia at single-cell resolution', 'https://doi.org/10.1038/s41586-019-0924-x', 'Crossref-verified DOI; cited by Lodi for microglia'),
    ('Muller2017', 'Single-cell profiling of human gliomas reveals macrophage ontogeny as a basis for regional differences in macrophage activation in the tumor microenvironment', 'https://doi.org/10.1186/s13059-017-1362-4', 'Crossref-verified DOI; cited by Lodi for microglia'),
    ('QianCR2020', 'A pan-cancer blueprint of the heterogeneous tumor microenvironment revealed by single-cell profiling', 'https://doi.org/10.1038/s41422-020-0355-0', 'Crossref-verified DOI; same study as V7 Qian2020 (PDF link); cited by Lodi for alveolar EC / macrophages'),
    ('Gillich2020', 'Capillary cell-type specialization in the alveolus', 'https://doi.org/10.1038/s41586-020-2822-7', 'Crossref-verified DOI; cited by Lodi as "Astrid et al, Nature 2020" for alveolar EC (first author Astrid Gillich)'),
    ('ChenCell2021_unresolved', 'Cited by Lodi as "Chen et al, Cell 2021" for alveolar macrophages', 'NA', 'Not resolved to a specific paper in this session; not used as independent evidence'),
    ('McGinnis2019', 'DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors', 'https://doi.org/10.1016/j.cels.2019.03.003', 'Crossref-verified DOI; doublet scoring used by Lodi'),
    ('Korsunsky2019', 'Fast, sensitive and accurate integration of single-cell data with Harmony', 'https://doi.org/10.1038/s41592-019-0619-0', 'Crossref-verified DOI'),
    ('Stuart2019', 'Comprehensive Integration of Single-Cell Data', 'https://doi.org/10.1016/j.cell.2019.05.031', 'Crossref-verified DOI; Seurat CCA integration used by Lodi for B cells'),
]
known = {r['id'] for r in R['references']}
refs = list(R['references']) + [{'id': i, 'title': t, 'url': u, 'verified_scope': s} for i, t, u, s in NEWREF if i not in known]


def panel(id_, axis, genes, parent=(), groups=None, comparators=(), naming=(), context='', refs_=(), level='', interp='',
          layer='', locator='', lower=()):
    return {'id': id_, 'axis': axis, 'parent_any': list(parent), 'core_groups': groups or [list(genes)], 'positive_genes': list(genes),
            'contrast_review_panels': list(comparators), 'expected_lower_genes': list(lower), 'naming_genes': list(naming),
            'required_context': context, 'reference_ids': list(refs_), 'evidence_level': level, 'interpretation': interp,
            'auto_final_label': False, 'group_logic': 'Each core group is inspected separately; no universal AND/OR numeric classifier is claimed.',
            'negative_evidence_policy': 'Missing genes or dropout are not negative; competitors prompt review, never deletion.',
            'source_layer': layer, 'source_locator': locator}


new = []
# ---- Lodi major-type check (exact gene lists, Cancer type analysis.R L328-L336 comments / L344 vector) ----
L = 'Lodi2025;Lodi2025paper'
loc = f'{LODI_URL}/blob/{LODI_COMMIT}/Cancer%20type%20analysis.R#L328-L345'
for k, g, par in (('T', ['CD3D', 'CD3E', 'CD4', 'CD8A'], 'T'), ('B', ['CD79A'], 'B'), ('MacMono', ['CD68', 'LYZ', 'AIF1', 'ITGAM'], 'Myeloid'),
                  ('pDC', ['LILRA4', 'CXCR3'], 'DC'), ('cDC', ['CLEC9A', 'XCR1', 'CD1C', 'CCR7', 'CCL17', 'CCL19'], 'DC'),
                  ('Mast', ['MS4A2', 'TPSAB1', 'CPA3'], 'Mast'), ('EC', ['CLDN5', 'PECAM1', 'VWF'], 'Endothelial'),
                  ('Fibroblast', ['COL1A1', 'BGN', 'DCN'], 'Fibroblast'), ('Epithelial', ['EPCAM', 'KRT7', 'KRT18'], 'Epithelial')):
    new.append(panel(f'Lodi_major_{k}', 'author_major_check', g, parent=[par], groups=[[x] for x in g], refs_=L.split(';'),
                     level='author_code_panel_exact', layer='Lodi2025_code', locator=loc,
                     interp='Lodi canonical major-type marker set, plotted on every major and lineage object to check identity and cross-lineage contamination. Used for review, not as a score threshold.'))
# ---- cancer-type-specific populations identified per cancer type and removed before pan-cancer pooling ----
mm = f'{LODI_URL}/blob/{LODI_COMMIT}/Cell%20subtype%20annotations/Macrophage%3Amonocyte%20analysis.R'
ec = f'{LODI_URL}/blob/{LODI_COMMIT}/Cell%20subtype%20annotations/EC%20analysis.R'
for id_, g, par, ctx, rf, lc in (
        ('Kupffer_like', ['CD163', 'MAFB', 'VSIG4', 'LYVE1', 'CD209', 'PLTP', 'SELENOP', 'MERTK'], 'Macrophage', 'liver (HCC, iCCA)', 'Lodi2025;MacParland2018;Aizarani2019;Cappuyns2023', mm + '#L4'),
        ('Microglia_like', ['TMEM119', 'P2RY12', 'P2RY13', 'SLC2A5', 'CX3CR1', 'FCGR3B'], 'Macrophage', 'brain (glioma)', 'Lodi2025;Darmanis2017;Masuda2019;Muller2017', mm + '#L5'),
        ('Alveolar_macrophage_like', ['FABP4', 'MARCO', 'PPARG', 'INHBA', 'ALDH2'], 'Macrophage', 'lung (NSCLC)', 'Lodi2025;QianCR2020;ChenCell2021_unresolved', mm + '#L6'),
        ('Alveolar_EC_like', ['EDNRB', 'IL1RL1', 'HPGD', 'CLDN18', 'AGER', 'SFTPC'], 'Endothelial', 'lung (NSCLC)', 'Lodi2025;QianCR2020;Gillich2020', ec + '#L4'),
        ('Sinusoidal_EC_like', ['CLEC4G', 'CLEC4M', 'FLT1', 'PECAM1', 'VWF', 'CD34', 'MGP', 'CCL14', 'DNASE1L3', 'SPARCL1', 'CLEC1B', 'LIFR', 'TM4SF1', 'FCN2', 'PTGDS', 'CLEC14A', 'S100A13', 'C7'],
         'Endothelial', 'liver (HCC, iCCA)', 'Lodi2025;Aizarani2019;MacParland2018', ec + '#L5')):
    new.append(panel(id_, 'tissue_specific_population', g, parent=[par], comparators=[par], context=ctx, refs_=rf.split(';'),
                     level='author_code_panel_exact', layer='Lodi2025_code', locator=lc,
                     interp='Lodi: identified within the individual cancer type and removed before pan-cancer subclustering so tissue-restricted cells do not form or distort pan-cancer subtypes. AGER/SFTPC in the alveolar EC list are alveolar context genes (possible ambient/epithelial); review with EC parent evidence.'))
# ---- state signatures (explicit in code) and named external signatures ----
ct = f'{LODI_URL}/blob/{LODI_COMMIT}/Cancer%20type%20analysis.R'
new.append(panel('Lodi_hypoxia_15', 'author_state_signature', ['VEGFA', 'SLC2A1', 'PGAM1', 'ENO1', 'LDHA', 'TPI1', 'P4HA1', 'MRPS17', 'CDKN3', 'ADM', 'NDRG1', 'TUBB6', 'ALDOA', 'MIF', 'ACOT7'],
                 parent=['*'], refs_=['Lodi2025'], level='author_code_panel_exact', layer='Lodi2025_code', locator=ct + '#L217-L218',
                 interp='AddModuleScore hypoxia score; regressed when subclustering macrophages/monocytes and B cells (Methods).'))
new.append(panel('Lodi_epithelial_differentiation_21', 'author_state_signature', ['EPCAM', 'KRT16', 'KRT6A', 'KRT6B', 'KRT6C', 'CLDN4', 'KLK5', 'KLK6', 'KLK7', 'KLK8', 'KLK9', 'KLK10', 'KLK11', 'KRT17', 'KRT75', 'S100A7', 'S100A8', 'S100A9', 'CLDN1', 'CLDN7', 'SPRR1B'],
                 parent=['Epithelial'], refs_=['Lodi2025'], level='author_code_panel_exact', layer='Lodi2025_code', locator=ct + '#L249',
                 interp='AddModuleScore epithelial differentiation score (display/QC feature plot); not an identity panel. S100A8/9 overlap neutrophil/monocyte programs.'))
for id_, desc, lc in (('Lodi_IFN_BROWNE', 'MSigDB c2 BROWNE_INTERFERON_RESPONSIVE_GENES (v7.0 in code, v6.2 in Methods); regressed in T/B/EC/Mac subclustering', '#L211'),
                      ('Lodi_EMT_SARRIO', 'MSigDB c2 SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN; feature plot only', '#L243'),
                      ('Lodi_stress', 'stress.genes.csv (not in the public repository; available from authors on request); regressed in T/B/EC subclustering', '#L232-L235'),
                      ('Lodi_IG', 'all genes matching ^IGH|^IGL|^IGK; IG score feature plot; IGHV/IGKV/IGLV excluded from B-cell HVGs', '#L222-L230')):
    new.append(panel(id_, 'author_state_signature_external', [], parent=['*'], groups=[[]], refs_=['Lodi2025'], level='author_code_named_external_set',
                     layer='Lodi2025_code', locator=ct + lc, interp=desc + '. Gene list not copied; must be obtained from the named source and versioned before use.'))
# ---- Lodi published DEG profiles (marker_genes.txt): top 15 positive genes per cluster by avg_log2FC with pct.1 >= 0.25 ----
mg = pd.read_csv(f'{lodi}/Master_files/marker_genes.txt', sep='\t', index_col=0)
mg = mg[mg['cluster'].notna() & mg['gene'].notna()]
for c, d in mg.groupby('cluster'):
    d = d[(d['avg_log2FC'] > 0) & (d['pct.1'] >= 0.25)].sort_values('avg_log2FC', ascending=False).head(15)
    new.append(panel('LodiDEG_' + str(c).replace(' ', '_').replace('+', 'pos').replace('/', '_'), 'author_DEG_profile', list(d['gene']), parent=[],
                     groups=[list(d['gene'])], refs_=['Lodi2025'], level='author_published_DEG_table', layer='Lodi2025_marker_genes',
                     locator=f'{LODI_URL}/blob/{LODI_COMMIT}/Master_files/marker_genes.txt (cluster={c})',
                     interp=f'Lodi FindAllMarkers output for final subtype "{c}" (top 15 by avg_log2FC, pct.1>=0.25). Reference profile for comparison; Ig/HLA/ambient genes included as published; not an identity classifier.'))
# ---- v1 additions (modules absent from V7) ----
V1 = 'TME_de_novo_v1'
for id_, axis, g, par, ctx, rf, naming, interp in (
        ('Epi_colonocyte_like', 'branch', ['CA1', 'CA2', 'SLC26A3', 'CEACAM7', 'GUCA2A'], 'Epithelial', 'colon/rectum', 'Smillie2019', ['CA1'], 'Absorptive colonocyte direction; tumour cells may retain or lose it.'),
        ('Epi_goblet_like', 'branch', ['MUC2', 'TFF3', 'SPINK4', 'CLCA1', 'FCGBP'], 'Epithelial', 'intestine', 'Smillie2019', ['MUC2'], 'TFF3 also on lymphatic EC; MUC2 required for naming.'),
        ('Epi_intestinal_stem_like', 'branch', ['LGR5', 'OLFM4', 'ASCL2', 'SMOC2'], 'Epithelial', 'intestine', 'Smillie2019', ['LGR5'], 'CRC tumour cells often stem-like; malignancy is a separate axis.'),
        ('Epi_gastric_pit_mucous_like', 'branch', ['MUC5AC', 'TFF1', 'GKN1', 'GKN2', 'MUC6'], 'Epithelial', 'stomach', 'ZhangP2019', ['MUC5AC'], 'Pit/mucous neck direction.'),
        ('Epi_gastric_chief_parietal_like', 'branch', ['PGC', 'PGA3', 'LIPF', 'ATP4A', 'CBLIF'], 'Epithelial', 'stomach', 'ZhangP2019', ['LIPF'], 'CBLIF is the current symbol of GIF (alias must be versioned).'),
        ('Epi_enteroendocrine_like', 'branch', ['CHGA', 'CHGB', 'SYP', 'SCG2', 'PCSK1N'], 'Epithelial', 'GI tract or pancreas', 'Smillie2019;Baron2016', ['CHGA'], 'Also PanNET tumour cells: never names non-malignant without the malignancy axis.'),
        ('Epi_pancreatic_acinar_like', 'branch', ['PRSS1', 'CPA1', 'CTRB1', 'CELA3A', 'PNLIP'], 'Epithelial', 'pancreas', 'Baron2016', ['PRSS1'], 'Strongly ambient-prone; require same-cell support and low ambient fraction.'),
        ('Epi_pancreatic_ductal_like', 'branch', ['CFTR', 'KRT19', 'SPP1', 'MMP7', 'SOX9'], 'Epithelial', 'pancreas', 'Baron2016', ['CFTR'], 'Shares KRT19/SOX9 with cholangiocytes; organ context required.'),
        ('Epi_islet_like', 'branch', ['INS', 'GCG', 'SST', 'PPY'], 'Epithelial', 'pancreas', 'Baron2016', ['INS', 'GCG'], 'Hormone transcripts are ambient-prone.'),
        ('Epi_breast_luminal_like', 'branch', ['KRT18', 'ESR1', 'PGR', 'ANKRD30A', 'AZGP1'], 'Epithelial', 'breast', 'Kumar2023', ['ANKRD30A'], 'Hormone-sensing luminal direction.'),
        ('Epi_breast_basal_myoepithelial_like', 'branch', ['KRT5', 'KRT14', 'TP63', 'ACTA2', 'OXTR'], 'Epithelial', 'breast', 'Kumar2023', ['KRT14'], 'ACTA2 shared with CAF/SMC; epithelial parent required.'),
        ('Epi_prostate_basal_club_like', 'branch', ['KRT5', 'KRT15', 'TP63', 'PIGR', 'LTF'], 'Epithelial', 'prostate', 'Henry2018', ['KRT15'], 'Basal (KRT5/KRT15/TP63) and club (PIGR/LTF) directions; review separately.'),
        ('Epi_urothelial_like', 'branch', ['UPK1B', 'UPK2', 'UPK3A', 'KRT13', 'GATA3'], 'Epithelial', 'bladder / upper urinary tract', 'YuZ2019', ['UPK1B'], 'GATA3 also in T cells and breast luminal cells.'),
        ('Epi_renal_distal_collecting_podocyte_like', 'branch', ['UMOD', 'SLC12A1', 'AQP2', 'NPHS1', 'NPHS2', 'PODXL'], 'Epithelial', 'kidney', 'Lake2023', ['UMOD', 'AQP2', 'NPHS2'], 'Complements V7 Epi_renal_tubular_like (proximal). PODXL also on endothelium.'),
        ('Epi_thyroid_follicular_like', 'branch', ['TG', 'TPO', 'TSHR', 'PAX8'], 'Epithelial', 'thyroid', 'Pu2021', ['TG'], 'PTC tumour cells keep follicular genes; malignancy separate.'),
        ('Epi_mullerian_secretory_like', 'branch', ['PAX8', 'OVGP1', 'MSLN', 'WT1'], 'Epithelial', 'ovary / fallopian tube / endometrium', 'HuZ2020', ['OVGP1'], 'MSLN/WT1 shared with mesothelium.'),
        ('Epi_endometrial_like', 'branch', ['PAEP', 'MSX1', 'SOX17', 'MMP26'], 'Epithelial', 'endometrium', 'WangW2020', ['PAEP'], 'Cycle-phase dependent.'),
        ('Epi_squamous_suprabasal_like', 'branch', ['KRT4', 'KRT13', 'SPRR1A', 'SPRR3', 'CRNN'], 'Epithelial', 'oesophagus / head and neck / cervix', 'Madissoon2019', ['KRT4'], 'Mucosal squamous differentiation; complements V7 Epi_basal.'),
        ('Epi_keratinocyte_like', 'branch', ['KRT1', 'KRT10', 'LOR', 'FLG'], 'Epithelial', 'skin', 'ChengJB2018', ['KRT10'], 'Cornified epidermis direction (BCC/melanoma skin samples).'),
        ('Osteo_chondro_like', 'parent', ['RUNX2', 'SP7', 'ALPL', 'COL2A1', 'ACAN', 'IBSP'], '', 'bone (osteosarcoma)', 'ZhouY2020', ['SP7', 'COL2A1'], 'Malignant osteoblastic/chondroblastic cells share these; normal osteoblast naming needs malignancy-negative evidence.'),
        ('Osteoclast_like', 'branch', ['ACP5', 'CTSK', 'MMP9', 'DCSTAMP'], 'Macrophage', 'bone (osteosarcoma)', 'ZhouY2020', ['DCSTAMP'], 'CTSK/MMP9 also in fibroblasts and macrophages; macrophage parent required.'),
        ('Adipocyte', 'parent', ['ADIPOQ', 'PLIN1', 'LEP', 'CIDEA'], '', '', 'Emont2022', ['ADIPOQ'], 'Rare after dissociation; often ambient.'),
        ('Skeletal_muscle', 'parent', ['ACTA1', 'MYH2', 'CKM', 'DES'], '', '', 'DeMicheli2020', ['ACTA1', 'CKM'], 'DES shared with smooth muscle; Lodi lists muscle cells as a cancer-specific major type.'),
        ('Fib_apCAF_like', 'state', ['CD74', 'HLA-DRA', 'HLA-DRB1', 'SLPI'], 'Fibroblast', 'tumour', 'Elyada2019', ['HLA-DRA'], 'MHC-II fibroblast state; rule out myeloid doublets first.'),
        ('Fib_PI16_universal_like', 'state', ['PI16', 'DPT', 'CD34'], 'Fibroblast', '', 'Buechler2021', ['PI16'], 'Universal adventitial fibroblast direction.'),
        ('Mac_CXCL9_like', 'state', ['CXCL9', 'CXCL10', 'GBP1', 'STAT1'], 'Macrophage', 'tumour', 'Bill2023', ['CXCL9'], 'CXCL9:SPP1 polarity axis; overlaps IFN_response (Lodi "IFN Mac").')):
    new.append(panel(id_, axis, g, parent=[par] if par else [], comparators=['Mesothelial'] if par == 'Epithelial' else [], naming=naming, context=ctx,
                     refs_=rf.split(';'), level='curated_literature_informed_module_reference_verified', layer=V1,
                     locator='methods/annotation_gene_panel_v1.tsv (TME de novo v1, 2026-11-05)', interp=interp))
ids = [p['id'] for p in R['panels']]
dup = [p['id'] for p in new if p['id'] in ids]
assert not dup, dup
for p in R['panels']:
    p.setdefault('source_layer', 'V7_panel_repair_v1'); p.setdefault('source_locator', 'resources/annotation/v7_panel_repair_v1/registry.json')
R2 = dict(R)
R2['resource_version'] = 'tme_annotation_panel_v2'
R2['status'] = 'evidence_registry_not_validated_cell_classifier'
R2['derived_from'] = {'V7_panel_repair_v1': 'unchanged 112 modules', 'Lodi2025_code': f'{LODI_URL}@{LODI_COMMIT}', 'TME_de_novo_v1': 'methods/annotation_gene_panel_v1.tsv'}
R2['references'] = refs
R2['panels'] = R['panels'] + new
R2['aliases'] = dict(R['aliases'], GIF='CBLIF')
json.dump(R2, open(f'{out}/registry.json', 'w'), ensure_ascii=False, indent=1)
rows = []
for p in R2['panels']:
    rows.append({'id': p['id'], 'axis': p['axis'], 'source_layer': p['source_layer'], 'parent_context': ';'.join(p['parent_any']) or 'NA',
                 'core_components': ' / '.join(' + '.join(g) for g in p['core_groups']) or 'NA', 'positive_genes': ';'.join(p['positive_genes']) or 'NA',
                 'expected_lower_genes': ';'.join(p['expected_lower_genes']) or 'NA', 'comparators': ';'.join(p['contrast_review_panels']) or 'NA',
                 'naming_genes': ';'.join(p['naming_genes']) or 'NA', 'required_context': p['required_context'] or 'NA', 'evidence_level': p['evidence_level'],
                 'interpretation': p['interpretation'], 'references': ';'.join(p['reference_ids']), 'source_locator': p['source_locator']})
T = pd.DataFrame(rows)
T.to_csv(f'{out}/all_panels.tsv', sep='\t', index=False)
pd.DataFrame(refs).to_csv(f'{out}/references.tsv', sep='\t', index=False)
tr = []
for p in new:
    for g in p['positive_genes']:
        st = {'author_code_panel_exact': 'author_code_exact_gene_list', 'author_published_DEG_table': 'author_DEG_table_exact'}.get(
            p['evidence_level'], 'module_reference_verified_gene_level_pending')
        tr.append({'panel': p['id'], 'gene': g, 'direction': 'positive', 'source_layer': p['source_layer'], 'source_locator': p['source_locator'],
                   'references': ';'.join(p['reference_ids']), 'status': st})
pd.DataFrame(tr).to_csv(f'{out}/gene_reference_trace_v2_additions.tsv', sep='\t', index=False)
voc = pd.read_csv(f'{lodi}/Master_files/Minor_2.0_Intermediatecelltype_annotation.csv').drop_duplicates()
voc['source_locator'] = f'{LODI_URL}/blob/{LODI_COMMIT}/Master_files/Minor_2.0_Intermediatecelltype_annotation.csv'
voc.to_csv(f'{out}/lodi_subtype_vocabulary.tsv', sep='\t', index=False)
with open(f'{out}/ALL_PANELS.md', 'w') as f:
    f.write('# TME annotation panel registry v2\n\nGenerated by `scripts/annotation_v2/build_panel_registry_v2.py`. Machine-readable: [registry.json](registry.json). '
            'Evidence registry, not a cell classifier.\n\n')
    cnt = T.groupby(['source_layer', 'axis']).size().reset_index(name='modules')
    f.write('| source layer | axis | modules |\n|---|---|---|\n' + ''.join(f'| {r.source_layer} | {r.axis} | {r.modules} |\n' for r in cnt.itertuples()) + '\n')
    for layer, d in T.groupby('source_layer', sort=False):
        f.write(f'## {layer}\n\n| id | axis | parent | positive genes | naming | context | evidence | references |\n|---|---|---|---|---|---|---|---|\n')
        for r in d.itertuples():
            f.write(f'| {r.id} | {r.axis} | {r.parent_context} | {r.positive_genes.replace(";", ", ")} | {r.naming_genes} | {r.required_context} | {r.evidence_level} | {r.references} |\n')
        f.write('\n')
print(json.dumps({'panels_total': len(R2['panels']), 'v7': len(R['panels']), 'added': len(new),
                  'added_by_layer': dict(collections.Counter(p['source_layer'] for p in new)), 'references': len(refs),
                  'genes_total': len({g for p in R2['panels'] for g in p['positive_genes']})}, indent=1))
