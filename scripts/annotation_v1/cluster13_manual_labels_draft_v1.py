import sys, pandas as pd
D, O = sys.argv[1:3]
d = pd.read_csv(D + '/subclusters.tsv', sep='\t').set_index('sub_r1')
# sub: label, class, confidence, residual call, evidence
L = {
10:('breast_luminal_progenitor_SAA1_LTF','organ_epithelium','medium','not_flagged','SAA1,CCL28,LTF,MMP7,SLPI,GABRP,SFRP1,KRT15; TP63 1%, ACTA2 15% argue against basal-myoepithelial; auto prostate_basal_club wrong; 98% breast'),
35:('prostate_basal_KRT15_SLC14A1','organ_epithelium','medium','not_flagged','KRT15,SLC14A1,AQP3,TRIM29,KRT5 72%,TP63 42%; 90% prostate; top source 0.855 so source_restricted; auto endo_capillary wrong'),
7:('prostate_luminal_club_like_WFDC2_TMPRSS2','organ_epithelium','low','identity_with_flag','PPP1R1B,TMPRSS2,WFDC2,RARRES1,ZG16B; CNV corr 0.21 but scATOMIC tumour 0.58, uncertain 0.74; auto Th17 wrong'),
20:('alveolar_AT2_SFTPB_NAPSA','organ_epithelium','high','not_flagged','SFTPB,NAPSA,SFTA2,SFTPD,SFTPA1; 6/7 sources; 99% NSCLC'),
19:('colonic_goblet_MUC2_TFF3','organ_epithelium','high','identity_with_flag','MUC2,TFF3,FCGBP,SPINK4; 8/8 sources; flag driven by uncertain 0.63 and CNV corr 0.31'),
21:('airway_club_SCGB3A1_BPIFB1','organ_epithelium','medium','identity_with_flag','SCGB3A1,BPIFB1,WFDC2,SLPI; 8/8 sources; uncertain 0.80'),
43:('breast_luminal_HR_ANKRD30A_ESR1','tumour_like','medium','likely_residual_tumour','ANKRD30A,ESR1,GATA3,TRPS1; scATOMIC tumour 0.86, uncertain 0.89 although CNV corr 0.19'),
58:('pancreatic_duct_reactive_MMP7_KRT7','organ_epithelium','low','not_flagged','MMP7,KRT7,SPINK1,MUC1,HLA-DR; source share 0.956; uncertain 0.30, CNV corr 0.50'),
59:('pancreatic_duct_CFTR_SLC4A4','organ_epithelium','medium','identity_with_flag','CFTR,SLC4A4,ANXA4,FXYD2,DCDC2 ahead of PRSS1/PNLIP (acinar transcripts read as ambient, interpretation); auto acinar not supported; source share 0.993'),
9:('prostate_luminal_KLK3_MSMB','organ_epithelium','medium','identity_with_flag','KLK3,KLK2,MSMB,ACPP; MT 22% but specific markers so not low_quality; source share 0.968'),
29:('gastric_foveolar_MUC5AC_TFF1','organ_epithelium','medium','identity_with_flag','MUC5AC,TFF1,CTSE,VSIG1; auto chief_parietal not supported'),
68:('mesothelial_ITLN1_C1S','non_epithelial_parenchyma','medium','not_flagged','ITLN1,C1S,C1R,EFEMP1; MSLN 36%,UPK3B 31%; uncertain 0.02, CNV corr 0.13; 93% mesothelioma, source share 0.93; auto CD8_tissue_resident wrong'),
46:('retinal_pigment_epithelium_BEST1_RDH5','non_epithelial_parenchyma','medium','not_flagged','BEST1,RDH5,LRAT,TTR,RBP1,TYRP1,SERPINF1; 89% uveal melanoma; 36 cells; auto islet wrong'),
73:('cholangiocyte_like_DCDC2_FGFR2','organ_epithelium','low','not_flagged','DCDC2,BICC1,FGFR2,ANXA4,SLC4A4,CA2,SPP1; 30 cells, 1 source; auto hepatocyte not supported'),
61:('melanocytic_MLANA_MITF','tumour_like','medium','likely_residual_tumour','MLANA,MITF; 79% uveal melanoma; uncertain 0.87, CNV corr 0.70'),
42:('low_quality_42','low_quality','medium','not_flagged','median genes 398 (<500), MT 19%, contamination 0.16; source share 0.926'),
38:('unresolved_38_low_information_nuclear','ambient_or_doublet','low','not_flagged','median genes 585; long nuclear transcripts (ZSWIM6,QKI,RBMS3,DMD) with fibroblast and endothelial panels; auto endo_venous not supported'),
44:('unresolved_44_multilineage_ambient','ambient_or_doublet','medium','not_flagged','COL1A1 92%,COL1A2 99%,HBB 99%,EPCAM 78%,PTPRC 40%; 1 source, 1 unit; auto fibroblast not supported'),
39:('erythroid_contaminated_epithelial','ambient_or_doublet','medium','technical_not_a_tumour_call','erythroid z 6.65 above any L1; S100A9/CSTA epithelial background; 80% NSCLC'),
18:('squamous_KRT5_KRT6A_stromal_admixture','tumour_like','low','likely_residual_tumour','KRT5 91%,KRT6A 91%,COL1A1 56%,RGS5 23%; CSF3R 1%,FCGR3B 0%: no neutrophil evidence, auto doublet not supported; stromal doublet possible'),
65:('doublet_like_squamous_fibroblast','ambient_or_doublet','low','likely_residual_tumour','KRT14 90%,COL1A1 89%,COL3A1,SPARC; co-detection at cell level not measured, partial EMT alternative'),
33:('squamous_KRT6A_S100A8_S100A9','tumour_like','medium','likely_residual_tumour','KRT6A 74%,KRT5 69%,S100A8 80%,S100A9 84%,PI3,GJB2; CSF3R 1%,FCGR3B 0%,CXCR2 0%: no neutrophil evidence, auto doublet not supported; 53% HNSCC, CNV corr 0.58'),
31:('squamous_partial_EMT_FN1_INHBA','tumour_like','low','likely_residual_tumour','FN1,INHBA,TNC,LAMC2,COL17A1,ITGA6; 42% HNSCC; CNV corr 0.55'),
52:('basal_like_KRT17_FHL2','tumour_like','low','likely_residual_tumour','KRT17,FHL2,TNFRSF12A,TAGLN,CAV1; 73% BCC; scATOMIC tumour 0.80; auto myCAF not supported'),
37:('basal_like_KRT14_AREG','tumour_like','low','likely_residual_tumour','KRT14,KRT17,SFN,AREG,CXCL8; 83% BCC; scATOMIC tumour 0.63'),
54:('thyroid_follicular_FN1_HOPX','tumour_like','medium','likely_residual_tumour','TG,PAX8,FN1,HOPX; FN1 reads as PTC-associated (interpretation); CNV corr 0.65'),
69:('thyroid_like_CHI3L1_TIMP1','tumour_like','low','likely_residual_tumour','99% PTC; scATOMIC tumour 0.95; auto gamma_delta wrong'),
53:('melanoma_like_VIM_S100A1','tumour_like','low','likely_residual_tumour','85% cutaneous melanoma; uncertain 0.93; scATOMIC tumour 0.73; auto CD4_naive_CM wrong'),
27:('breast_tumour_like_PEG10_MDK','tumour_like','low','likely_residual_tumour','CNV score 0.0084, CNV corr 0.74; source share 0.976; auto PI16 fibroblast wrong'),
}
prog = {
23:('program_stress_IEG_low_depth','JUN,FOS,HSPA1A,FOSB,DNAJB1; median genes 852; auto plasma_IgG wrong'),
16:('program_stress_IEG','JUN,EGR1,FOS,DNAJB1,ATF3; auto enteroendocrine wrong'),
24:('program_low_depth_ambient','NEAT1 with discordant MLANA,ALB,PMEL; median genes 626; auto plasma_IgG wrong'),
4:('program_ribosomal_glycolytic','GAPDH,RPLP0,MIF,FTL; auto Kupffer wrong'),
3:('program_housekeeping_MHC','GAPDH,FTL,B2M,HLA-A,LGALS1; auto endometrial wrong'),
2:('program_OXPHOS_glycolytic','GAPDH,COX5B,TPI1,UQCRH; auto endometrial wrong'),
17:('program_MHC_interferon','HLA-A/B/C,CD74,HLA-DRA,IFI27,IFI6; auto macrophage_IL1B wrong'),
0:('program_translation_ribosomal','RPLP0,NPM1,EEF1B2,EEF2; auto enteroendocrine wrong'),
28:('program_generic_epithelial_LGALS3_CLDN4','LGALS3,ELF3,CLDN4,S100A6; 33% CRC; organ identity not established; auto colonocyte not supported'),
13:('program_proliferation','TOP2A,UBE2C,BIRC5,MKI67 panel; cycling z 6.06; 34/34 sources; auto endometrial wrong'),
30:('program_low_information_MHC_lncRNA','HLA-A/B,CCSER1,NEAT1,MAPK10; MT 11.5%; auto Treg wrong'),
14:('program_hypoxia_glycolysis','NDRG1,PGK1,SLC2A1,BNIP3; auto endo_tip wrong'),
}
small = {
57:('urothelial_tumour_like_AQP3_S100P','tumour_like'),40:('breast_luminal_tumour_like_GATA3_TRPS1','tumour_like'),
60:('breast_tumour_like_SFRP1_MGP','tumour_like'),62:('breast_tumour_like_KRT15_ribosomal','tumour_like'),
56:('breast_tumour_like_CRABP2','tumour_like'),66:('breast_tumour_like_CRABP2_PRRX1','tumour_like'),11:('breast_tumour_like_CTNNB1_SOX4','tumour_like'),
22:('breast_tumour_like_KRT81_GABRP','tumour_like'),25:('hepatocyte_tumour_like_APOA1_ALB','tumour_like'),63:('hepatocyte_tumour_like_SERPINA1_ALB','tumour_like'),
6:('hepatic_intestinal_tumour_like_SPINK1_CES2','tumour_like'),1:('thyroid_tumour_like_TG','tumour_like'),50:('basal_like_BCC_KRT17','tumour_like'),
36:('squamous_tumour_like_KRT6A_MYC','tumour_like'),47:('squamous_suprabasal_tumour_like_KRT13_SPRR3','tumour_like'),
55:('colon_epithelial_tumour_like_MUC2_LGALS4','tumour_like'),5:('gastrointestinal_epithelial_tumour_like_LGALS4_TFF1','tumour_like'),
67:('colon_epithelial_tumour_like_TFF3_CLDN3','tumour_like'),70:('gastric_epithelial_tumour_like_EPCAM','tumour_like'),
64:('pancreatic_tumour_like_MSLN','tumour_like'),72:('pancreatic_tumour_like_AKR1C2_S100P','tumour_like'),
34:('renal_pancreatic_tumour_like_PDZK1IP1','tumour_like'),71:('renal_tumour_like_NDUFA4L2','tumour_like'),45:('renal_tumour_like_SPP1_VCAN','tumour_like'),
12:('bladder_tumour_like_S100P_CCND1','tumour_like'),15:('prostate_tumour_like_FOLH1_FASN','tumour_like'),
41:('prostate_tumour_like_IEG','tumour_like'),8:('epithelial_tumour_like_ISG15_EPCAM','tumour_like'),
49:('program_housekeeping_PDAC','program_mixed'),48:('unresolved_48_15_cells','ambient_or_doublet'),
26:('low_quality_26_ambient','low_quality'),32:('low_quality_32','low_quality'),51:('low_quality_51','low_quality'),
}
rows = []
for s, r in d.iterrows():
    if s in L: lab, cl, conf, res, ev = L[s]
    elif s in prog: lab, ev = prog[s]; cl, conf, res = 'program_mixed', 'low', 'program_no_identity'
    elif s in small:
        lab, cl = small[s]; conf = 'low'
        res = {'tumour_like':'likely_residual_tumour','program_mixed':'program_no_identity'}.get(cl, 'technical_not_a_tumour_call' if r.residual_tumour_suspect else 'not_flagged')
        ev = f'top DE {",".join(r.top_de.split(",")[:6])}; {r.non_malignant} non-malignant cells; uncertain {r.uncertain_share}, CNV corr {r.cnv_corr_median}'
    else: sys.exit(f'unlabelled {s}')
    rows.append(dict(sub_r1=s, cells=r.cells, non_malignant=r.non_malignant, uncertain_share=r.uncertain_share, sources=r.sources, top_source_share=r.top_source_share, cnv_corr_median=r.cnv_corr_median, scatomic_tumour_share=r.scatomic_tumour_share, auto_label=r.auto_label, source_restricted=r.source_restricted, residual_tumour_suspect_auto=r.residual_tumour_suspect, manual_label_draft=lab, manual_class=cl, manual_confidence=conf, residual_review=res, evidence=ev))
m = pd.DataFrame(rows).sort_values('non_malignant', ascending=False)
assert len(m) == 74 and m.non_malignant.sum() == 68374 and m.cells.sum() == 290546
m.to_csv(O + '/manual_labels_draft.tsv', sep='\t', index=False)
f = m[m.residual_tumour_suspect_auto]
f[['sub_r1','cells','non_malignant','uncertain_share','cnv_corr_median','scatomic_tumour_share','top_source_share','manual_label_draft','manual_class','residual_review','evidence']].to_csv(O + '/residual_tumour_suspects_for_user.tsv', sep='\t', index=False)
print(m.groupby('manual_class')[['cells','non_malignant']].sum().assign(n=m.groupby('manual_class').size()).sort_values('non_malignant', ascending=False))
print(m.groupby('residual_review')[['cells','non_malignant']].sum().assign(n=m.groupby('residual_review').size()))
print('flag', f.non_malignant.sum(), len(f), 'SR', int(m.source_restricted.sum()))
print(pd.crosstab(m.residual_tumour_suspect_auto, m.residual_review))
print(m[m.manual_class=='organ_epithelium'].non_malignant.sum())
