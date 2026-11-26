# Phase B pre-specification (written 2026-11-26 before any phase-A cluster or association output was viewed)

Question: after curated T-state annotation, do Treg and tumour-reactive exhausted CD8 (Tex) co-vary across tumours, and is this specific?

## States (assigned at cluster level from markers + CellTypist, reviewed by hand)
CD4: Treg (incl. cycling Treg), CD4 naive, CD4 Tcm/memory, CD4 Tfh/CXCL13, CD4 Th17/other effector.
CD8: Tex (CXCL13/HAVCR2/PDCD1/LAYN; incl. cycling Tex), CD8 GZMK Tem, CD8 TRM (ZNF683/ITGAE w/o Tex), CD8 Temra (FGFBP2/CX3CR1), CD8 naive.
Other: MAIT, gdT, NK-like, ambiguous/doublet-like. Stress/cycling/ISG are state flags, not identities.

## Primary measures (per analysis unit)
- Treg_CD4 = Treg / all CD4 T; Tex_CD8 = Tex / all CD8 T (lineage-internal denominators; avoids CD4:CD8 closure).
- Units with <50 CD4 or <50 CD8 cells excluded. Patient-level mean where verified patient keys exist.

## Primary test
Spearman(Treg_CD4, Tex_CD8) within each source, tumour tissue, baseline/pre-treatment only:
discovery GSE286827 pre, GSE294109, GSE305078, (GSE256326 descriptive, n=5);
replication GSE169246 tumour-site Pre (non-blood), GSE176021 tumour (post-neoadjuvant; noted).
Support = positive rho in every evaluable source with n>=10, and Fisher-z pooled replication-only estimate >0 with 95% CI excluding 0.

## Specificity / controls (descriptive, reported in full)
- Treg_CD4 vs each other CD8 state fraction (GZMK Tem, TRM, Temra, naive): is Tex the strongest partner?
- Tfh/CXCL13 CD4 vs Tex (expected positive, tumour-reactive positive control).
- Same correlation in blood and adjacent normal/LN (expected weak; tumour context control).
- Partial for stress-flag fraction and per-unit T count; CLR within-T version.
No result here is used to change state definitions after viewing.
