# Bulk readout v2 and longitudinal single-cell check — pre-specification
Written 2026-11-26, before any v2 bulk outcome or single-cell pre/post comparison was computed.
Motivation: clean CoVarNet rerun gave a Tex-centred module with two independent CD4 partners (Treg, CXCL13+ Tph).
The v1 bulk readout (pCM02: regulatory FOXP3/IL2RA/CTLA4/IKZF2/TNFRSF18 x generic cytotoxic GZMB/PRF1/CTSW/NKG7/CD8A/CD8B)
does not match these states. v1 results are kept unchanged; v2 is a new, separately reported family.

## 1. Arm genes (fixed rule, single-cell only, no outcome data)
Candidates:
- Treg: FOXP3, IL2RA, IKZF2, TNFRSF18, TNFRSF4, CCR8, IL1R2, CTLA4, BATF, LAYN
- Tph: CXCL13, CXCR5, MAF, TOX2, PDCD1, ICOS, BTLA, IL21, CD200, SH2D1A
- Tex: HAVCR2, LAG3, ENTPD1, LAYN, PDCD1, TOX, CXCL13, TNFRSF9, KIR2DL4, PHLDA1, CTLA4, ITGAE, CD8A, CD8B
Rule (computed on curated T states + all non-T majors of the 4 whole-tumour CM02 sources, immune atlas v1):
gene kept for arm X if (a) detected in >=25% of cells of state X; (b) mean log-normalised expression in state X
>= 2x that of each of the other two arm states; (c) >= 3x the mean of every non-T major lineage present
(Myeloid, B, Plasma, NK, Mast, Epithelial, Endothelial, Fibroblast, Pericyte_SMC). Arms with <3 kept genes are reported unavailable.
Pan-T control: CD3D, CD3E, TRAC (unchanged). Scoring unchanged from v1: within-cohort per-gene percentile, arm mean.

## 2. Readouts
- R1 Treg–Tex joint = min(Treg arm, Tex arm); R2 Tph–Tex joint = min(Tph arm, Tex arm). Per SD, adjusted for pan-T.
- Single-arm models (each arm + pan-T) reported descriptively.
- Context benchmarks (descriptive): Ayers 18-gene T-cell-inflamed GEP, CYT (GZMA, PRF1), and v1 pCM02 recomputed.
- Increment: likelihood-ratio of R1 over (pan-T + GEP) — descriptive.

## 3. Cohorts and endpoints
Pre-treatment samples only; one sample per patient (mean if repeated).
- Held-out validation (never used in this project): Liu 2019 and Gide 2019 melanoma anti-PD-1 (cBioPortal iAtlas harmonised).
- Previously analysed (re-analysis, not validation): IMvigor210, SU2C-MARK, Riaz 2017, Hugo 2016.
Endpoints: objective response CR/PR vs SD/PD (primary); OS (secondary); PFS where available.
## 4. Primary hypothesis and criterion
H1: R1 (Treg–Tex) per SD, adjusted for pan-T, is associated with response in the held-out cohorts:
inverse-variance fixed-effect pooled log-OR over Liu+Gide with 95% CI excluding 0 (direction not pre-assumed; report as is).
H2 (secondary): same for R2 (Tph–Tex).
Multiplicity: BH over all R1/R2 x {response, OS, PFS} x cohort estimable slots (non-estimable slots counted as P=1).

## 5. Single-cell longitudinal check (curated states, cells_phaseB_states)
- GSE286827 (HNSCC, neoadjuvant ICI): paired pre vs post per patient — Treg/CD4, Tph/CD4, Tex/CD8 (Wilcoxon signed-rank);
  edge strength pre vs post: Spearman Treg~Tex and Tph~Tex in each phase, difference with patient bootstrap (2000) CI.
- GSE169246 (TNBC, paclitaxel +/- atezolizumab): same, Pre vs Post tumour, split by treatment arm; by author response where available (descriptive).
No new state definitions; any result is exploratory.

## Amendment A (2026-11-26, after gene selection and clinical-field inspection; before any score–outcome or pre/on computation)
1. Arm selection result: Treg = FOXP3, IL2RA, IKZF2, CCR8, IL1R2; Tex = LAG3, CD8A, CD8B; Tph arm empty (shared CXCL13/PDCD1/ICOS with Tex,
   CXCR5 in B cells) -> R2 unavailable in bulk, reported as such. First selection run wrongly compared against all non-T lineages
   (incl. 1–10-cell Glial/Neuron groups); corrected to the 9 listed lineages; both outputs kept.
   Known limitation: Tex arm is dominated by CD8A/CD8B (rule did not require specificity against other CD8 states).
2. Added (exploratory, not part of H1/H2): bulk pre vs on-treatment change in Riaz 2017 and Gide 2019 paired patients:
   Treg arm, Tex arm, R1, pan-T; Wilcoxon signed-rank; stratified by response (CR/PR vs SD/PD), descriptive.
3. Descriptive sensitivity: Liu 2019 + prior anti-CTLA4 indicator; Gide 2019 + combination (ipilimumab) indicator.
