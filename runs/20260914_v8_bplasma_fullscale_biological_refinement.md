# Run: 20260914_v8_bplasma_fullscale_biological_refinement

## TASK

Which B/Plasma identities and states are supported strongly enough by full-cohort development evidence to enter V8?

## INPUT

Frozen V7 `Bcell`/`Plasma`/`B_Plasma_boundary` universe: 693,261 cells, 1,312 original samples and 54 GSE; frozen 359-gene full-library `log1p(CP10K)` expression; preregistered identity/state/contamination gates; prospectively sealed 12-GSE held-out split.

## EXECUTION

The fixed runner processed every original sample with at least 50 scoped cells using per-sample PCA, neighbors, Leiden clustering and marker-panel evidence. Full run `20260914_bplasma_full_v1` completed from 22:10:15 to 22:18:17 +08:00. Taxonomy review used only `heldout_gse=False` outputs. Genuine existing scATOMIC predictions were exactly mapped for two non-heldout GSE274229 samples; a separate exact source audit was run on two GSE188711 samples. No V7/global-V8 writeback occurred.

## RESULT

Full execution reached 981 success and 331 explicit `<50-cell` skips, 0 failed, 685,843 processed cells and 4,378 clusters. Development-only review contained 783 success samples, 39 GSE, 557,601 unique cells and 3,501 clusters; 11/11 implementation contracts passed.

Bcell passed in 667 clusters/213,675 cells/33 GSE/377 samples with canonical top-10 coherence 604/667. Plasma passed in 576/85,574/27/278 with coherence 574/576. B-naive passed in 364 clusters/92,435 cells/28 GSE/252 samples with coherence 313/364. Plasma-IgA covered 540/576 Plasma clusters (93.75%), triggering the preregistered >=90% universal-state parent fallback. Competing-lineage fallback retained 252,944 development cells.

scATOMIC mapped exactly 1,015/1,015 cells: 47 B support, 148 broad and 820 competing; Plasma was not evaluable because the model has no Plasma category. GSE188711 mapped exactly 2,338/2,338 source cells but exposed no cell-level biological annotation, so source concordance was recorded as not assessable rather than inferred.

## DECISION

`FROZEN_V8_B_PLASMA_TAXONOMY`: freeze Bcell, Plasma and `Bcell_naive`. Plasma-IgA falls back to Plasma; IFN/cycling remain overlays; boundary/competing calls remain unresolved. Final development downstream states are Bcell-naive 92,435; Bcell parent 121,268; Plasma parent 85,574; B/Plasma unresolved 5,380; competing-lineage unresolved 252,944. This is a lineage taxonomy freeze only; held-out predictions remain sealed, global V8 remains NOT_FROZEN and V7 remains default.

## NEXT

Run the preregistered 12-GSE held-out reference benchmark without changing frozen rules.

External evidence: `${PROJECT_ROOT}/Project_v3/v8_upgrade/biological_freeze_sprint/b_plasma_20260914/FROZEN_V8_B_PLASMA_TAXONOMY`.
