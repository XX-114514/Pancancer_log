# V8 Myeloid/DC third-GSE biological validation

## TASK

Does a rule-frozen independent third GSE plus a genuinely rerun classifier support promotion of the proposed Myeloid/DC taxonomy?

## INPUT

Frozen V7 candidate universe; GSE161529/GSE131907 reviewed evidence; independent prostate GSE274229 samples `GSM8445684_S5_18_GEX` and `GSM8445716_S37_32_GEX`; full-gene counts for scATOMIC.

## EXECUTION

Rules were frozen before third-GSE unblinding (SHA-256 `53db08b20aa70a20414b36b0165ffc673f4bb943c248b93e80c22affc758c649`). The targeted Scanpy review completed 2/2 samples. scATOMIC 2.0.3 attempt 1 failed sparse conversion; the single permitted retry predicted 17,109/17,109 cells. Exact cell-ID mapping covered all 2,117 target cells.

## RESULT

Third GSE produced 13 clusters: C1QC macrophage 1,488 cells/7 clusters/2 samples; cDC2 95 cells/1 cluster; parent fallback 492 cells; unresolved 42 cells. No SPP1/TREM2 or pDC cluster passed. scATOMIC supported 90/95 cDC2, 994/1,488 C1QC-marker macrophage and 443/492 parent-fallback cells.

## DECISION

`FROZEN_V8_MYELOID_TAXONOMY`: accept FCN1 monocyte, Macrophage and cDC2 identities; accept C1QC and inflammatory states; SPP1/TREM2 and pDC fall back to parents. This is a lineage-taxonomy freeze, not global V8 cell-level freeze or held-out accuracy.

## NEXT

Run T/NK refinement with the same frozen workflow template.

External evidence: `${PROJECT_ROOT}/Project_v3/v8_upgrade/biological_freeze_sprint/myeloid_dc_20260913/FROZEN_V8_MYELOID_TAXONOMY`.
