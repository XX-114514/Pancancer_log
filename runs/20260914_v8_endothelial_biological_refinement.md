# Run: 20260914_v8_endothelial_biological_refinement

## TASK

Which Endothelial identities and states are supported strongly enough to enter V8?

## INPUT

Frozen V7 Endothelial candidate universe only; preregistered fixed rules; full-library log1p(CP10K) expression; four discovery samples from GSE222315/GSE256136 and two independent-validation samples from GSE274229; genuine scATOMIC predictions for S37.

## EXECUTION

`20260914_endothelial_3gse_v1` ran per-original-sample Leiden clustering and fixed identity/state/contamination gates. All 6/6 samples completed between 09:00:43 and 09:01:21 +08:00. State recurrence was evaluated only among identity-positive clusters. Exact source-cell IDs were checked for all 8,951 cells; the source objects contain no biological annotation fields. Existing genuinely rerun scATOMIC predictions were exactly mapped to all 1,547 S37 target cells. No V7 or global V8 writeback occurred.

## RESULT

The run produced 35 clusters/8,951 cells. Endothelial identity passed in 7,881 cells/25 clusters/3 GSE/6 samples; 1,070 cells/10 clusters triggered the fixed competing-lineage contamination fallback. In S37, scATOMIC called 1,538/1,547 target cells `Endothelial Cells`, 2 CAF/EndMT-like and 7 broad/ambiguous, all confident.

Arterial raw evidence recurred in 3 GSE/5 samples/1,978 cells with 382 validation cells and marker coherence in 5/5 positive clusters; one 398-cell capillary-arterial boundary cluster was excluded to parent, leaving 1,580 frozen arterial assignments. Venous recurred in 3 GSE/5 samples/1,718 cells with 1,408 validation cells and marker coherence in 6/6 clusters. Capillary technically recurred but all 6/6 clusters overlapped angiogenic-tip and only 3/6 had top-10 CA4/RGCC coherence. Angiogenic-tip was broad/overlapping; lymphatic had only 35 cells in one sample.

## DECISION

`FROZEN_V8_ENDOTHELIAL_TAXONOMY`: freeze Endothelial parent, `Endothelial_arterial`, and `Endothelial_venous`. Final review assignments are 4,583 parent Endothelial, 1,580 arterial, 1,718 venous and 1,070 unresolved. Capillary, angiogenic-tip and lymphatic remain parent/overlay evidence. This is a lineage freeze only; global V8 remains NOT_FROZEN and V7 remains default.

## NEXT

Run lightweight B/Plasma validation with the same fixed lineage template.

External evidence: `${PROJECT_ROOT}/Project_v3/v8_upgrade/biological_freeze_sprint/endothelial_20260914/FROZEN_V8_ENDOTHELIAL_TAXONOMY`.
