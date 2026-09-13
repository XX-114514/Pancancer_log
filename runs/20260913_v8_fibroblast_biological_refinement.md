# Run: 20260913_v8_fibroblast_biological_refinement

## TASK

Which Fibroblast identities and states are supported strongly enough to enter V8?

## INPUT

Frozen V7 Fibroblast/stromal candidate universe; preregistered fixed rules; original-sample full-library log1p(CP10K) marker expression; four discovery samples from GSE131907/GSE188737 and two independent-validation samples from GSE184880; source-provided labels.

## EXECUTION

20260913_fibroblast_3gse_v1 ran sample-specific clustering and fixed identity/state gates. All 6/6 samples completed. The 10,147 review cells were then exact-unique mapped to source labels, and state recurrence was recomputed only within identity-gate-positive clusters. No V7 writeback occurred.

## RESULT

The run produced 30 clusters. Fibroblast identity passed for 7,620 cells in 16 clusters across 3 GSE/6 samples. Source major labels supported Fibroblast/stromal in 10,060/10,147 cells overall and in 7,619/7,620 accepted identity cells. myCAF was the only discriminative state accepted: 2,748 cells/10 clusters/3 GSE/4 samples, including 1,055 validation cells.

The remaining states fall back: iCAF had zero validation cells; ECM-remodeling passed all 16/16 identity clusters and was not discriminative from identity; antigen-presentation-like was restricted to one GSE; IFN was two-GSE/two-sample and co-occurred with myCAF; cycling had zero accepted clusters. The 2,527 contamination-gate cells remain unresolved_stromal; source support does not override the frozen gate.

## DECISION

FROZEN_V8_FIBROBLAST_TAXONOMY: freeze Fibroblast parent and Fibroblast_myCAF as the only Fibroblast-specific downstream state. Final review assignments are 2,748 Fibroblast_myCAF, 4,872 Fibroblast, and 2,527 unresolved_stromal. This is a lineage freeze only; global V8 remains NOT_FROZEN and V7 remains default.

## NEXT

Run lightweight Endothelial and B/Plasma validation with the same fixed biological template.

External evidence: ${PROJECT_ROOT}/Project_v3/v8_upgrade/biological_freeze_sprint/fibroblast_20260913/FROZEN_V8_FIBROBLAST_TAXONOMY.
