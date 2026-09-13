# Run: 20260913_v8_tnk_biological_refinement

## TASK

Does a rule-frozen, three-GSE T/NK review support a conservative V8 lineage taxonomy and stable downstream state?

## INPUT

Frozen V7 candidate universe; original/raw expression by GSE/sample; source T/NK/other labels; fixed lineage-specific marker workflow; genuine scATOMIC 2.0.3 predictions from the prior real run.

## EXECUTION

T/NK acceptance, parent-fallback and unresolved rules were frozen before validation review (2026-09-13T16:03:57+0800; rules TSV SHA-256 `3396c64cab20df27503c84577f23075974b171b4a777a5e861764ba1917baa08`). The primary `20260913_tnk_3gse_v3` run completed 6/6 samples across three GSE. Source-to-review mapping was exact. Genuine scATOMIC evidence was mapped separately from the prior real run; no adapter/interface/schema result was promoted as execution.

## RESULT

The primary run covered 14,185 cells and 35 clusters. Exact source mapping was 14,185/14,185; source labels were T/NK/other = 12,534/1,258/393. Review outcomes were parent 10,664, boundary 2,845, resolved-unresolved 180 and conflict 103. No fine identity achieved cross-GSE acceptance; CD4 naive-memory reached only 1 GSE/1 sample/1,271 cells. The accepted conservative identity space is Tcell, NK, T_or_NK and unresolved. Tissue-resident-like is accepted only as a state (3 GSE/5 samples/6,574 cells); other fine identities/states fall back or remain unresolved.

The genuinely executed scATOMIC subset mapped 4,719/4,719 cells across two samples: 4,413 confident, 3,705 supporting, 203 contradictory, 17 ambiguous, 17 auxiliary and 471 not evaluable.

## DECISION

`FROZEN_V8_TNK_TAXONOMY`: freeze this lineage taxonomy and conservative state boundary only. V1/V2 are excluded from the accepted evidence because they contained two contamination-overlap bugs; V2→V3 decision columns were identical after correction. This is a lineage taxonomy freeze, not a global V8 cell-level freeze; V7 remains the default and global V8 remains `NOT_FROZEN`.

## NEXT

Run Fibroblast refinement using the same frozen biological template.

External evidence: `${PROJECT_ROOT}/Project_v3/v8_upgrade/biological_freeze_sprint/tnk_20260913/FROZEN_V8_TNK_TAXONOMY`.
