# Codex Task: v8-tnk-biological-refinement

- Date: 2026-09-13
- Status: completed
- Related project stage: Biological Freeze Sprint P0-3
- Related run ID: `20260913_tnk_3gse_v3`
- Requested by: Project user

## TASK

Record the completed T/NK biological refinement and freeze the lineage taxonomy only.

## INPUT

The central `project-records` repository, its current V8/Myeloid records, and the
external Project_v3 T/NK execution artifacts. V7 remains the project default;
the V8 software candidate is not a biological reference.

## EXECUTION

Read `AGENTS.md`, `README.md`, `STATUS.md`, `ROADMAP.md`, `TODO.md`, the related
annotation methods, the latest Myeloid run/log, and the T/NK execution record.
Checked the central worktree before editing. The T/NK rules had already been
frozen before validation review at `2026-09-13T16:03:57+0800` (rules TSV
SHA-256 `3396c64cab20df27503c84577f23075974b171b4a777a5e861764ba1917baa08`).
This records the real `20260913_tnk_3gse_v3` execution; it does not rerun the
analysis or treat adapters, schemas, tests, or historical salvage as execution.

## RESULT

The primary run completed 6/6 samples across three GSE, with 14,185 cells and
35 clusters. Exact source mapping was 14,185/14,185; source T/NK/other counts
were 12,534/1,258/393; parent/boundary/resolved-unresolved/conflict counts were
10,664/2,845/180/103. No fine identity passed cross-GSE recurrence; CD4
naive-memory was limited to 1 GSE/1 sample/1,271 cells. Tissue-resident-like
was supported as a state across 3 GSE/5 samples/6,574 cells.

The independent classifier evidence was a genuine scATOMIC 2.0.3 run mapped
one-to-one for 4,719/4,719 cells in two samples: 4,413 confident, 3,705
supporting, 203 contradictory, 17 ambiguous, 17 auxiliary and 471 not
evaluable. V1/V2 remain traceable but are excluded because of two
contamination-overlap bug fixes; V2→V3 decision columns were identical.

## DECISION

Record `FROZEN_V8_TNK_TAXONOMY`: accept conservative Tcell, NK and T_or_NK
parent identities plus unresolved rejection, and tissue-resident-like on the
state axis only. Other fine identities/states fall back or remain unresolved.
This is not a global V8 freeze, writeback, or accuracy claim.

## NEXT

Run Fibroblast refinement using the same fixed biological template.

## Files inspected

- `AGENTS.md`, `README.md`, `STATUS.md`, `ROADMAP.md`, `TODO.md`
- `methods/project-v3-unified-workflow.md`
- `methods/pancancer-5m-v7-annotation-cnv.md`
- `runs/20260913_v8_myeloid_dc_third_gse_validation.md`
- `logs/codex/2026/2026-09-13_v8-biological-freeze-sprint-myeloid.md`
- External `tnk_20260913/EXECUTION_RECORD.md`, `BUGFIX_LOG.md`, V3 summary,
  and freeze sentinel

## Commands executed

```bash
git status --short
sed -n '1,260p' methods/project-v3-unified-workflow.md
sed -n '1,260p' runs/20260913_v8_myeloid_dc_third_gse_validation.md
sed -n '1,260p' ${PROJECT_ROOT}/Project_v3/v8_upgrade/biological_freeze_sprint/tnk_20260913/EXECUTION_RECORD.md
bash scripts/validate_repository.sh
```

## Changes made

- Added the T/NK run record and this Codex task log.
- Updated the authoritative top biological-status section with the T/NK freeze
  and retained the global `NOT_FROZEN`/V7-default boundary.
- Advanced the nine-line Biological Freeze Sprint board from P0-3 `NEXT` to
  `DONE/FROZEN` and P0-5 to `NEXT`.
- Appended the T/NK run to `inventories/runs.tsv`.

## Validation

The repository validator was run after edits and completed with zero errors and
zero warnings. Final `git status --short` and `git diff --stat` were reviewed;
no commit was created.

## Failed attempts

None for this records-only task.

## Unresolved issues

Global V8 identity is not frozen. Held-out benchmark, Fibroblast and light
Endothelial/B-Plasma validation, V8 malignancy rebuild, and final robustness
checks remain outstanding. V7 remains the default annotation for downstream
analysis.

## Recommended next action

Execute the Fibroblast biological gate; do not add V8 framework or governance
work unless a biological execution bug blocks it.

## Proposed commit message

```text
docs: record V8 T/NK biological taxonomy freeze
```
