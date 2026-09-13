# Codex Task: v8-fibroblast-biological-refinement

- Date: 2026-09-13
- Status: completed
- Related project stage: Biological Freeze Sprint P0-5
- Related run ID: 20260913_fibroblast_3gse_v1
- Requested by: Project user

## TASK

Execute and record the Fibroblast biological gate; freeze only identities/states supported by marker, recurrence, validation and source-label evidence.

## INPUT

The central records repository and external Project_v3 Fibroblast preregistration, fixed V7 candidate universe, 359-gene full-library expression object, six selected original samples, completed cluster evidence and source-label artifacts.

## EXECUTION

Read the records repository instructions and current authoritative files. Bound the runner to pre-unblinding rule hashes, compile-checked it, and executed the fixed six-sample run in the scanpy environment. Recomputed recurrence only among identity-positive clusters, exact-mapped all 10,147 cells to source labels, and generated cluster/cell freeze decisions plus a hash-bearing sentinel. No threshold was changed after validation unblinding and no V7 field was written.

## RESULT

All 6 samples across three GSE succeeded: 10,147 cells and 30 clusters. Fibroblast identity passed in 7,620 cells/16 clusters/3 GSE/6 samples. Exact source mapping was 10,147/10,147; source major labels supported Fibroblast/stromal in 10,060 cells and explicit competing lineage in 86, with one mixed. myCAF alone was accepted for downstream use (2,748 cells/10 clusters/3 GSE/4 samples; 1,055 validation cells). Final assignments are 2,748 Fibroblast_myCAF, 4,872 parent Fibroblast, and 2,527 unresolved_stromal.

## DECISION

Record FROZEN_V8_FIBROBLAST_TAXONOMY. iCAF, ECM-remodeling, antigen-presentation-like, IFN and cycling remain parent/overlay evidence and are not CoVarNet nodes. This does not freeze global V8; V7 remains default.

## NEXT

Run lightweight Endothelial and B/Plasma validation as P1-1.

## Files inspected

- AGENTS.md, README.md, STATUS.md, ROADMAP.md, TODO.md
- Latest T/NK run/log and inventories/runs.tsv
- External preregistration, runner, run summary, cluster/state tables, source-label mapping and Fibroblast freeze sentinel

## Commands executed

- Compiled the Fibroblast runner in scanpy.
- Ran the fixed six-sample biological validation.
- Ran and repaired the source-label evidence mapper.
- Generated and asserted the frozen taxonomy, recurrence, cluster and cell tables.
- Ran the central repository validator.

## Changes made

- Added the Fibroblast run record and this Codex task log.
- Updated authoritative status and the nine-line sprint board.
- Appended the external run to inventories/runs.tsv.
- External Project_v3 artifacts contain the taxonomy, recurrence, cluster and cell decisions, execution record and freeze sentinel.

## Validation

External assertions passed for 6/6 samples, 10,147 cells, 30 clusters, 100% exact source mapping, and downstream-state counts 2,748/4,872/2,527. Repository validation and final diff review were run after this edit.

## Failed attempts

The workspace patch helper failed because the host kernel could not create its namespace; narrow fallback edits were used. The first source-label join used untransformed global IDs and mapped 0 cells; it was corrected to the literal source-cell suffix and reached 10,147/10,147. Two evidence-collapse bugs were corrected before decision: competing fallback needed precedence over V7 Fibroblast, and source literal CAFs required plural-aware matching. A delegated freeze-artifact builder stalled without writing files and was interrupted; the same fixed decision was generated directly and validated.

## Unresolved issues

Global V8 is not frozen. Endothelial/B-Plasma validation, held-out benchmark, frozen-identity malignancy rebuild, and the single final robustness audit remain outstanding.

## Recommended next action

Execute P1-1 Endothelial and B/Plasma validation; do not resume software-governance work.

## Proposed commit message

analysis: freeze V8 Fibroblast taxonomy after biological validation
