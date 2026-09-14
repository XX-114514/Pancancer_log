# Codex Task: v8-endothelial-biological-refinement

- Date: 2026-09-14
- Status: completed
- Related project stage: Biological Freeze Sprint P1-1 Endothelial
- Related run ID: 20260914_endothelial_3gse_v1
- Requested by: Project user

## TASK

Execute and record the lightweight Endothelial biological gate; freeze only identities/states supported by preregistered marker, recurrence, validation and independent-classifier evidence.

## INPUT

Central project records; frozen V7 Endothelial universe; 359-gene full-library expression object; preselected six-sample cohort; existing genuine scATOMIC GSE274229 predictions.

## EXECUTION

Verified the runner uses the V7 Endothelial universe only, removed unavailable EMCN/ENG from the identity panel, preregistered immutable marker/recurrence/fallback rules, compiled the runner and executed all six samples. Reviewed identity-positive recurrence, state overlap, top-marker coherence, source metadata availability and exact scATOMIC overlap. Generated deterministic cell, cluster, recurrence and taxonomy freeze artifacts. No threshold changed after unblinding and no global annotation was written.

## RESULT

All 6/6 samples succeeded: 8,951 cells and 35 clusters. Endothelial identity passed for 7,881 cells/25 clusters across 3 GSE/6 samples. Exact scATOMIC mapping covered 1,547/1,547 S37 cells; 1,538 supported Endothelial identity, 2 were competing and 7 broad. Source objects had no biological cell-label fields, so source labels were recorded as unavailable rather than inferred. Final assignments are 4,583 Endothelial, 1,580 Endothelial_arterial, 1,718 Endothelial_venous and 1,070 unresolved.

## DECISION

Record `FROZEN_V8_ENDOTHELIAL_TAXONOMY`. Capillary, angiogenic-tip and lymphatic are parent/overlay evidence and are not CoVarNet nodes. Global V8 remains NOT_FROZEN; V7 remains default.

## NEXT

Execute B/Plasma lightweight validation, completing P1-1.

## Files inspected

- `AGENTS.md`, `README.md`, `STATUS.md`, `ROADMAP.md`, `TODO.md`
- Prior Fibroblast run/log and run inventory
- External Endothelial preregistration, runner, run summary, cluster/state tables, marker tables and source/scATOMIC mapping

## Commands executed

- Syntax-checked and ran the fixed Endothelial six-sample runner in the scanpy environment.
- Calculated identity-positive recurrence, validation, state overlap and top-marker coherence.
- Exact-mapped original source IDs and existing genuine scATOMIC predictions.
- Generated and asserted frozen taxonomy, recurrence, cluster and cell decision tables.
- Ran the central repository validator and reviewed the final diff.

## Changes made

- Added the Endothelial run record and this Codex task log.
- Updated authoritative status, roadmap, changelog and the nine-line sprint board.
- Appended the run to `inventories/runs.tsv`.
- External Project_v3 artifacts contain preregistration, predictions/evidence, final taxonomy and freeze sentinel.

## Validation

External assertions passed for 6/6 samples, 8,951 unique cells, 35 clusters, 1,547/1,547 exact scATOMIC overlap and downstream counts 4,583/1,580/1,718/1,070. Repository validation and final diff review were run after editing.

## Failed attempts

The workspace patch helper failed twice because the host kernel could not create its bwrap namespace; narrow fallback edits were used. The first finalizer assertion used an incorrect parent arithmetic expectation (2,455 instead of 4,583); no biological output was accepted from that failed invocation, the invariant was corrected, and the deterministic finalizer passed. One delegated summary text misstated 9,951 cells; the per-sample tables and run summary consistently show 8,951, which is the recorded denominator.

## Unresolved issues

Global V8 is not frozen. B/Plasma validation, held-out benchmark, frozen-identity malignancy rebuild and the final robustness pass remain outstanding.

## Recommended next action

Execute B/Plasma validation; do not resume software-governance work.

## Proposed commit message

analysis: freeze V8 Endothelial taxonomy after biological validation
