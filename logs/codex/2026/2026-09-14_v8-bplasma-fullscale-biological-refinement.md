# Codex Task: v8-bplasma-fullscale-biological-refinement

- Date: 2026-09-14
- Status: completed
- Related stage: Biological Freeze Sprint P1-1 B/Plasma
- Related run: `20260914_bplasma_full_v1`

## TASK

Replace the proposed lightweight B/Plasma gate with a feasible full-cohort run, produce actual biological evidence and stage the 0914 progress in the central records repository.

## INPUT

Central project records; frozen V7 B/Plasma universe; frozen 359-gene expression object; preregistered rules; prospectively sealed held-out GSE list; existing genuine scATOMIC predictions.

## EXECUTION

Counted the exact candidate universe and available marker genes, narrowed panels to structurally available canonical genes, froze rules and the held-out split before predictions, ran a two-sample smoke, corrected one genuine parent-state gating bug, reran the smoke, then processed all 1,312 original samples. Reviewed only development outputs, mapped scATOMIC and source IDs exactly, generated a deterministic taxonomy/cell/cluster finalizer and recorded the lineage freeze.

## RESULT

Full run: 981 success, 331 explicit skips, 0 failed, 685,843 cells and 4,378 clusters. Development review: 783 samples, 39 GSE, 557,601 cells and 3,501 clusters; 11/11 contracts PASS. Bcell and Plasma identities plus B-naive state passed; Plasma-IgA fell back under the frozen universal-state rule. The finalizer independently verified 557,601 unique cell IDs and exact downstream-state count closure. Genuine scATOMIC mapping was 1,015/1,015; source mapping was 2,338/2,338, with source biological labels unavailable.

## DECISION

Create `FROZEN_V8_B_PLASMA_TAXONOMY`. Mark P1-1 DONE/FROZEN and activate P1-2 held-out benchmark. Do not unseal benchmark labels during taxonomy tuning, write back V7, or resume schema/audit framework work.

## Files inspected

- Repository instructions, current status, roadmap, task board, V7 method and data/reproducibility policies.
- External preregistration, runner, manifests, development-only evidence, source/classifier mappings and freeze artifacts.

## Commands and validation

- Compiled and ran the fixed runner in the `scanpy` environment.
- Executed exact source/scATOMIC mapping and development-only marker/recurrence review.
- Compiled and ran the deterministic finalizer; independently counted 557,601 unique cell IDs and five final states.
- Ran the central repository validator and reviewed the staged diff before commit.

## Failed attempts

The first smoke exposed Plasma-IgA state acceptance within Bcell identities; the parent-specific gate was fixed before the full run and smoke v2 had zero cross-identity state calls. The finalizer's first invocation used one incorrect relative path; it produced no accepted freeze outputs, the path was corrected, and the complete rerun passed. A delegated finalizer did not start an execution process and was interrupted; it produced no files. The preferred patch helper could not create its host bwrap namespace, so narrow fallback writes were used. Two ordinary GitHub push attempts—default SSH and the existing explicit ed25519 identity—both failed with `Permission denied (publickey)`; no third environment repair was attempted.

## Unverified / unresolved

Held-out predictions remain sealed. GSE188711 contains no source biological cell label, and scATOMIC has no Plasma class; these limitations are explicit. Global V8, V8 malignancy and final robustness remain incomplete.

## NEXT

Run the frozen 12-GSE held-out reference benchmark.

## Proposed commit

`analysis: freeze V8 B/Plasma taxonomy after full-cohort validation`
