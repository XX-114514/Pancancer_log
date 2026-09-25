# Run: 20260923_v7broad_covarnet_liana

## Metadata

- Run ID: `20260923_v7broad_covarnet_liana`
- Run type: production CoVarNet + per-sample LIANA + historical comparison
- Start time: 2026-09-22; LIANA attempt 2 started 2026-09-23 01:29:20 +08:00
- End time: 2026-09-23 03:03:05 +08:00
- Timezone: Asia/Shanghai
- Execution status: completed; 1,297/1,297 sample statuses terminal
- Current scientific status: `SUPERSEDED_HISTORICAL_AUDIT_ONLY`
- Current storage status: tracked evidence recoverable from Git; ignored native LIANA outputs absent from worktree
- Code repository: `${PROJECT_ROOT}/Project_v3/annotation_rebuild`
- Completion commit: `3c2234e` (`run: complete V7 CM LIANA comparison`)
- Current containing commit checked on 2026-09-25: `af589d9`; `3c2234e` is an ancestor
- Execution worktree status: not recorded
- Environment: `scanpy` Python 3.10; LIANA 1.7.1; CoVarNet R environment `scrna_r`

## Objective

Use the then-current `V7BROAD_OPERATIONAL_V1` labels to fit a K=9 CoVarNet,
run native-gene per-sample LIANA, compare the new CM structure and communication
edges with the historical V5 branch, and produce reproducible tables and figures.

## Inputs

| Input ID | Logical path | Version/checksum | Role |
| --- | --- | --- | --- |
| V7 broad operational labels | `${PROJECT_ROOT}/Project_v3/annotation_rebuild/20260922_v7_broad_operational_v1/` | Git history before/at `3c2234e` | Cell-state labels used by this run |
| Native full-gene parents | `${PROJECT_ROOT}/Project_v3/annotation_rebuild/20260922_fullgene_v2/` | external run manifest | Per-GSE expression input for LIANA |
| Historical V5 LIANA | `${PROJECT_ROOT}/Project_v3/runs/20260818_pancancer_5M_merge_prep/10_final_v5_malignancy_covarnet_liana/` | historical V5 branch | Matched-CM comparison |

## Input qualification discovered after completion

`20260923_annotation_master_v1/README_CN.md` subsequently documented that the
operational input incorrectly treated `strict_major_eligible_v7` as permission to
name a known V7 broad class. The corrected master separates naming from strict
comparison eligibility and explicitly states that the old operational labels and
their CoVarNet/LIANA downstream are historical audit only. Therefore this run does
not complete the current `REL-003` final-identity downstream requirement.

## Parameters

| Parameter | Value | Source/rationale |
| --- | --- | --- |
| CoVarNet K | 9 | Fixed comparison with prior CM analysis |
| Minimum cells per sample | 100 | CoVarNet run contract |
| LIANA method | consensus/RRA | Workflow manifest |
| `expr_prop` | 0.1 | Workflow manifest |
| minimum cells per group | 20 | Workflow manifest |
| edge rank filter | specificity and magnitude ranks both <=0.05 | Aggregate config |
| primary membership | `core90` | Positive finite W reaching 90% mass, bounded to 2--10 states |
| sensitivity memberships | `positive_top10`, `historic_top10` | Kept separate; no union |

## Execution summary

- CoVarNet completed on 4,055,992 eligible cells, 1,297 samples and 17 states.
- LIANA attempt 1 exited before sample processing because exact DataFrame equality
  was too strict for serialized floating values; no mixed per-sample result was retained.
- LIANA attempt 2 used four workers and reached 1,268 success plus 29 explicit
  `skipped_lt2_groups`, with zero failed or running statuses.
- Aggregation generated separate membership results, CM-abundance associations,
  V5 comparison tables and seven PDF/300-dpi PNG figure pairs.

## Key results at completion

- `core90`: 28,929 qualified consensus edges.
- `positive_top10`: 40,449 edges.
- `historic_top10`: 86,540 edges; historical reproduction only.
- 6,657 sufficiently replicated core90 edge-abundance associations.
- Exact ordered broad sender/receiver plus exact ligand/receptor Jaccard across nine
  matched V7/V5 CM pairs ranged from 0.043 to 0.239: partial shared backbone, not
  nine one-to-one identical modules.
- The completion workflow passed 18 tests; seven gzip result tables and seven
  PDF/PNG figure pairs were validated at completion.

## Current output state audited 2026-09-25

| Artifact | Current state | Validation |
| --- | --- | --- |
| CoVarNet output | physically present at `${PROJECT_ROOT}/Project_v2/CoVarNet/outputs/Covarnet_R_v7broad_marker_v1_20260922/covarnet_analysis_v1/` | 35 files, 315,143,761 bytes; core output contract present |
| CoVarNet W | physically present | SHA-256 `26d0d305d3da2a20b9529eea88015bec06914a6ef6da1d528fcdd435a5338b5a` |
| Entire tracked workflow/result tree | absent from worktree, present in current Git tree | 81 files, 19,609,822 bytes |
| Aggregate LIANA result subtree | absent from worktree, present in current Git tree | 28 files, 12,477,538 bytes |
| Aggregate manifest | present in Git object | SHA-256 `d9b7f51c5a13ff4153863113387ab7c691c61d51a049d707d759c989ef459a1e` |
| Interpretation report | present in Git object | SHA-256 `fcb4aee96b006da63a2b8285ffac61db165a9706457baa05939732745501f4fc` |
| Native per-sample LIANA and metadata caches | absent from worktree and never tracked | historical `.gitignore` excluded `outputs/native_liana_v1/` and `outputs/native_liana_metadata_v1/` |

## Warnings

- This completed run used a broad-label gate later identified as incorrect. It is
  suitable for historical sensitivity/audit, not as the current V7-derived release.
- The tracked lightweight aggregate can be restored exactly from Git. The ignored
  native per-sample results require a new versioned rerun if they are needed.
- Current worktree deletion is uncommitted. Do not commit it as cleanup without an
  explicit archive/retirement decision.

## Interpretation

The execution itself completed and its lightweight aggregate evidence is intact in
Git. Scientific supersession and physical deletion are separate facts: the later
master invalidates this branch as the current downstream input, while Git still
preserves the historical analysis needed for audit and comparison.

## Next action

1. Decide whether to restore only the tracked historical evidence for convenient
   inspection or leave it absent and use `git show`/`git archive` on demand.
2. If current V7-major downstream is required, rerun CoVarNet/LIANA from the
   corrected annotation master under a new run ID; do not revive this run as final.
3. Do not commit the current bulk deletion until the repository owner explicitly
   approves the retirement scope.
