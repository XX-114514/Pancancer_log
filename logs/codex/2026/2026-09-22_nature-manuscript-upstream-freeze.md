# Codex Task: nature-manuscript-upstream-freeze

- Date: 2026-09-22
- Status: completed_with_open_scientific_gates
- Related project stage: P7 report and publication
- Related run ID: `20260922_nature_manuscript_upstream_freeze`
- Requested by: user

## User request

在 `project-records` 中确定并冻结当前注释上游，按 Nature 投稿标准准备发表文件。

## Initial state

- Records repository was clean, on `main` at `52dcf6f`, six commits ahead of origin.
- V7 was recorded as the frozen global identity/CNV reference; V8 was explicitly
  not globally frozen; communication figures were V5-derived.
- External Project_v3 contained a newer 2026-09-19 V9 candidate that had not been
  indexed in this repository.

## Plan

1. Audit version authority and external artifacts without rerunning analysis.
2. Lock the defensible upstream with path-safe logical IDs, byte sizes and hashes.
3. Quarantine incomplete candidates and preserve downstream version boundaries.
4. Create a Nature Article submission-preparation package and validate the repository.

## Files inspected

- Repository instructions, status, roadmap, task board, data/version/reproducibility policies.
- V7 release, methods, configuration, run records, external objects and final audits.
- V8 lineage run records and external freeze sentinels.
- V9 final report, candidate manifest and manual-review template.
- Communication evidence release and current figure/table assets.
- Nature official formatting, initial-submission, reporting and declaration pages.

## Commands executed

- Read-only `find`, `rg`, `sed`, `stat`, `du`, `sha256sum`, `git status` and `git log` checks.
- Deterministic raw-Zarr directory-root calculation over 13,696 files.
- Repository and release validators after editing.

## Findings

### Observed facts

- V7 is the only `FROZEN_PROJECT_REFERENCE` global identity/CNV release.
- V8 has five lineage-taxonomy freezes but no global writeback.
- V9 is `FROZEN_CANDIDATE`; 15 L3 clusters need manual review and 632 samples have
  suggested-but-unexecuted CNV reruns.
- V7's prior release did not lock the raw-count Zarr, V2/V5 parents, gene coordinates
  or dynamically loaded CNV engines; this task added their hashes.
- Historical exact environment versions and a Project_v3 code commit remain unavailable.
- Current communication/spatial/clinical figures are V5-derived, not V7/V9-derived.

### Interpretation

Freezing V7 as the manuscript upstream is scientifically defensible; promoting V8
or V9 now would convert incomplete development evidence into an unsupported claim.

### Assumptions

The target venue is the flagship *Nature* Article format, not another Nature
Portfolio journal.

### Recommendations

Close the release README blockers and create a new submission-candidate release.

## Changes made

- Added a manuscript upstream freeze release with V7 dependency locks, V8 hashes,
  V9 quarantine, manuscript/checklist/declaration templates and figure/SI inventories.
- Added a formal run record and updated current-state/release/run/artifact indexes.

## Validation

- `bash scripts/validate_repository.sh`: 0 errors and 0 warnings.
- `bash scripts/validate_releases.sh`: 0 release-validation errors.
- Path-safety, file-size, TSV-shape, relative-link and `git diff --check` scans passed.

## Failed attempts

- Sandboxed reads failed because the host lacks usable user namespaces; approved
  read-only commands were rerun outside that wrapper.
- The preferred patch helper failed to open existing files because its bwrap
  namespace is unavailable; narrow system `patch` edits were reviewed instead.
- The installed Git 1.8.3.1 does not support `git -C`; commands were rerun with an
  explicit working directory.
- The V9 manifest was first addressed at the V9 root; the actual file is under
  `preregistration/`. No file was modified by the failed read.

## Unresolved issues

- Scientific and governance gates listed in the release README.
- Historical environment lock, publication-ready 54-GSE source/licence inventory,
  repository DOI/tag and remote synchronization.

## Recommended next action

Complete V9/V7 version decision and version-consistent downstream reruns before
populating the final summary paragraph and main figures.

## Proposed commit message

```text
release: freeze manuscript upstream and add Nature submission package
```
