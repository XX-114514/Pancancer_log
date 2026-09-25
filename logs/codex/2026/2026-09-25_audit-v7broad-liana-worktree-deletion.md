# Codex Task: audit missing V7 broad CM/LIANA files and commit project records

- Date: 2026-09-25
- Status: completed_and_verified
- Related project stage: P5/P7 downstream governance
- Related run ID: `20260923_v7broad_covarnet_liana`
- Requested by: project owner

## User request

Check whether the prior CM/LIANA result files were accidentally deleted and commit
the corresponding update in the project-records repository.

## Initial state

- The project-records repository was clean at `fabd2d9`.
- The external analysis repository was dirty and contained many unrelated untracked
  production outputs.
- The target `20260923_v7broad_cm_liana_v1` directory was absent from the analysis
  worktree.

## Plan

1. Verify physical files, Git status, completion commit and current Git objects.
2. Separate tracked aggregate evidence, ignored native intermediates and external
   CoVarNet outputs.
3. Inspect the corrected annotation master and cleanup policy before attributing cause.
4. Record the run, incident, current status, inventories and follow-up without
   restoring or committing the external deletion.
5. Validate and commit only the project-records changes.

## Files inspected

- Project-records `AGENTS.md`, `README.md`, `STATUS.md`, `ROADMAP.md`, `TODO.md`,
  `CHANGELOG.md`, relevant methods/docs, recent release/run/Codex records and templates.
- `${PROJECT_ROOT}/Project_v3/annotation_rebuild/.git` history/status for commit
  `3c2234e` and the current HEAD.
- Corrected master `README_CN.md` and `CLEANUP_CANDIDATES_CN.md`.
- External CoVarNet output contract and selected hashes.

## Commands executed

```bash
git status --short
git show --stat --oneline --summary 3c2234e
git ls-tree -rl HEAD -- 20260923_v7broad_cm_liana_v1
git ls-files -d
git cat-file -e HEAD:<tracked-result-path>
git merge-base --is-ancestor 3c2234e HEAD
find <logical-result-roots> -type f
sha256sum <selected-small-artifacts>
bash scripts/validate_repository.sh
```

## Findings

### Observed facts

- All 81 tracked files in the target directory are absent from the worktree and are
  unstaged deletions; 120 tracked files are absent repository-wide.
- No target deletion is staged. The current Git tree still contains all 81 target
  files (19,609,822 bytes); the aggregate subtree contains 28 files (12,477,538 bytes).
- Completion commit `3c2234e` exists and is an ancestor of the current annotation HEAD.
- The historical `.gitignore` excluded native LIANA results and metadata caches;
  those directories are absent and have no Git objects.
- The external CoVarNet directory remains present with 35 files / 315,143,761 bytes
  and all checked core contract files nonempty.
- The corrected master explicitly superseded the old operational gate and listed the
  target directory among manual cleanup candidates. It also states agents did not
  automatically delete the listed directories.

### Interpretation

- The tracked reports, figures, scripts and aggregate tables were deleted from the
  filesystem but not lost from Git.
- The native per-sample LIANA intermediates are genuinely unavailable unless another
  copy exists outside the inspected scope or the run is repeated.
- The matching cleanup list makes manual cleanup plausible, but the executor and
  intent cannot be established from inspected evidence.
- This historical run cannot serve as the current V7-derived release because its
  annotation input was later corrected.

### Assumptions

- No assumption was made about who executed the deletion.
- No uninspected backup is assumed to exist.

### Recommendations

- Keep `REL-003` open and use the corrected master for any current downstream rerun.
- Decide explicitly whether to restore tracked historical evidence or commit a formal
  retirement. Preserve current bulk deletions unstaged until that decision.

## Changes made

- Added a formal run record and deletion incident.
- Added the 2026-09-25 daily log.
- Updated STATUS, TODO, CHANGELOG and run/artifact inventories.
- Did not modify or restore the external annotation repository.

## Validation

- `git diff --check` passed.
- `bash scripts/validate_repository.sh` completed with 0 errors and one pre-existing
  class of path-scan warning on logical-path inventory files; new records contain no
  real username, `/data4` path, credential or patient-level payload.
- TSV column validation passed for both modified inventories; no candidate file
  exceeded 5 MiB and no broken local Markdown link was detected.
- External evidence checks are detailed in the run and incident records.

## Failed attempts

- The initial relative-path check used the Project_v3 root incorrectly and returned
  no target; it was corrected by entering the nested annotation repository.
- The installed legacy Git does not support `git -C` or `git branch --show-current`;
  subsequent checks used an explicit working directory and `git rev-parse`.

## Unresolved issues

- Deletion executor and exact deletion time remain unverified.
- Native per-sample LIANA outputs are not Git-recoverable.
- Restore-versus-retire decision remains open.

## Recommended next action

Resolve `REL-003A`; if current downstream is needed, start a new corrected-master run.

## Proposed commit message

```text
records: audit superseded V7 broad LIANA deletion
```
