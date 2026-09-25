# Incident: V7 broad CM/LIANA worktree directory absent after cleanup proposal

## Summary

On 2026-09-25, the tracked directory
`${PROJECT_ROOT}/Project_v3/annotation_rebuild/20260923_v7broad_cm_liana_v1/`
was found absent from the filesystem. Git reports all 81 tracked files in that
directory as unstaged deletions. The files remain in the current Git object tree.

## Impact

- Historical scripts, aggregate tables, reports and figures are not directly
  accessible from the worktree but are exactly recoverable from Git.
- The ignored native per-sample LIANA outputs and metadata cache are absent and were
  never committed, so they are not recoverable from Git.
- The external CoVarNet result directory remains physically present and passed its
  core output-contract existence check.
- This run was already superseded scientifically by the corrected annotation master;
  no current final downstream release is lost.

## Detection

The user asked whether the prior result files had been deleted. `find` returned no
target directory and `git status --short` showed `D` for every tracked target file.

## Timeline

| Time | Event | Evidence |
| --- | --- | --- |
| 2026-09-23 03:03 +08:00 | LIANA aggregate completed | completion commit `3c2234e` and Git-stored manifest |
| 2026-09-23 11:06 +08:00 | corrected master documented the old branch as superseded and listed a manual cleanup command | `20260923_annotation_master_v1/README_CN.md`; `CLEANUP_CANDIDATES_CN.md` |
| 2026-09-25 | missing worktree directory detected | 81 target tracked deletions; 120 tracked deletions repository-wide |

## Evidence

- `3c2234e` is an ancestor of the current annotation repository HEAD.
- Current Git tree contains 81 target files totaling 19,609,822 bytes.
- Aggregate subtree contains 28 files totaling 12,477,538 bytes.
- Git-object SHA-256 values for the aggregate manifest and interpretation report are
  recorded in `runs/20260923_v7broad_covarnet_liana.md`.
- The cleanup document states that agents did not delete directories and that the
  listed `rm -rf` operations required manual execution after acceptance.

## Root cause

**Unconfirmed.** The cleanup document contains the exact command that would remove
the target and several other currently absent legacy directories, so execution of
that cleanup list is a plausible explanation. There is no inspected evidence that
identifies who executed it or proves whether the deletion was intentional.

## Contributing factors

- The directory was marked superseded after a strict naming-gate correction.
- Large native LIANA results were intentionally ignored by Git.
- Cleanup was not followed by an archive manifest or a committed retirement change.

## Resolution

No restoration or deletion commit was performed during this audit. Git objects and
the external CoVarNet outputs were left unchanged.

## Validation

- `git cat-file -e HEAD:<path>` passed for the aggregate manifest, interpretation
  report and core90 consensus table.
- All required CoVarNet tables checked in the external output directory were nonempty.
- Project-records validation completed with 0 errors and one logical-path inventory
  warning.

## Prevention

- Before manual cleanup, record tracked/untracked scope, artifact hashes and the
  decision to retain or retire the run.
- Never treat ignored per-sample results as recoverable merely because aggregate
  tables are committed.
- Commit approved retirements separately; do not leave broad unstaged deletions in
  a long-lived analysis worktree.

## Follow-up tasks

- `REL-003A`: decide tracked-evidence restoration versus explicit retirement, and
  rerun from the corrected master if a current downstream release is required.
