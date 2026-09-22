# Codex Task: push V8/V9 annotation and audit records

- Date: 2026-09-22
- Status: ready_for_commit_and_push
- Related run: `20260922_v8_v9_github_release_audit`

## User request

Confirm and push the recent V8/V9 annotation and latest audit files to GitHub.

## Initial state

- Local `main` was at `52dcf6f`, six commits ahead of `origin/main` at `ad07778`.
- Those six commits contain the five V8 lineage-taxonomy freezes plus the prior
  synchronization audit.
- The prior Nature upstream-freeze package was present but uncommitted.
- No files were staged.

## Plan

1. Verify the configured private GitHub destination and authenticated remote head.
2. Audit recent external V8/V9 files for status, paths, size and disclosure risk.
3. Create append-only lightweight releases for the safe annotation/audit evidence.
4. Validate and review the entire staged range before a normal push.

## Key findings

- Remote `main` was verified at `ad0777835c0ddee516643c02574c18b8c843b2c8`.
- V8 remains globally `NOT_FROZEN`; T/NK held-out evidence is compromised and
  Myeloid/DC held-out outputs existed before the split was frozen.
- V9 remains `FROZEN_CANDIDATE`, with 15 open L3 reviews and 632 suggested CNV
  reruns that were not executed.
- The V9 final report and review template are safe lightweight copies. The original
  freeze manifest is not copied because 78 keys contain internal absolute paths.

## Changes prepared

- Added V8 lineage-audit and V9 candidate append-only releases.
- Added run/log records and release/artifact indexes.
- Preserved the earlier Nature manuscript upstream-freeze package in the same
  reviewed publication-record commit.

## Validation

- All eight copied V8/V9 files are byte-identical to their reviewed sources.
- JSON parsing, targeted credential/path scans and the 5-MiB file gate passed.
- `scripts/validate_repository.sh`: 0 errors and 0 warnings.
- `scripts/validate_releases.sh`: 0 release-validation errors.
- Staged and full remote-to-worktree `git diff --check` passed. The exact V9 TSV
  source retains four structurally meaningful trailing empty fields under a
  file-specific Git whitespace attribute.
- Post-push remote verification remains pending.

## Failed attempts

- `gh` is not installed.
- Default Git SSH and legacy Git's unsupported `GIT_SSH_COMMAND` path selected the
  wrong identity and failed. A temporary `ssh-agent` with the previously verified
  Ed25519 identity authenticated successfully; no remote setting was changed.
- Adding the file-specific whitespace attribute initially displaced existing
  `.gitattributes` rules. Staged review caught it; all original rules were restored
  before commit and only the narrow V9 TSV exception remains appended.

## Unresolved scientific issues

- V8/V9 promotion gates remain open; this GitHub publication does not close them.
- V9 manual decisions and CNV reruns are not performed by this task.

## Proposed commit message

```text
release: publish V8 and V9 annotation audit snapshots
```
