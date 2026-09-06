# Integrate V7 and current progress for Git

## User request

Maintain a Git repository under the PanCancer project, integrate the prior
Pancancer_log history and current progress, and make it suitable for version
control, Prism, manuscript writing and supplementary submission.

## Initial observed state

- `project-records` was already a Git repository with the
  `XX-114514/Pancancer_log.git` remote.
- Local branch `main` was one commit ahead of `origin/main`.
- `STATUS.md` had one pre-existing uncommitted blank-line insertion.
- Adjacent Project_v1-v5 trees total more than one terabyte and contain nested Git
  repositories, environments and multi-gigabyte scientific objects.
- Git LFS was not installed.

## Decision and actions

- Preserved the existing Pancancer_log history instead of creating a competing
  parent-level repository.
- Created append-only lightweight releases for V7 and the 2026-09-06 communication
  evidence chain.
- Copied only reviewed reports, aggregate tables, figures, config and path-safe
  source files.
- Ignored local script copies containing workstation-specific paths; recorded
  their logical provenance instead.
- Added Git maintenance, versioning and Prism/supplementary guidance.
- Updated run/release/artifact indexes and current status without deleting
  historical records.

## Validation plan

- Check all candidate files for secrets, absolute/internal paths and size.
- Check TSV column consistency and Markdown links.
- Compile copied Python code where applicable.
- Compare copied release files against source hashes.
- Review staged diff before commit.

## Failed attempts

The preferred patch helper intermittently failed because the host kernel does not
support its sandbox namespace. A system patch was used only for additive changes.
An attempted path-refactor patch was rejected by the workspace no-deletion rule;
no source analysis script was altered.

## Suggested commit message

`release: add V7 reference and manuscript evidence snapshots`
## Validation result

- Repository validation: 0 errors, 0 warnings.
- Release hash validation: 0 errors; all copied V7 and evidence-chain files matched
  their source copies.
- JSON syntax: passed.
- TSV column consistency: passed.
- V7 Python publication copies: byte-compiled successfully.
- Releasable-file internal path scan: passed.
- Largest tracked candidate: 1,180,313 bytes, below the 5 MiB policy limit.
- Sample-level V7 tables were reviewed as public accession-derived method summaries
  without clinical or direct identity fields.
- Remote push was intentionally not performed; the user requested future upload.
