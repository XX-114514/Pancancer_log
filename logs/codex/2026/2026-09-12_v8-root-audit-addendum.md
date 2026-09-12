# Codex Task: v8-root-audit-addendum

- Date: 2026-09-12
- Status: validated; local commit pending
- Related project stage: P4/P7 release governance; V8 development candidate only
- Related run ID: `20260912_v8_full_metadata_candidate_v7_rerun2`
- Requested by: user

## User request

Create a small, Project_v3-external V8 root-audit addendum in `project-records`.
Independently recompute paths, sizes, and SHA-256 values for five named root-audit
and independent-reference artifacts. Preserve the existing V8 release manifest;
state the development-only root gate, V7 default, 167 recorded tests, blocked
independent reference, and absence of a biological freeze. Validate and commit
locally without pushing.

## Initial state

- `project-records` HEAD was `137249d` (`release: index V8 development software
  candidate`) and `git status --short` was clean.
- The existing V8 release is `NOT_FROZEN`; V7 is the frozen/default identity and
  CNV reference.
- The existing V8 `external_artifacts.tsv` is a published candidate index. This
  task therefore adds a distinct, dated addendum rather than modifying that file.

## Plan

1. Read repository governance, the V8 release/run/log context, and the referenced
   historical commit.
2. Inspect the five named Project_v3 files read-only and directly recompute their
   SHA-256 values and byte sizes.
3. Add a dated external-evidence addendum and global artifact-inventory rows while
   leaving the original V8 release manifest unchanged.
4. Run repository, release, external-manifest, and Git diff validation; review and
   commit only the scoped records files without pushing.

## Files inspected

- `AGENTS.md`, `README.md`, `STATUS.md`, `ROADMAP.md`, `TODO.md`, relevant
  policies, validators, release files, run record, inventories, and V8 Codex logs.
- Historical release commit `137249d`.
- External Project_v3 root audit, root gate, independent-reference report, V2
  Census metadata contract, and V2 candidate inventory named in the request.

## Commands executed

```bash
git status --short
git show --stat 137249d
sha256sum <five named Project_v3 root-audit/reference artifacts>
stat --format='%s %n' <the same five external artifacts>
bash scripts/validate_repository.sh
bash scripts/validate_releases.sh
<read-only five-row external size/SHA-256 manifest validator>
git diff --check
```

## Findings

### Observed facts

- The external root-release gate reports
  `PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`, V7 as default, automatic
  promotion and formal downstream mutation both false, and `pytest_passed: 167`.
- The independent-reference contract and inventory report
  `BLOCKED_NO_LOCAL_CELLXGENE_CENSUS_CLIENT`; they do not report an executed
  reference mapping, model training, prediction, or truth set.
- The five directly recomputed byte sizes and SHA-256 values are recorded in
  `ROOT_AUDIT_EXTERNAL_ARTIFACTS_20260912.tsv` and match the observed root-audit
  references supplied for the first two files.

### Interpretation

- The new evidence strengthens provenance for a development software candidate,
  not a biological V8 release. V7 remains the project default.

### Assumptions

- None. Project_v3 test execution is reported from the external root gate rather
  than inferred or rerun by this records-indexing task.

### Recommendations

- Keep V7 frozen/default. Consider a V8 biological promotion only in a new,
  versioned run after independent reference/truth and downstream-impact gates are
  completed and explicitly reviewed.

## Changes made

- Added the dated V8 root-audit addendum and its five-artifact external index.
- Added five corresponding provenance rows to `inventories/artifacts.tsv`.
- Kept the existing V8 release `README.md`, `VERSION.json`, and
  `external_artifacts.tsv` unchanged.

## Validation

- `bash scripts/validate_repository.sh`: PASS, 0 errors and 0 warnings.
- `bash scripts/validate_releases.sh`: PASS, 0 release-validation errors.
- Read-only external-manifest validator: PASS; all five indexed paths, byte sizes,
  and SHA-256 values matched the Project_v3 files.
- `git diff --check`: PASS. The original V8 release manifest was not modified.

## Failed attempts

- The default sandbox shell could not create a user namespace on this host.
  Scoped host execution was used for read-only evidence inspection after approval.

## Unresolved issues

- Independent scANVI reference mapping, held-out accuracy evidence, qualified
  orthogonal malignancy evidence, and biological V8 promotion remain unavailable.

## Recommended next action

- Create the scoped local commit; do not push.

## Proposed commit message

```text
docs: add V8 root-audit provenance addendum
```
