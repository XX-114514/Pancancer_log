# Codex Task: integrate-v8-development-release

- Date: 2026-09-12
- Status: completed; local commit pending at log creation and updated after validation
- Related project stage: P4/P7 release governance; V8 development candidate only
- Related run ID: `20260912_v8_full_metadata_candidate_v7_rerun2`
- Requested by: user

## User request

Create a lightweight, auditable project-records release/run index for the 2026-09-12
V8 development software candidate. Retain V7 as the frozen default, do not modify
Project_v3, do not copy the 251,914,833-byte pickle, and record only external
logical paths, SHA-256 values and status. Include the full candidate rerun2,
statistics/report, scATOMIC salvage report/JSON, run policy v4 and pilot replay
attestation; update status/planning/inventories and make a local commit without
pushing.

## Initial state

- `project-records` was at commit `0e6f3c9` with two pre-existing, untracked V8
  task logs dated 2026-09-11. They were preserved for review rather than replaced.
- V7 was the recorded `FROZEN_PROJECT_REFERENCE` and the default identity/CNV
  reference. The communication release remained V5-derived.
- The task scope was this independent records repository only; Project_v3 evidence
  was inspected read-only.

## Plan

1. Read repository rules, current authority documents, templates, relevant V7
   release/run records, prior V8 logs and Git state.
2. Inspect only the named Project_v3 V8 evidence and recompute required hashes.
3. Build an external-only V8 release, run record and structured inventory entries.
4. State the non-promotion and non-biological boundaries in current status and
   planning documents.
5. Validate, inspect the diff, commit reviewed files locally, and do not push.

## Files inspected

- Parent and repository `AGENTS.md`; README, STATUS, TODO, ROADMAP, CHANGELOG;
  data, version, reproducibility, workflow, directory and Git-maintenance policies.
- Release/run/log templates, release and run indexes, artifact inventory, V7 release
  manifest/run record, relevant V7 method, validation scripts and recent V8 logs.
- External V8 rerun2 summary, statistics, audit report, policy v4, scalable-pilot
  replay attestation, scATOMIC retry2 salvage report/JSON, and final reproducible
  scATOMIC supportive-audit v2 JSON.

## Commands executed

```bash
git status --short
git log --oneline -12
find "${PROJECT_ROOT}/Project_v3" -type f -name '*full*candidate*'
sha256sum "${PROJECT_ROOT}/Project_v3/runs/20260912_v8_full_metadata_candidate_v7_rerun2/annotation_v8_full_candidate_v7.pkl.gz" \
  "${PROJECT_ROOT}/Project_v3/runs/20260912_v8_full_metadata_candidate_v7_rerun2/run_summary_v7.json" \
  # plus the named statistics/report/policy/attestation/salvage files
python -c '<read-only JSON status and contract extraction>'
bash scripts/validate_repository.sh
bash scripts/validate_releases.sh
git diff --check
```

## Findings

### Observed facts

- Rerun2 records `PASS_FULL_METADATA_SOFTWARE_CANDIDATE_V7` and `NOT_FROZEN`.
  It records V7 as default, full-candidate-only mode, no automatic promotion, no
  formal downstream mutation, no independent accuracy claim and no new classifier
  execution claim.
- Direct SHA-256 recomputation matched the candidate pickle
  `23fbf7e21768687bd4df89558462db04be710afd2b9a5a7c37db70e6ed69567c` and
  summary `4dfa482c82adcc28356d0aecbc1f1782be395be3f812b1717e89798d678d998f`.
  It also produced the named small-evidence checksums in the external manifest.
- The statistics audit records 4,676,787 cells, 1,322 samples, 54 GSE, 43 cancers,
  matching frozen dimensions and index hash, and 128/128 V7 field contracts passed.
- The scATOMIC supporting manifest records a `PASS` read-only salvage with no
  algorithm rerun, network, install, overwrite or patient merge; it is not an
  accuracy claim or a scale-up result. The final reproducible supportive audit v2
  is `PASS_SUPPORTIVE_ONLY`, scope-unverified and has no candidate deserialization.
- During indexing, the external full-candidate audit report changed from 13,147
  bytes / `cca46a…253b2` to 19,021 bytes / `2a5f86…108cc`. Review of its current
  content confirmed an append-only v2 supportive-audit correction, not a changed
  candidate; candidate, summary, statistics, policy, attestation and salvage hashes
  remained stable. The canonical index now pins the current report and v2 JSON.
- The two pre-existing 2026-09-11 logs had no detected credential, absolute internal
  path or SSH-style remote pattern and were compatible with the inspected V8
  summary. They are historical records rather than a complete current task template;
  their original text is retained unchanged. Their recorded 48/48 test claim was
  not independently rerun in this indexing task.

### Interpretation

- The evidence supports a software-contract candidate and its provenance index, not
  a biological V8 release. V7 remains authoritative and downstream history remains
  append-only.

### Assumptions

- None. Fields not available in the inspected run summary (exact command, code
  commit, host and scheduler ID) are explicitly recorded as not recorded.

### Recommendations

- Require a separate versioned biological-validation run and explicit promotion
  decision before any V8 reference or downstream use is claimed.

## Changes made

- Added `releases/annotation_v8_development_20260912/` with README, VERSION,
  data review and external-only artifact manifest.
- Updated the report entry after its observed append-only v2 correction and added
  the final reproducible scATOMIC supportive-audit v2 JSON to the canonical index.
- Added the rerun2 lightweight run record; updated README, STATUS, TODO, ROADMAP,
  CHANGELOG, release/run/artifact inventories and release index.
- Included the reviewed pre-existing 2026-09-11 V8 logs without rewriting them.
- Did not create, modify, copy or delete any Project_v3 artifact.

## Validation

- `bash scripts/validate_repository.sh`: PASS, 0 errors and 0 warnings; it also
  passed placeholder, privacy/path, file-size, TSV and local-link checks.
- `bash scripts/validate_releases.sh`: PASS, 0 release-validation errors.
- `VERSION.json` parsed successfully; `git diff --check` passed.
- A read-only manifest validator rehashed all 9 indexed external files, including
  the 251,914,833-byte candidate and final scATOMIC supportive-audit v2 JSON; every
  observed byte size and SHA-256 matched the current manifest.
- The initial report-hash mismatch correctly detected its concurrent append-only v2
  correction; the manifest was updated and the final all-nine recheck passed.
- The staged diff was reviewed and `git diff --cached --check` passed before the
  local commit; no push is planned.

## Failed attempts

- The standard sandbox shell could not create a user namespace on this host; scoped
  read-only commands were therefore run with approved host execution.
- `jq` was unavailable. A read-only Python JSON parser was used instead; no external
  artifact was written or changed.

## Unresolved issues

- Independent identity/malignancy accuracy, new classifier execution, qualified
  orthogonal evidence and biological V8 promotion are not established by this task.
- The historical 2026-09-11 V8 logs are not template-complete; this is recorded as a
  provenance-format limitation, not silently rewritten.

## Recommended next action

- Keep V7 frozen/default and undertake only a new, append-only, independently
  validated V8 biological decision workflow when its prerequisites are available.

## Proposed commit message

```text
release: index V8 development software candidate
```
