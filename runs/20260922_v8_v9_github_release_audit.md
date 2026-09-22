# Run: 20260922_v8_v9_github_release_audit

## Metadata

- Run ID: `20260922_v8_v9_github_release_audit`
- Run type: lightweight release audit and GitHub synchronization preparation
- Date and timezone: 2026-09-22, Asia/Shanghai
- Status: `COMPLETED_AND_PUSHED`
- Content commit: `834af7afc577ef2e197cfb6d2c8f4936ebd6f43c`
- Analysis computation: none; external artifacts inspected and copied read-only
- Global identity reference: V7

## Objective

Publish the recent V8 lineage audit qualifications and the latest safe V9
candidate annotation records without uploading cell-level data, internal absolute
paths or incomplete results as frozen truth.

## Observed inputs

- Five V8 lineage manifests and the generated 2026-09-16 V8 status.
- V8 held-out-exposure, contamination, clustering-sensitivity, marker-panel and
  DecontX audit summaries.
- V9 final report dated 2026-09-19, freeze manifest v2 and open L3 review template.

## Release decision

- Copy the six compact V8 status/manifest files and index related audits externally.
- Copy the V9 final report and open review template.
- Do not copy the V9 freeze manifest because 78 keys contain internal absolute paths.
- Do not copy matrices, per-cell outputs, detailed runtime logs or path-bearing reports.
- Preserve V7 as the only global identity/CNV reference.

## Outputs

- `releases/annotation_v8_lineage_audit_20260917/`
- `releases/annotation_v9_candidate_20260919/`

## Validation

- All eight copied source files were byte-identical to their external sources.
- Repository validation completed with 0 errors and 0 warnings; release validation
  completed with 0 errors.
- JSON syntax, TSV shape, relative links, path/credential scans, the 5-MiB gate and
  staged/full-range diff checks passed.
- A normal push updated private `origin/main` from `ad07778` to `834af7a`; a fresh
  authenticated `ls-remote` returned the same full SHA as local `HEAD`.

## Promotion boundary

V8 is still not globally frozen. V9 remains a candidate with 15 open L3 reviews,
632 suggested-but-unexecuted CNV reruns and unresolved validation/stability gates.
