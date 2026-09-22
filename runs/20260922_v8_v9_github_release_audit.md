# Run: 20260922_v8_v9_github_release_audit

## Metadata

- Run ID: `20260922_v8_v9_github_release_audit`
- Run type: lightweight release audit and GitHub synchronization preparation
- Date and timezone: 2026-09-22, Asia/Shanghai
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

## Validation plan

- Recompute source/copy hashes.
- Run repository and release validators.
- Review staged diff, file sizes, path safety and GitHub target before commit/push.

## Promotion boundary

V8 is still not globally frozen. V9 remains a candidate with 15 open L3 reviews,
632 suggested-but-unexecuted CNV reruns and unresolved validation/stability gates.
