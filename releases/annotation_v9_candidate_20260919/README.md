# PanCancer V9 annotation candidate snapshot

- Release ID: `annotation_v9_candidate_20260919`
- Candidate evidence date: 2026-09-19
- Repository snapshot date: 2026-09-22
- Status: `FROZEN_CANDIDATE_NOT_PROMOTED`
- Current global project reference: V7

This append-only release publishes the latest safe, lightweight V9 annotation
report and open L3 review template. It is a candidate evidence snapshot, not a
frozen V9 biological reference and not a substitute for V7.

## Included files

- [V9 final candidate report](V9_FINAL_REPORT_20260919_CN.md), which supersedes the
  2026-09-18 report.
- [Open L3 review template](manual_review_decisions_L3_pilot_TEMPLATE_v2.tsv), with
  15 unresolved aggregate cluster decisions and blank reviewer/date fields.
- [External artifact index](external_artifacts.tsv).
- [Data review](DATA_REVIEW.md).

The copied report and review template are byte-identical to their Project_v3
sources. The 203,565-byte original freeze manifest remains external because 78
source-hash keys contain internal absolute paths; its size and original SHA-256
are retained in the artifact index.

## Open promotion gates

1. Complete and sign the 15 L3 cluster decisions.
2. Decide and document the 632 suggested CNV reruns; none were executed in V9.
3. Resolve or explicitly accept the preregistered L2/stability failures.
4. Run a formal promotion audit and create a new global release if approved.
5. Rebuild identity-dependent downstream analyses on the promoted version.

Until those gates close, V9 labels must not replace V7 in manuscript denominators,
formal figures, CNV truth or downstream version claims.
