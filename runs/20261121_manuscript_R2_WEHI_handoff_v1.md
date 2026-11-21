# Manuscript R2 to WEHI handoff v1

The existing accepted manuscript archive and a private Claude plotting handoff were transferred through the user's existing R2 bucket. The source archive SHA256 and ZIP CRC were rechecked; the archive and every sidecar were then read back from R2. A separate R2 receipt/object-list read confirmed the cloud payload. No raw expression matrices, credentials or new scientific fits were introduced.

The receiver uses immutable copies, checksums, safe extraction, per-file verification and the plotting-input verifier. It creates a real receiver receipt only after those stages pass; no model fit or figure generation is run by the receiver. Exact archive checksums, personal bucket paths and receiving commands remain in private handoff metadata.

At this checkpoint R2 is ready for the intended receiver. A WEHI receipt has not been observed, so transfer to the bucket is not reported as completed receipt or plotting on WEHI. The user reaffirmed the R2-intermediate method; existing private migration configuration is reused without credential disclosure. Receiver execution is available via the supplied bootstrap command in the intended project directory.

Validation: shell syntax, Python AST, synthetic repeatability/traversal fixtures, archive/manifest hashes, R2 full readback and independent object/receipt reads passed. Existing unrelated repository broken links and a warning remain; current-scope checks pass. Original science files, bucket items and unrelated local changes were preserved.

[Task trace](../logs/codex/2026/2026-11-21_manuscript-R2-WEHI-handoff.md) · [Receiver code](../scripts/manuscript_wehi_handoff_20261121/README.md) · [Technical acceptance](../inventories/publication/2026-11-21/manuscript_wehi_handoff_v1/technical_acceptance.json)
