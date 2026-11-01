# CT75 finding: Cloudflare synchronization

## Request and starting state

User requested publishing the current finding from the research server to Cloudflare and an explanation of synchronization. Scientific delivery version: 2026-11-01 qualified TME CT75 v1, complete before publication. Historical CM02 clinical evidence belongs to its separate version. Preserve unrelated dirty worktree files.

## Plan and evidence

- Read repository instructions, current status/roadmap/todo, existing CT75 deliverable and source manifests; inspect Cloudflare inventory.
- Existing Cloudflare login refreshed using Wrangler. Read-only account inventory found no Pages projects and no Workers scripts; a separate Pages project is appropriate.
- Connect through the configured Tailscale jump; read-only verification completed with actual exit 0 at 2026-11-01T22:18:46 UTC / 2026-11-02T06:18:46 UTC+8. All seven server export manifest SHA256 values match the local frozen review and the server asset bytes.
- Stage only reviewed figures and aggregate tables. Prepared 24 main PNG figures, nine SVG counterparts, eleven aggregate TSV tables, original report HTML, publication index/sync explanation, README, bounded ZIP and SHA256 manifest: 50 files. No patient/sample/cell rows, raw matrices, learned parameters or credentials included.
- Frozen original HTML and scientific assets are unchanged. The publication index adds synchronization/download links; explanatory copy distinguishes the current CT75 finding from historical CM02 and describes RNA compatibility at its actual evidence level.

## Validation and pre-deployment checkpoint (historical)

Local asset hashes and all local HTML links pass. Server read-only result, private transport receipt, preparation script and bounded publication acceptance are kept in the task workspace. Initial sandbox Wrangler inventory failed to refresh authentication because of network/write restrictions; authorized ordinary execution restored login. No research computation was restarted.

Cloudflare project creation/deployment and remote HTTP/hash verification are in progress. This is a fixed-version manual publication, not a recurring background sync.

## Next step

Complete publication, record production URL/deployment ID and verify each served file against the prepared SHA256 manifest. Publish this narrowly scoped technical checkpoint without staging unrelated local history.

Suggested commit: docs: record CT75 Cloudflare publication checkpoint

## Completed publication and independent verification

- Production URL: https://pancancer-ct75-review.pages.dev/
- Immutable deployment URL: https://4436c5b3.pancancer-ct75-review.pages.dev
- Production deployment ID: `4436c5b3-9b1f-40f5-9b70-2b3176a88e88`.
- Cloudflare reports production deployment success at 2026-11-01T22:20:17.909976Z UTC / 2026-11-02T06:20:17.909976 UTC+8 / 2026-11-02T09:20:17.909976 Melbourne.
- All 50 served files independently returned HTTP 200 and matched exact prepared bytes and SHA256, including the 9,842,221-byte download ZIP. The ZIP SHA256 is `5c1e9fc8981b2f9fb2dcfdfa7c191c7a10ec71db8d4685a3804b996b08dccfc8`; public manifest SHA256 is `08ee8d9c8e8895bf7d21d5f80e166313ed615c0356f2605601789c10a0207ea6`.
- The live HTTPS page was opened and visually inspected; title, scope counts and synchronization/download entries display correctly. Screenshot retained in the private delivery workspace. This verifies the newly published public page; the earlier blocked local-file browser attempt remains separate historical evidence.
- Wrangler initially delegated new-project creation to a Worker deployment without an entry point. Nothing was deployed. A direct Pages project creation succeeded, then ordinary Pages upload completed; no existing project was changed.
- Repository validation completed with exit 0, zero errors and one inherited warning. Task diff whitespace check passes. Python source compilation passed using a task-local temporary cache after the default system cache directory was not writable. A preliminary diff command used the parent non-repository directory; corrected invocation in the repository passed.

This publication is a fixed-version server-to-verified-local-staging-to-Cloudflare snapshot. No recurring task was configured. Future server changes require a new explicit validated publication. Scientific fits, authoritative labels, raw data and unrelated repository changes remain untouched. This task log is the narrowly scoped public synchronization record; larger private execution histories are not included in this GitHub checkpoint.
