# Initial manuscript code publication and integration benchmark audit

Date: 2026-11-21. Times below are UTC. This record contains technical provenance and scope, with no scientific measurements, private input maps or individual records.

## Task trajectory

- User requested publication of the initial manuscript plotting scripts to `XX-114514/Pancancer_code`, with a README. The accepted November 20 draft snapshot was kept separate from subsequent revisions.
- User additionally requested a server search for existing scIB comparisons and an explanation of integration-method choice. Existing scripts, aggregate tables, terminal summaries and the current atlas coordinate provenance were reviewed read-only. No integration, annotation or scientific model was rerun.
- Packaging checks completed: Python parsing and command help for every source script; frozen-source SHA256 checks; documentation links; matching external input-manifest verification. The original full rendering acceptance remains separate from these publication checks.
- Command-line Git push returned an authentication-related HTTP 403. The connected GitHub interface successfully published the reviewed source files and documentation.
- Automatic review rejected publication of the internal figure-input map because it exposed private input paths, filenames and coverage metadata. The public map was replaced with a generic renderer/access contract. The complete input map remains in the private companion delivery. The rejected content was not published through an alternative route.
- Public code release verified by an independent clone and comparison of every frozen source hash. Cloud head: `5098c7733ac2fcd6526bfe3f8ef788b01e2b2c4a`. Repository: https://github.com/XX-114514/Pancancer_code
- By 11:06 UTC, historical aggregate score arithmetic had been checked, comparative figures rendered and visually inspected, and a private interpretation report assembled. Historical scIB comparisons, current neighborhood diagnostics and state-preservation diagnostics were kept distinct. Representation labels, scaling modes and benchmark scopes were not conflated.

## Boundaries and next use

The public code repository contains scripts, dependencies, bilingual reproduction instructions, figure mapping and code provenance. Scientific input tables, arrays, manuscript content and learned weights remain outside public Git. The benchmark review bundle is a private server artifact; its recorded scores are not a newly fitted benchmark or independent validation of the final annotation subset. Method-choice prose must match the population and metric actually evaluated.
