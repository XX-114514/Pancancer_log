# English homepage and six companion pages

## Request and starting state

The user requested an English version of the repository homepage README and its corresponding subpages. The source is the public 2026-11-02 paper-style homepage, based on commit `e1900f8` and the 0906 analysis release. The existing English README still described the 2026-10-26 CM02 immune atlas.

An existing local checkout contained substantial unrelated work. A separate shallow checkout of the confirmed public repository was created for this task; its initial working tree was clean. The original checkout was not modified.

## Plan and inspected material

Translate the current homepage and six numbered `docs/homepage/` pages, preserve the previous English homepage, add reciprocal language navigation, and verify content and links before ordinary commit/push.

Inspected: `AGENTS.md`, `README.md`, `README.en.md`, the six Chinese companion pages, current `STATUS.md`, `ROADMAP.md`, and `TODO.md` sections, `docs/CODEX_CLOUD_DEVELOPMENT.md`, `docs/RUNTIME_SYNC.md`, `docs/MANUSCRIPT_NARRATIVE_PRINCIPLES.md`, the 2026-11-02 homepage publication task log, the related daily comparison/publication record, and `scripts/validate_repository.sh`. The homepage source log records the relevant 0906 release, annotation, and figure provenance. No scientific execution was needed.

## Changes

- Replaced `README.en.md` with the English counterpart of the current homepage, including all seven figure captions, ten main sections, references, and version boundaries.
- Added English counterparts of all six numbered companion pages, retaining tables, captions, parameters, and limitations.
- Preserved the former English homepage as `README_20261026_CM02.en.md`; only its navigation line differs from the previous file.
- Updated the Chinese homepage and six Chinese pages with English navigation; corrected the historical Chinese homepage to link to its corresponding English archive.
- Kept the original figures and scientific results. Supporting historical methods/reports remain in their source language, explicitly noted on the English homepage.

## Checks and findings

- All relative links and heading anchors in the seven translated pages resolve in the checkout; image counts match the Chinese originals (7 homepage images and 8 companion-page images, including two reused supplements).
- The archived English body matches the previous `README.en.md` except for its navigation line.
- Numerical comparison and manual review retain cohort sizes, thresholds, effect sizes, intervals, and version identifiers. Chinese units are converted explicitly: 1,121万 → 11.21 million and 38.2万 → 382,000.
- `git diff --check` passed. Reviewed the Chinese-page diffs: changes are confined to language/archive navigation.
- Full repository validation completed with 1 pre-existing error and 1 pre-existing warning. The error is the link from `docs/MALIGNANCY_TME_V1.md` to absent `inventories/tme_v1/`; both the link and missing target were verified against the untouched base commit. The warning lists existing internal-path patterns in unrelated files. No translated or newly added page is implicated. Secret-pattern, size, data-format, and TSV checks passed. The unrelated broken link was not changed.

## Commands and failed attempts

Used read-only GitHub fetches, `git clone --depth 1`, `git status --short`, `git diff`, a local Python link/anchor/archive check, and `bash scripts/validate_repository.sh`. The first clone attempt failed because DNS/network access was unavailable in the default sandbox; the approved network-enabled retry succeeded. An optional process-status inspection was unavailable in the sandbox and did not affect document checks.

## Scope and remaining work

This is a translation and documentation-navigation change. Scientific analyses, source-paper claims, external websites, and downstream clinical results were not recomputed or independently revalidated. No patient records, expression matrices, or new scientific data were added. Historical methods and reports linked as supporting references are outside the seven-page translation scope.

Next: commit only these documentation changes and verify ordinary publication against the remote branch. The unrelated baseline broken link remains open. Suggested commit message: `docs: translate current homepage and six companion pages into English`.

## Time record

Documentation review evidence timestamp: 2026-11-02 22:29:20 UTC / 2026-11-03 06:29:20 UTC+8 / 2026-11-03 09:29:20 Australia/Melbourne. The earlier exact task start time was not recorded; the dates above are review timestamps, not scientific run times.
