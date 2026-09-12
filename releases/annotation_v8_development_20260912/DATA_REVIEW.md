# V8 development-index data review

Review date: 2026-09-12.

This release contains only a Markdown scope statement, JSON metadata and a TSV
index of logical external paths, byte sizes, SHA-256 values and status labels. It
does not copy candidate pickles, AnnData objects, expression matrices, cell IDs,
sample-level tables, raw scATOMIC output, full logs, credentials or absolute
internal paths from Project_v3.

The release is therefore a lightweight provenance index rather than a data or
biological-results package. Source reports and all heavy artifacts remain outside
Git and are identified only through `external_artifacts.tsv`.
