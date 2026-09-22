# Data and privacy review

This release contains six small Markdown/JSON source snapshots plus repository
documentation and a checksum index.

- No cell-level tables, expression matrices, H5AD/RDS/pickle/Zarr objects, raw
  sequencing data or complete runtime logs are included.
- The copied files contain relative Project_v3 paths, aggregate GSE/sample/cell
  counts and hashes; no literal internal absolute server path was detected.
- Targeted scans found no credential-like strings or direct patient identifiers.
- Diagnostic JSONs, audit scripts, large per-cell outputs and path-bearing reports
  remain external and are not copied into Git.

Private GitHub storage does not change these minimum-disclosure rules.
