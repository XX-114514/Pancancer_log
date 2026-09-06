# Git maintenance and release workflow

## Repository role

The maintained repository is `${PANCANCER_ROOT}/project-records`, preserving the
history of the existing `Pancancer_log` remote. It is the project control plane
for status, methods, provenance, lightweight releases and manuscript assets. The
adjacent Project_v1-v5 analysis trees are not Git payloads.

## Storage contract

Track:

- Markdown methods, decisions, run records and task logs;
- code and configuration after portability/privacy review;
- aggregate TSV/CSV/JSON tables;
- publication PDF/PNG figures;
- checksums and logical paths for external objects.

Do not track raw human data, cell/patient-level matrices, H5AD/RDS/Zarr/pickle
objects, complete runtime environments, scheduler output, credentials or files
larger than 5 MiB. Git LFS is not installed on the current host, so large objects
must remain external.

## Normal update

```bash
cd "${PANCANCER_ROOT}/project-records"
git status --short
bash scripts/validate_repository.sh
git diff --stat
git diff
git add <explicitly-reviewed-files>
git diff --cached --stat
git diff --cached
git commit -m "type: concise description"
git push origin main
```

Push only after confirming `origin` is the intended private repository and
authentication is available. Never use force push or rewrite frozen history.

## Release rule

A frozen release is append-only. Do not overwrite a V7 file after publication.
Any biological/statistical change creates a new run and version (V8 or later).
A corrected narrative that does not change data may be added as a dated erratum.

Each release should contain:

1. README with scope and version boundaries;
2. machine-readable version metadata;
3. small reports, source tables and figures;
4. external artifact manifest with hashes;
5. source-code provenance;
6. repository validation evidence.

## Before Prism or supplementary export

Use [the Prism and supplementary guide](PRISM_SUPPLEMENTARY_GUIDE.md), verify that
all plotted statistics match the source table, and cite the exact release ID.
Do not mix V7 annotation statistics with V5-derived CoVarNet/LIANA results without
an explicit version statement.
