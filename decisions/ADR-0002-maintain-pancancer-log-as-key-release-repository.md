# ADR-0002: maintain Pancancer_log as the key lightweight release repository

- Date: 2026-09-06
- Status: accepted

## Context

The PanCancer analysis parent contains more than one terabyte of data across
Project_v1-v5, several nested Git repositories and runtime environments. A separate
`project-records` repository already preserves the history and private remote of
`Pancancer_log`.

## Decision

Continue `project-records` as the single maintained Git repository. Integrate
current progress through versioned lightweight releases and external-object
manifests instead of initializing a competing parent-level monorepo.

## Consequences

- Existing history and remote remain intact.
- Heavy and sensitive objects cannot be accidentally staged with `git add .`.
- Prism and supplementary materials have stable, reviewable inputs.
- Reproduction of heavy analyses still requires access to external objects whose
  logical paths and hashes are recorded.
- V7 identity/CNV and V5-derived downstream communication results remain explicitly
  version-separated.
