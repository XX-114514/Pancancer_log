# PanCancer V8 lineage-taxonomy audit addendum

- Release ID: `annotation_v8_lineage_audit_20260917`
- Evidence dates: 2026-09-16 to 2026-09-17
- Repository snapshot date: 2026-09-22
- Status: `AUDIT_ADDENDUM_NOT_GLOBAL_FREEZE`
- Current global project reference: V7

This append-only release preserves the compact audit layer generated after the
five V8 lineage taxonomy freezes. It does not contain a cell-level V8 annotation
object and does not promote V8 to the global reference.

## Audit conclusion

- Global V8 remains `NOT_FROZEN`; no lineage is `HELDOUT_VALIDATED`.
- Myeloid/DC remains `PILOT_FROZEN`, but all 12 later-designated held-out GSEs had
  outputs before the held-out split was frozen.
- T/NK is `PILOT_FROZEN_HELDOUT_COMPROMISED`; GSE116256 contributed to the frozen
  boundary decision and cannot serve as its held-out evaluation cohort.
- Fibroblast and Endothelial are `PILOT_FROZEN` with clean held-out exposure, but
  still carry clustering/panel limitations recorded in their manifests.
- B/Plasma is `DEVELOPMENT_FROZEN`; 45.4% of development cells fell back to a
  competing lineage and held-out outputs are sealed by policy rather than access
  control.

## Included snapshots

- [Generated V8 status](V8_STATUS_SNAPSHOT_20260916.md)
- [Myeloid/DC manifest](manifests/Myeloid_DC.json)
- [T/NK manifest](manifests/T_NK.json)
- [Fibroblast manifest](manifests/Fibroblast.json)
- [Endothelial manifest](manifests/Endothelial.json)
- [B/Plasma manifest](manifests/B_Plasma.json)
- [External artifact index](external_artifacts.tsv)
- [Data review](DATA_REVIEW.md)

The six copied source artifacts are byte-identical to the external Project_v3
files. Related diagnostic audits remain external and are recorded by logical path,
size and SHA-256 to avoid splitting code/data authority.

## Publication boundary

These records qualify the five lineage-taxonomy freezes; they do not invalidate
the already recorded taxonomy tables, create a unified V8 writeback, establish
independent accuracy or authorize V8-derived downstream analyses.
