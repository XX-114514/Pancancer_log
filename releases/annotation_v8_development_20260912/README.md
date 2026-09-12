# PanCancer V8 development software-candidate index

- Release ID: `annotation_v8_development_20260912`
- Candidate run: `20260912_v8_full_metadata_candidate_v7_rerun2`
- Status: `NOT_FROZEN`
- Scope: full-cohort metadata software candidate, not a biological annotation release
- Current project reference: [V7 frozen identity/CNV release](../annotation_v7_20260905/README.md)

This is a lightweight, append-only evidence index. It records externally retained
Project_v3 artifacts by logical path, byte size, SHA-256 and observed status; it
does not copy the 251,914,833-byte candidate pickle, cell-level data, full runtime
logs or Project_v3 source files into Git.

## Observed candidate status

The rerun2 summary records `PASS_FULL_METADATA_SOFTWARE_CANDIDATE_V7`. The
candidate has the frozen V7 universe (4,676,787 cells, 1,322 samples, 54 GSE and
43 cancers), the V7 index hash, and 128/128 frozen V7 field contracts passing.
The external statistics audit records `PASS` for the candidate's internal
integrity, schema and fresh gzip round-trip checks.

These are software-contract and coverage results. They are not independent
identity or malignancy accuracy, a new classifier execution, an independent truth
assessment, or a biological-completion claim.

## Non-promotion boundary

`run_policy_v4.json` and the run summary explicitly retain V7 as the default
project reference. The candidate is `full_candidate_only`; automatic reference
promotion and formal downstream mutation are disabled. Therefore V7 remains
frozen, and existing V5/V7-derived downstream analyses must not be overwritten or
described as V8-derived.

The scATOMIC retry2 evidence is a supporting, single exact-sample, read-only
postprocess salvage. Its manifest records no algorithm rerun, network, install,
overwrite or patient merge. It is not an accuracy or scale-up result. The canonical
supportive comparison is the final reproducible `full_candidate_scatomic_lung_t20_supportive_audit_v2.json`,
which is `PASS_SUPPORTIVE_ONLY`, scope-unverified and performs no candidate-pickle
deserialization. Its v1 predecessor remains external historical context rather
than the canonical index target.

## External evidence

See [external_artifacts.tsv](external_artifacts.tsv) for the canonical logical
paths, sizes, directly re-computed SHA-256 values and status labels. The candidate
run record is [here](../../runs/20260912_v8_full_metadata_candidate_v7_rerun2.md).
The release metadata is [VERSION.json](VERSION.json), and this index's data-policy
review is [DATA_REVIEW.md](DATA_REVIEW.md).

## Next gate

Any V8 biological-reference decision requires a new, versioned validation run with
independent truth/reference or other qualified evidence, a review of downstream
impact, and an explicit promotion decision. It cannot be inferred from this
software candidate.
