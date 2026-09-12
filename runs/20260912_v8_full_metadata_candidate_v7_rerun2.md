# Run 20260912_v8_full_metadata_candidate_v7_rerun2

## Metadata

- Run ID: `20260912_v8_full_metadata_candidate_v7_rerun2`
- Run type: full-metadata software candidate / rerun2
- Start time: 2026-09-11T16:35:55.947777+00:00
- End time: 2026-09-11T18:16:03.047490+00:00
- Timezone: UTC, as recorded in the external run summary
- Status: `PASS_FULL_METADATA_SOFTWARE_CANDIDATE_V7`; release status `NOT_FROZEN`
- Exit status: not separately recorded; external summary status is PASS
- Host and scheduler job ID: not recorded in the inspected summary
- Code repository, branch and commit: not recorded in the inspected summary
- Code worktree status: not recorded in the inspected summary
- Environment: Python 3.10.19; pandas 2.3.3; numpy 1.26.4

## Objective

Produce a reproducible, full-cohort V8 metadata software candidate from the frozen
V7 checkpoint while preserving V7 contracts and applying policy gates. This is not
a biological V8 reference run.

## Inputs

| Input ID | Logical path | Version/checksum | Role |
| --- | --- | --- | --- |
| V7 frozen checkpoint | `${V7_RUN_ROOT}/checkpoints/final_annotation_v7.pkl.gz` | `e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e` | Immutable V7 backbone |
| V7 cell index | `${V7_RUN_ROOT}/cell_index` | `3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16` | Frozen index contract |
| Run policy v4 | `${PROJECT_ROOT}/Project_v3/v8_upgrade/config/run_policy_v4.json` | `ae3df80b5a963ae1bc3bc21cd73b3d711cc5edea19c98dda12a93956e2049462` | Candidate-only guardrails |
| Scalable pilot replay attestation | `${PROJECT_ROOT}/Project_v3/runs/20260912_v8_full_metadata_gate_v1/scalable_pilot_replay_attestation_v1.json` | `2d676686d640e0067688c4ac30298df26986fda1d94073ec43e4230a8616031b` | Required precondition gate |

## Input validation

- The external summary records one V7 load, the frozen checkpoint SHA-256 and
  frozen cell-index SHA-256 above.
- The external statistics audit reports matching candidate dimensions and index
  hash, all recorded stages passing, and 128/128 frozen V7 field contracts passing.

## Parameters and policy boundary

| Parameter | Value | Source/rationale |
| --- | --- | --- |
| default project reference | V7 | run summary and policy v4 |
| full candidate only | true | policy v4 |
| full expansion allowed | false | policy v4 |
| automatic reference promotion | false | policy v4 |
| formal downstream mutation | false | policy v4 |
| independent accuracy claimed | false | run summary and policy v4 |
| new classifier execution claimed | false | run summary and policy v4 |
| overwrite allowed | false | run summary and policy v4 |

## Commands

The exact execution command line was not retained in the inspected summary and is
not reconstructed here. Policy v4 identifies the authorized driver module as
`v8_upgrade.integration.run_migration_v7`; its use is bounded by the recorded
policy and pilot-attestation gates.

## Outputs

| Artifact ID | Logical path | Size/checksum | Status |
| --- | --- | --- | --- |
| Full metadata candidate | `${PROJECT_ROOT}/Project_v3/runs/20260912_v8_full_metadata_candidate_v7_rerun2/annotation_v8_full_candidate_v7.pkl.gz` | 251,914,833 bytes; `23fbf7e21768687bd4df89558462db04be710afd2b9a5a7c37db70e6ed69567c` | external; `NOT_FROZEN` |
| Run summary | `${PROJECT_ROOT}/Project_v3/runs/20260912_v8_full_metadata_candidate_v7_rerun2/run_summary_v7.json` | 140,529 bytes; `4dfa482c82adcc28356d0aecbc1f1782be395be3f812b1717e89798d678d998f` | external; PASS summary |
| Statistics audit | `${PROJECT_ROOT}/Project_v3/runs/20260912_v8_full_metadata_candidate_v7_rerun2/full_candidate_statistics_v1.json` | 112,577 bytes; `bf185a25b485b11a0d0f0b9ccb366e8cf19b22de1c46675f8f6745b773a49e49` | external; PASS audit |
| Audit report | `${PROJECT_ROOT}/Project_v3/v8_upgrade/reports/08_FULL_METADATA_CANDIDATE_AUDIT_CN.md` | 19,021 bytes; `2a5f868eba8071e71231567a5ad70b9b166bdf17de1e6e6d9b2557ebf2108cc5` | external; scope/limitations plus append-only v2 correction |

The complete external index, including the scATOMIC salvage, policy and
attestation, is [the V8 artifact manifest](../releases/annotation_v8_development_20260912/external_artifacts.tsv).

## Output validation

- Candidate dimensions match the V7 universe: 4,676,787 cells, 1,322 samples, 54
  GSE and 43 cancers.
- The reported fresh gzip round-trip, prewrite/postwrite schema checks, V7 index
  hash, candidate SHA-256 and 128 field contracts all passed.
- On 2026-09-12 this record-repository task directly recomputed the listed
  candidate, summary, statistics, report, policy, attestation and salvage-file
  SHA-256 values; they match the recorded evidence values.

## Key results

- The output is an auditable full-cohort metadata software candidate with 171
  columns after collision-safe sidecar merging; the 128 frozen V7 fields are
  preserved under the recorded bridge contract.
- The supporting scATOMIC retry2 salvage is a `PASS` exact-order, read-only
  postprocess result. It did not rerun scATOMIC and does not establish accuracy or
  scale-up capacity.
- The final reproducible scATOMIC supportive audit v2 is `PASS_SUPPORTIVE_ONLY`,
  scope-unverified and has zero candidate-pickle deserializations. The external
  audit report's append-only v2 correction was rehashed in this index; it does not
  alter the candidate pickle or its recorded summary/statistics hashes.

## Warnings and interpretation boundary

- `NOT_FROZEN` is intentional. Schema/integrity coverage is not biological
  accuracy, independent truth, new classifier evidence or a V8 clinical/scientific
  completion claim.
- V7 remains the default project reference. Existing downstream artifacts must not
  be overwritten, promoted or labeled V8-derived from this run.
- The candidate has no qualified independent accuracy claim or formal downstream
  mutation authorization.

## Errors

No rerun2 execution error is recorded in the inspected summary. The related
scATOMIC retry2 historical failure remains preserved outside this candidate; the
new supporting salvage does not rewrite that failure history.

## Reproducibility notes

The external policy, attestation and checksums are indexed in the V8 release. The
source repository/commit, worktree state, host, scheduler ID and exact command are
not recorded in the inspected summary, so no values are inferred here.

## Next action

Create a new versioned biological-validation run only after independent reference
or truth evidence, classifier/orthogonal evidence where applicable, and a formal
downstream-impact review are available. A promotion decision must be explicit and
append-only.
