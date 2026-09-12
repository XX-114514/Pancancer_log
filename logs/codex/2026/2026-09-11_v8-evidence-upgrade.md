# V8 evidence-system upgrade

## User request

Actually upgrade the V7 annotation system to an auditable V8, retaining V7 as the
conservative backbone; add separated identity/state/malignancy schema, classifier
evidence workflows, independent benchmark safeguards and orthogonal-evidence gates.

## Initial observed state

Record repository clean at `0e6f3c9`. V7 is the frozen project reference;
communication release is V5-derived. No existing V8 annotation upgrade was found
in the scoped initial filename check. Full baseline and dependency audit is at
`${PROJECT_ROOT}/Project_v3/v8_upgrade/reports/00_baseline_and_dependencies.md`.

## Plan and ownership

Root owns shared contracts, integration and final audit. Three Luna MAX workers
own disjoint identity, benchmark and malignancy directories. Read-only audits
precede implementation. Pilot precedes expansion. Scientific validation is a
separate gate from passing software tests.

## Inspected evidence

Repository rules; README/STATUS/ROADMAP/TODO; version/data policies; V7 method and
run record; previous Git-release task log; V7 freeze/reference/config/full method;
annotation, CNV, pipeline and resource skills; official method documentation.

## Commands and initial results

- `git status --short`: clean; `git log -3 --oneline`: latest `0e6f3c9`.
- `df -h`, `quota -gs`, `free -h`, GPU utilization query: sufficient for bounded
  metadata/pilot work; existing GPU utilization prevents assuming exclusive use.
- Resource detector wrote a new local JSON snapshot; no old snapshot overwritten.
- Default sandbox execution failed due to host user-namespace support; scoped
  escalated read commands succeeded.
- Additive patch created V8 README and baseline report.

## Current validation and limitations

Implementation and biological validation are in progress, not yet accepted.
No new classifier output, independent accuracy, Numbat run or paired-DNA
validation is claimed at this stage. V7 and formal downstream objects unchanged.

## Suggested commit message

`feat: add audited V8 annotation evidence upgrade candidate`
