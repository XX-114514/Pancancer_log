# Annotation and analysis version policy

## Current authority

V7 is the default cell-identity and malignancy reference from 2026-09-06 onward.
V5/V6 remain historical inputs for already completed downstream analyses.

## Version boundaries

- Identity version: major and fine cell-type definitions.
- Malignancy version: candidate/reference construction, CNV engines and call rule.
- Downstream version: CoVarNet, LIANA, spatial and clinical analyses.
- Manuscript release: curated set of tables, figures and prose.

These axes must be recorded independently. A downstream analysis does not inherit
a newer identity version unless it is actually rerun.

## Increment rule

Create V8 or later when any of these change: cell universe, gene mapping,
candidate/reference definition, marker logic, CNV parameters, malignancy rule or
canonical field semantics. Cosmetic prose or figure-layout changes do not require
a biological annotation version bump, but they do require a dated release update.

## Freeze rule

Frozen artifacts are append-only and verified by hash. Never replace an object
under the same version label. Record corrections in a new run, release or erratum.
