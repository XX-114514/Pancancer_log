# Working title

**A single-cell atlas of cellular states across 43 human cancers**

Title length: 62 characters including spaces. Final title must remain at or below
75 characters and should avoid abbreviations and unnecessary punctuation.

## Authors and affiliations

[AUTHOR LIST — confirm order, equal-contribution notes, affiliations, present
addresses and corresponding author before submission.]

## Summary paragraph

[TARGET: no more than 200 words; fully referenced; no unexplained abbreviations.
Write 2–3 sentences of broad context, the specific gap, “Here we show” or an
equivalent statement of the main result, and 2–3 sentences explaining how the
work advances the field. Do not finalize until the identity version and
CoVarNet/LIANA rerun decision are closed.]

## Main text

### A pan-cancer single-cell resource

[Report the frozen discovery cohort and inclusion/exclusion logic. The current
global annotation reference comprises 4,676,787 cells from 1,322 samples, 54 GSE
and 43 cancers. Reconcile this universe with the earlier 31-GSE discovery object
before writing a single study-flow denominator.]

### Conservative cell-state annotation

[Describe V7 as the global reference. Separate major identity, fine identity,
confidence, boundary/provisional states and strict eligibility. V8 lineage-only
taxonomy evidence belongs in Extended Data unless a global V8 release is completed.]

### Sample-resolved malignant evidence

[Describe per-sample candidate/reference construction and CNV inference on raw
counts. Report inferCNVpy and Copykat_python separately. Explicitly state that
skipped samples are not negative. Do not generalize the unresolved
`genome='hg20'` parameter until audited.]

### Cellular programmes across cancers

[Insert CoVarNet results only from one declared identity/downstream version. The
existing K=9 communication branch is V5-derived. If rerun on the final annotation,
replace this paragraph and all linked panels under a new release ID.]

### Spatial and clinical validation

[Describe each cohort, endpoint, sample unit, fixed-mapper rule and statistical
model. Distinguish exploratory direction from multiplicity-adjusted evidence.
Retain the known limitations for sort-matched subsets and legacy clinical labels.]

### Discussion

[Lead with the principal biological advance, then compare with prior pan-cancer
atlases. Discuss annotation uncertainty, dataset heterogeneity, non-independent
evidence, CNV non-evaluability and the V5/V7 downstream version boundary. End with
the reusable resource and testable implications without clinical overstatement.]

## Main references

[TARGET: normally no more than 50 main-text references. One publication per
number; include titles. Add Methods-only references in a separate sequence as
required by Nature formatting.]

## Tables

[Main tables, if any. Each requires a short title sentence and editable content.]

## Figure legends

[TARGET: less than 250 words per legend. Begin with a title sentence; define exact
n, biological/technical replicates, statistics, tails, error bars and panels.]

## Methods

### Study design and dataset inclusion

[Source: METHODS_EVIDENCE_MAP.md. State preregistration/held-out boundaries and
dataset roles.]

### Expression-matrix processing and quality control

[Specify count semantics, per-sample QC, doublet handling, gene mapping,
normalization views, integration methods and random seeds.]

### Cell identity annotation

[Describe V7 global annotation. If V8 is not globally released, describe lineage
taxonomy work as a separately bounded development analysis.]

### Malignancy inference

[Describe per-sample raw-count inferCNVpy and Copykat_python inputs, thresholds,
minimum-cell gates, mappings, status semantics and evidence combination.]

### CoVarNet and cell–cell communication

[Use the final declared downstream release only. Do not merge V5- and V7/V8-derived
outputs under a single method label.]

### Spatial and clinical analyses

[List cohorts, endpoints, exclusions, fixed projections and multiplicity control.]

### Statistics and reproducibility

[For every analysis give exact n, unit of analysis, test, tail, covariates,
multiple-testing correction, software/version and seed. Define error bars.]

### Ethics

[For each source cohort: public accession or controlled source, original approval
and consent statements, secondary-use basis, and whether identifiable data were
accessed. Authors must verify this text against source studies.]

### Data availability

[Insert the reviewed statement from submission/DATA_CODE_AVAILABILITY_DRAFT.md.]

### Code availability

[Insert the reviewed statement from submission/DATA_CODE_AVAILABILITY_DRAFT.md.]

## Methods references

[Add Methods-only references.]

## Acknowledgements

[TBD]

## Funding

[TBD — separate statement.]

## Author contributions

[TBD — use author initials and specific contributions; verify all authors agree.]

## Competing interests

[TBD — include an explicit declaration even if none.]

## Additional information

Supplementary Information is available for this paper.

Correspondence and requests for materials should be addressed to [NAME].

## Extended Data legends

[Up to ten multi-panel Extended Data display items.]
