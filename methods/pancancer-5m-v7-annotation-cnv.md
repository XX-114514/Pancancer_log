# PanCancer 5M V7 annotation and CNV method

- Status: frozen project reference
- Run ID: `20260905_v7_cnv_rerun_freeze_v1`
- Full specification: 
  [V7 reference and annotation logic](../releases/annotation_v7_20260905/methods/V7_REFERENCE_AND_ANNOTATION_LOGIC_CN.md)

V7 freezes the uncertainty-resolution V2 identity axis and reruns inferCNVpy plus
Copykat_python per original sample from raw non-negative integer counts.

The broad annotation combines genuine source-provided per-cell labels, full-library
marker scores and CellTypist support. Single negative markers no longer hard-veto
a lineage; coherent competing programs, cross-compartment conflicts and sample-level
doublet evidence determine provisional, boundary or ambiguous states. Structurally
missing genes are masked by GSE, aliases are resolved, and `KIM1` maps to `HAVCR1`.

Subtype annotation first uses lineage-specific expanded panels. Remaining clusters
receive de novo marker descriptions only when effect size, detection difference and
cross-sample recurrence pass the frozen thresholds. Descriptive
`*_unresolved_GENE` labels are not treated as mature cell types.

CNV candidate/reference sets are rebuilt within each original sample using
cancer-specific lineage policies, including CNS, melanocytic, neuroendocrine,
mesenchymal, myeloid and germ-cell malignancies. Skips are not negatives.

Primary malignancy is a per-sample CNV-cluster call. The strict two-method
intersection and sensitive union are parallel evidence definitions rather than
synonyms of the primary cluster call. Exact filters and denominators are frozen in
the linked method document and release tables.
