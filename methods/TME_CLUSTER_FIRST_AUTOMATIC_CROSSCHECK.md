# Cluster-first working names and automatic reference crosscheck

## Purpose

Name the accepted numeric graph clusters from their own native expression
evidence, then compare those names with scope-appropriate automatic references.
Keep broad lineage names where a finer identity is not resolved. Cluster-level
working names and individual-cell identity evidence remain distinct fields.

## Marker-first names

Bind each name to the exact cluster membership and retain that membership SHA in
a separate field. Use an informative marker from the cluster's own accepted
statistics: positive mean-expression enrichment and detection enrichment over
the remaining fitted cells. When an accepted inferential DEG screen is
available on the same membership, distinguish it from descriptive expression
specificity. Borrowing a parent marker requires an explicit provenance field.

The display format is `Celltype_MARKER_number`. A single marker is for display;
the lineage interpretation uses a multigene expression pattern. Preserve the
numeric graph cluster, display number and annotation version separately.

## Official CellTypist crosscheck

Use the installed official prediction interface with cached, fingerprinted
models. Normalize the complete adopted native RNA counts to CP10K and log1p
before reference matching. Retain accepted ambient subtraction. Unresolved
native gene symbols keep their RNA contribution to the denominator, with
explicit identifiers that cannot match a reference feature. Unmeasured model
features are reported as coverage gaps rather than measured negative genes.

Run predictions in bounded blocks and retain private per-cell proposals. Reuse
the existing graph clusters for official majority voting. Summarize coarse
families separately from fine reference labels; retain heterogeneous votes and
ties as review cases. Reference probabilities are classifier outputs, not
calibrated identity accuracy.

Immune references diagnose immune populations. They cannot directly assign an
epithelial or stromal identity. Organ and compartment references follow a
frozen applicability scope. Their predictions remain available for discrepancy
review. Healthy-organ similarity does not determine tumor malignancy.

## Acceptance

Independently recompute native axes, exact group membership, all source and
cluster counts, model coverage masks and majority votes. Retain cells outside
the fit pool with their baseline annotations and an explicit absent cluster.
Preserve previous identity, state, malignancy, QC and source-qualification
fields. A completed reference crosscheck is a separate milestone from final
CM-input eligibility.

Sources: [official prediction API](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html),
[official guide](https://www.celltypist.org/tutorials/onlineguide),
[official implementation](https://github.com/Teichlab/celltypist).
