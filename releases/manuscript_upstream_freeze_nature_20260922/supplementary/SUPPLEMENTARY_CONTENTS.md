# Supplementary and Extended Data plan

Nature prefers integral small figures/tables as Extended Data. Supplementary
Information should contain only essential material that is too large, specialized
or impractical for the main article/Extended Data.

## Extended Data candidates

Use `figures/FIGURE_INVENTORY.tsv` as the controlling list. Keep the combined number
of Extended Data figures and tables at or below ten unless the editor instructs
otherwise. Each item must be interpretable independently and have exact source data.

## Main Supplementary Information PDF

Provisional order:

1. Supplementary Methods that do not fit the online Methods.
2. Supplementary Notes on annotation uncertainty and version boundaries.
3. Large study-flow and cohort provenance descriptions.
4. Supplementary References.

Do not duplicate Extended Data or deposit large matrices in the PDF.

## Separate Supplementary Data

Large, machine-readable tables should be separate TSV/XLSX/CSV files with stable
column dictionaries, release IDs and checksums. Raw/cell-level data require a
repository and governance review rather than Git or a giant SI file.

## SIGuide

Before submission create `SIGuide.doc` listing every SI file, its title and a
description of no more than 50 words. Treat SI as final at acceptance because
Nature does not copyedit it and post-acceptance changes are restricted.
