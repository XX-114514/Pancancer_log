# Portable manuscript plotting

These renderers take an externally supplied, reviewed plotting package. No private input tables or manuscript text are stored in this directory.

```bash
python code/run_all_figures.py --package-root . --output outputs/replay_new
```

The package must supply `inputs/r4`, `inputs/immune`, `inputs/integrated_evidence`, `inputs/R4_figure_tables` and `inputs/figure_configuration.json`. `run_all_figures.py` executes the R4, composition/lymph-node and story renderers in sequence, then assembles the manuscript figure IDs and records output SHA256 hashes. The output directory must be new. The package figure-input map documents exact inputs and the display-array dictionaries.

The portable adaptation retains frozen numerical tables, coordinates, membership and correction families. Plotting edits adjust labels, layout and figure order. Original raw matrices, identifiers and learned weights are not inputs to this plotting workflow. Statistical fitting is outside this command.

Validated rendering environment: Python 3.10, NumPy 1.26, pandas 2.3, matplotlib 3.10 and NetworkX 3.4. The bilingual document builder and supplementary-PDF packager remain in `../integrated_preprint_20261120/` and consume an external manuscript-content JSON.
