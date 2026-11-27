Opening atlas figure renderer

Consumes an external frozen display package, source registry and resource configuration; outputs a five-panel English atlas figure and companion identity key. Run render_figure1.py --help for required paths. No coordinate arrays, registry, expression tables or manuscript text are distributed here. Existing source directories are read-only; the output directory must be new.

Requires NumPy, pandas and matplotlib (tested with matplotlib 3.10.8). Set MPLCONFIGDIR to a writable task directory. Keeps scatter layers rasterized and labels/diagram/dotplot as vector elements. Exports plotted aggregate source tables and hash/label receipts into the private output. No embedding refit, downsampling, new differential expression or clinical analysis.
