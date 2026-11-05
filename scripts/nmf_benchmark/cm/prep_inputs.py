"""Build the shared CM-benchmark inputs from the V7 primary_strict CoVarNet inputs.

Usage: python -I prep_inputs.py <covarnet_input_dir> <run_dir>
Writes <run_dir>/inputs/: X_minmax.tsv (nodes x samples, CoVarNet min-max of within-parent
frequency, same formula as freq_normalize), counts.tsv (nodes x samples, raw cell counts),
metadata.tsv, X_minmax.h5ad / counts.h5ad (samples x nodes, for cNMF), manifest.json.
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad

src, run = Path(sys.argv[1]), Path(sys.argv[2])
out = run / "inputs"; out.mkdir(parents=True, exist_ok=True)
fq = pd.read_csv(src / "frequency.tsv", sep="\t", index_col=0)
major = fq.pop("majorCluster")
meta = pd.read_csv(src / "metadata.tsv", sep="\t").set_index("sampleID").loc[fq.columns]
cnt = pd.read_csv(src / "node_sample_counts.tsv.gz", sep="\t", index_col=0)
cnt = cnt.loc[fq.index]  # node_sample_counts also holds the 11 nodes CoVarNet dropped
assert (cnt.pop("majorCluster") == major).all()
cnt = cnt[fq.columns]
raw = fq.to_numpy(float)
lo, span = raw.min(1, keepdims=True), np.ptp(raw, 1, keepdims=True)
assert (span > 0).all()
X = pd.DataFrame((raw - lo) / span, index=fq.index, columns=fq.columns)
X.to_csv(out / "X_minmax.tsv", sep="\t")
cnt.to_csv(out / "counts.tsv", sep="\t")
meta.to_csv(out / "metadata.tsv", sep="\t", index_label="sampleID")
major.rename("majorCluster").to_frame().to_csv(out / "nodes.tsv", sep="\t")
for name, m in [("X_minmax", X), ("counts", cnt.astype(float))]:
    ad.AnnData(m.T.to_numpy(np.float64), obs=pd.DataFrame(index=m.columns), var=pd.DataFrame(index=m.index)).write_h5ad(out / f"{name}.h5ad")
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
json.dump({"source_dir": str(src), "frequency_sha256": sha(src / "frequency.tsv"),
           "n_nodes": int(X.shape[0]), "n_samples": int(X.shape[1]), "n_sources": int(meta.source.nunique()),
           "X_minmax_sha256": sha(out / "X_minmax.tsv"), "counts_sha256": sha(out / "counts.tsv")},
          open(out / "manifest.json", "w"), indent=1)
print(X.shape, cnt.shape, meta.source.nunique())
