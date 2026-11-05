"""One Python NMF method, one seed, K CMs on the shared CM-benchmark inputs.

Usage: python -I run_py_method.py <run_dir> <method> <seed> <K>
methods: sklearn_fro (cd, Frobenius), sklearn_kl (mu, KL), cNMF_X (cNMF on X_minmax), cNMF_counts.
Writes <run_dir>/methods/<method>/seed_<seed>/{W.tsv (nodes x K), H.tsv (K x samples), info.json}
"""
import json, shutil, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.decomposition import NMF

run, method, seed, K = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
X = pd.read_csv(run / "inputs/X_minmax.tsv", sep="\t", index_col=0)
out = run / "methods" / method / f"seed_{seed}"; out.mkdir(parents=True, exist_ok=False)
t0, info = time.time(), {"method": method, "seed": seed, "K": K}
if method.startswith("sklearn"):
    kw = dict(solver="cd", beta_loss="frobenius") if method == "sklearn_fro" else dict(solver="mu", beta_loss="kullback-leibler")
    best = None
    for r in range(30):  # best-of-30 random restarts, as CoVarNet nsNMF nrun=30
        m = NMF(K, init="random", max_iter=2000, tol=1e-6, random_state=seed * 100 + r, **kw)
        Wt = m.fit_transform(X.to_numpy())  # nodes x K (rows of X are nodes)
        if best is None or m.reconstruction_err_ < best[2]:
            best = (Wt, m.components_, m.reconstruction_err_, m.n_iter_)
    W, H = best[0], best[1]
    info |= {"objective": float(best[2]), "n_iter": int(best[3]), "input": "X_minmax", "package": f"scikit-learn {__import__('sklearn').__version__}", **kw}
elif method.startswith("cNMF"):
    from cnmf import cNMF
    name = "X_minmax" if method == "cNMF_X" else "counts"
    h5 = str(run / f"inputs/{name}.h5ad")
    genes = out / "nodes.txt"; genes.write_text("\n".join(X.index) + "\n")
    work = out / "cnmf"
    c = cNMF(output_dir=str(work), name="cm")
    # X arm: X itself is passed as the "TPM" so no per-sample renormalisation; counts arm: cNMF's own TPM
    c.prepare(counts_fn=h5, tpm_fn=h5 if name == "X_minmax" else None, components=[K], n_iter=200,
              seed=seed, genes_file=str(genes), beta_loss="frobenius")
    c.factorize(worker_i=0, total_workers=1)
    c.combine()
    for dt in (0.1, 0.5, 2.0):  # local-density outlier filter; loosen only if a cluster empties
        try:
            c.consensus(k=K, density_threshold=dt, show_clustering=False)
            break
        except Exception as e:
            print(f"consensus failed at density_threshold={dt}: {e}")
    usage, scores, spectra_tpm, _ = c.load_results(K=K, density_threshold=dt, norm_usage=False)
    W = spectra_tpm.loc[X.index].to_numpy()  # load_results returns spectra as nodes x K, input units
    H = usage.loc[X.columns].to_numpy().T
    info |= {"input": name, "n_iter": 200, "density_threshold": dt, "beta_loss": "frobenius", "package": f"cnmf {__import__('importlib.metadata').metadata.version('cnmf')}"}
else:
    raise SystemExit(f"unknown method {method}")
info["seconds"] = time.time() - t0
cols = [f"F{i + 1:02d}" for i in range(K)]
pd.DataFrame(W, index=X.index, columns=cols).to_csv(out / "W.tsv", sep="\t")
pd.DataFrame(H, index=cols, columns=X.columns).to_csv(out / "H.tsv", sep="\t")
json.dump(info, open(out / "info.json", "w"), indent=1)
print(method, seed, "done", round(info["seconds"], 1), "s")
