"""Compare CMs found by each NMF method against the frozen V7 nsNMF reference.

Usage: python -I compare_methods.py <run_dir> <reference_W.tsv> <reference_H.tsv>
Each factor is matched to a reference CM by Hungarian assignment on node-membership cosine
(W rows normalised to sum 1, so per-node scale and the X-vs-counts input arms are comparable).
Writes <run_dir>/compare/*.tsv and summary.md.
"""
import itertools, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment

run, refW, refH = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
out = run / "compare"; out.mkdir(exist_ok=True)
meta = pd.read_csv(run / "inputs/metadata.tsv", sep="\t", index_col=0)
nodes = pd.read_csv(run / "inputs/nodes.tsv", sep="\t", index_col=0).majorCluster
TOP = 10

def load(d):
    W = pd.read_csv(d / "W.tsv", sep="\t", index_col=0).loc[nodes.index]
    H = pd.read_csv(d / "H.tsv", sep="\t", index_col=0)[meta.index]
    return W, H, json.load(open(d / "info.json")) if (d / "info.json").exists() else {}

def membership(W):  # node -> factor share; unit-norm columns for cosine
    M = W.div(W.sum(1).replace(0, np.nan), axis=0).fillna(0)
    return M / np.linalg.norm(M, axis=0).clip(1e-300)

def top(W, raw):  # raw W for X-input methods (CoVarNet definition), membership for count inputs
    S = W if raw else W.div(W.sum(1).replace(0, np.nan), axis=0).fillna(0)
    return {c: set(S[c].nlargest(TOP).index) for c in S}

def hcorr(A, B):
    a = (A.T - A.T.mean()) / A.T.std(); b = (B.T - B.T.mean()) / B.T.std()
    return (a.T @ b) / (len(a) - 1)

def match(Wa, Wb):
    C = membership(Wa).T @ membership(Wb)
    r, c = linear_sum_assignment(-C.to_numpy())
    return C, list(zip(Wa.columns[r], Wb.columns[c], C.to_numpy()[r, c]))

def eta2(H):  # share of each factor's sample variance explained by source
    h = H.T.div(H.T.std().replace(0, np.nan), axis=1).fillna(0)
    tot = ((h - h.mean()) ** 2).sum()
    btw = h.groupby(meta.source).apply(lambda g: len(g) * (g.mean() - h.mean()) ** 2).sum()
    return btw / tot

RW = pd.read_csv(refW, sep="\t", index_col=0).loc[nodes.index]
RH = pd.read_csv(refH, sep="\t", index_col=0)[meta.index]
rtop = top(RW, True)
runs = {m.name: {d.name: load(d) for d in sorted(m.glob("seed_*")) if (d / "W.tsv").exists()}
        for m in sorted((run / "methods").iterdir())}
runs = {m: s for m, s in runs.items() if s}
rows, novel, stab, eta = [], [], [], []
for m, seeds in runs.items():
    for s, (W, H, info) in seeds.items():
        raw = info.get("input", "X_minmax").startswith("X_minmax")
        C, pairs = match(RW, W); hc = hcorr(RH, H); mt = top(W, raw)
        for rc, f, cos in pairs:
            rows.append(dict(method=m, seed=s, ref_CM=rc, factor=f, member_cos=cos, H_r=hc.loc[rc, f],
                             top10_jaccard=len(rtop[rc] & mt[f]) / len(rtop[rc] | mt[f]),
                             top10=";".join(W[f].sort_values(ascending=False).index[:TOP]) if raw else ";".join(sorted(mt[f], key=lambda n: -W.loc[n, f]))))
        best = C.max(0)
        for f in W.columns[best.to_numpy() < 0.5]:
            novel.append(dict(method=m, seed=s, factor=f, best_ref=C[f].idxmax(), best_cos=best[f], top10=";".join(sorted(mt[f], key=lambda n: -W.loc[n, f]))))
        e = eta2(H); eta.append(dict(method=m, seed=s, median_source_eta2=e.median(), max_source_eta2=e.max(),
                                     seconds=info.get("seconds"), input=info.get("input")))
    for (s1, a), (s2, b) in itertools.combinations(seeds.items(), 2):
        _, p = match(a[0], b[0]); stab.append(dict(method=m, seed_a=s1, seed_b=s2, mean_member_cos=np.mean([x[2] for x in p]), min_member_cos=np.min([x[2] for x in p])))
rows, novel, stab, eta = map(pd.DataFrame, (rows, novel, stab, eta))
for n, d in [("ref_cm_matches", rows), ("unmatched_factors", novel), ("seed_stability", stab), ("source_eta2_runtime", eta)]:
    d.to_csv(out / f"{n}.tsv", sep="\t", index=False)
rec = rows.pivot_table(index="ref_CM", columns="method", values="member_cos", aggfunc="mean").round(2)
jac = rows.pivot_table(index="ref_CM", columns="method", values="top10_jaccard", aggfunc="mean").round(2)
hr = rows.pivot_table(index="ref_CM", columns="method", values="H_r", aggfunc="mean").round(2)
summ = pd.DataFrame({"seeds": rows.groupby("method").seed.nunique(),
                     "mean_ref_member_cos": rows.groupby("method").member_cos.mean(),
                     "mean_ref_H_r": rows.groupby("method").H_r.mean(),
                     "mean_ref_top10_jaccard": rows.groupby("method").top10_jaccard.mean(),
                     "seed_stability_cos": stab.groupby("method").mean_member_cos.mean() if len(stab) else np.nan,
                     "median_source_eta2": eta.groupby("method").median_source_eta2.median(),
                     "unmatched_per_seed": novel.groupby("method").size().reindex(rows.method.unique(), fill_value=0) / rows.groupby("method").seed.nunique(),
                     "median_seconds": eta.groupby("method").seconds.median()}).round(3)
# null: node labels permuted in each method's first seed (100x) -> chance-level matching to the reference
rng = np.random.default_rng(0); null = {}
for m, seeds in runs.items():
    W = next(iter(seeds.values()))[0]; v = []
    for _ in range(100):
        P = pd.DataFrame(W.to_numpy()[rng.permutation(len(W))], index=W.index, columns=W.columns)
        v.append(np.mean([p[2] for p in match(RW, P)[1]]))
    null[m] = (np.mean(v), np.quantile(v, 0.95))
summ["null_member_cos_mean"] = pd.Series({m: v[0] for m, v in null.items()}).round(3)
summ["null_member_cos_q95"] = pd.Series({m: v[1] for m, v in null.items()}).round(3)
labels = pd.Series({c: " + ".join(RW[c].nlargest(4).index) for c in RW}, name="reference top-4 nodes (raw W)")
# cross-method agreement: mean Hungarian member cosine over all seed pairs
ms = list(runs); X = pd.DataFrame(np.nan, ms, ms)
for a, b in itertools.product(ms, ms):
    v = [np.mean([p[2] for p in match(wa[0], wb[0])[1]]) for (sa, wa), (sb, wb) in itertools.product(runs[a].items(), runs[b].items()) if not (a == b and sa == sb)]
    X.loc[a, b] = np.mean(v) if v else np.nan
X.round(3).to_csv(out / "cross_method_member_cos.tsv", sep="\t")
summ.to_csv(out / "method_summary.tsv", sep="\t")
with open(out / "summary.md", "w") as f:
    f.write("# CM benchmark summary\n\n## Method summary\n\n" + summ.to_markdown() +
            "\n\n## Reference CMs\n\n" + labels.to_frame().to_markdown() +
            "\n\n## Reference CM recovery (membership cosine, mean over seeds)\n\n" + rec.to_markdown() +
            "\n\n## Reference CM top-10 Jaccard\n\n" + jac.to_markdown() +
            "\n\n## Reference CM H Pearson r (matched factor)\n\n" + hr.to_markdown() +
            "\n\n## Cross-method agreement\n\n" + X.round(2).to_markdown() + "\n")
print(open(out / "summary.md").read())
