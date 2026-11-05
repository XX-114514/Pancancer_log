"""Per-panel V7 label transfer to Xenium cells; held-out-sample F1 gates which labels are spatially identifiable.
Exports per-section counts (MatrixMarket), coordinates and labels for blisa."""
import sys, json, hashlib, glob
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scipy.io as sio
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import f1_score, classification_report

run = Path(sys.argv[1]); sections = sys.argv[2:]
F1_MIN, MIN_TRAIN = 0.5, 50
ref = ad.read_h5ad(run / 'reference/v7_spatial_reference_counts.h5ad')
out = run / 'xenium'; out.mkdir(exist_ok=True)

def lognorm(M):
    M = sp.csr_matrix(M, dtype=np.float64); s = np.asarray(M.sum(1)).ravel(); s[s == 0] = 1
    return np.log1p(sp.diags(1e4 / s) @ M).toarray()

by_panel = {}
for f in sections:
    genes = tuple(ad.read_h5ad(f, backed='r').var_names)
    by_panel.setdefault(hashlib.sha1('|'.join(genes).encode()).hexdigest()[:8], []).append(f)
for key, files in by_panel.items():
    files = [f for f in files if not (out / Path(f).name.replace('.counts.h5ad', '') / 'cells.tsv.gz').exists()]
    if not files: continue
    xg = list(ad.read_h5ad(files[0], backed='r').var_names); genes = [g for g in xg if g in set(ref.var_names)]
    r = ref[:, genes]; keep = np.asarray(r.X.sum(1)).ravel() >= 10; r = r[keep]
    tr = r.obs.split.eq('train').to_numpy(); counts = r.obs.label[tr].value_counts()
    labels_ok = counts.index[counts >= MIN_TRAIN]; m = r.obs.label.isin(labels_ok).to_numpy()
    Xr = lognorm(r.X); sc = StandardScaler().fit(Xr[tr & m])
    clf = LogisticRegression(max_iter=1000, C=0.1, n_jobs=-1).fit(sc.transform(Xr[tr & m]), r.obs.label[tr & m])
    ho = (~tr) & m; pred = clf.predict(sc.transform(Xr[ho]))
    f1 = pd.Series(f1_score(r.obs.label[ho], pred, labels=clf.classes_, average=None), index=clf.classes_, name='heldout_F1')
    f1.to_frame().assign(n_train=counts.reindex(clf.classes_).values, n_genes=len(genes)).to_csv(out / f'panel_{key}_heldout_F1.tsv', sep='\t')
    print(key, len(genes), 'genes; labels', len(clf.classes_), 'identifiable', int((f1 >= F1_MIN).sum()), flush=True)
    for f in files:
        name = Path(f).name.replace('.counts.h5ad', '')
        s = ad.read_h5ad(f); s = s[s.obs.qc_pass.astype(bool).to_numpy()].copy()
        P = clf.predict_proba(sc.transform(lognorm(s[:, genes].X)))
        lab = clf.classes_[P.argmax(1)]; conf = P.max(1)
        ident = f1.reindex(lab).fillna(0).to_numpy() >= F1_MIN
        # unidentifiable fine labels fall back to their V7 major so blisa groups stay measurable
        maj = r.obs.drop_duplicates('label').set_index('label').major.astype(str)
        group = np.where(ident, lab, 'Unresolved_' + maj.reindex(lab).to_numpy())
        d = out / name; d.mkdir(exist_ok=True)
        pd.DataFrame({'cell': s.obs_names, 'x': s.obsm['spatial'][:, 0], 'y': s.obsm['spatial'][:, 1], 'label': lab,
                      'prob': conf, 'identifiable': ident, 'group': group, 'panel': key}).to_csv(d / 'cells.tsv.gz', sep='\t', index=False)
        sio.mmwrite(str(d / 'counts_genes_x_cells.mtx'), sp.csr_matrix(s.X).T.tocoo())
        pd.Series(s.var_names).to_csv(d / 'genes.tsv', index=False, header=False)
        print(name, s.n_obs, pd.Series(group).value_counts().head(8).to_dict(), flush=True)
