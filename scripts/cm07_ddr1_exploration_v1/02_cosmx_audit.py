#!/usr/bin/env python3
"""CM07 CosMx 复核。探索性分析, 不改写原结果。A: 原 T3b 的密度假象; B: FOV 内公平比较; C: 全基因阴性对照。"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors as NN, BallTree
RUN = Path(sys.argv[1]); OUT = RUN / "outputs"; OUT.mkdir(exist_ok=True)
PX = 0.18  # 假定值: 由最近邻与 FOV 尺寸推断, 见报告
rng = np.random.default_rng(20261106)
T = RUN.parent / "communication_evidence_chain_20260906/tables"
m = pd.read_csv(T / "cosmx_cells_meta.tsv.gz", sep="\t").set_index("cell_id")
m["compartment"] = m.compartment.str.strip(); m["cell_type"] = m.cell_type.str.strip()
cells = (RUN / "inputs/cosmx_cells.txt").read_text().split(); genes = (RUN / "inputs/cosmx_genes.txt").read_text().split()
X = sio.mmread(RUN / "inputs/cosmx_counts.mtx").tocsc()  # genes x cells
m = m.loc[cells]; print("cells", len(m), "genes", len(genes), flush=True)
gi = {g: i for i, g in enumerate(genes)}
libsize = np.asarray(X.sum(0)).ravel()
print("nCount vs sum corr", np.corrcoef(libsize, m.nCount)[0, 1])
Xc = X.tocsr()
L = Xc.multiply(1 / np.maximum(libsize, 1)).tocsr(); L.data = np.log1p(L.data * 1e4)  # genes x cells, log CP10K
Lc = L.tocsc()
def g(n): return np.asarray(Lc[gi[n]].todense()).ravel() if False else np.asarray(L[gi[n]].todense()).ravel()
FIB = {"Fibroblast", "Myofibroblast"}
IMM = {"CD8.T.cell","NKT.cell","NK.cell","CD4.T.cell","B.cell","MNP.c.dendritic.cell","MNP.d.Tissue.macrophage","MNP.a.classical.monocyte.derived","MNP.b.non.classical.monocyte.derived","Peritubular.capillary.endothelium.1","Peritubular.capillary.endothelium.2","Ascending.vasa.recta.endothelium","Descending.vasa.recta.endothelium","Glomerular.endothelium","pDC","Mast.cell","Plasma.cell"}
m["is_tum_orig"] = (m.compartment == "Tumor") & ~m.cell_type.isin(FIB | IMM)
m["is_cd8"] = m.cell_type.isin({"CD8.T.cell", "NKT.cell"}); m["is_cd8_strict"] = m.cell_type == "CD8.T.cell"
# 簇身份: 各 letter 簇的 marker
mk = ["CA9","NDUFA4L2","VEGFA","VHL","EPCAM","KRT8","KRT18","LRP2","CUBN","SLC34A1","COL1A1","COL1A2","VIM","PTPRC","CD68","PECAM1","DDR1","HIF1A"]
mk = [x for x in mk if x in gi]
pb = pd.DataFrame({x: pd.Series(g(x), index=m.index).groupby(m.cell_type).mean() for x in mk}).round(2)
pb["n"] = m.groupby("cell_type").size(); pb["tum_frac"] = m.groupby("cell_type").is_tum_orig.mean().round(2); pb["nCount"] = m.groupby("cell_type").nCount.median()
pb.sort_values("n", ascending=False).to_csv(OUT / "cluster_markers.tsv", sep="\t"); print(pb.sort_values("n", ascending=False).head(14).to_string(), flush=True)
rows_a, rows_b, rows_c = [], [], []
K = 12; R_um = 50; R_px = R_um / PX
col1a1 = g("COL1A1"); ddr1 = g("DDR1")
gene_keep = [i for i, n in enumerate(genes) if not n.lower().startswith(("neg", "sys")) ]
for tissue, sub in m.groupby("tissue"):
    idx = m.index.get_indexer(sub.index); xy = sub[["x","y"]].values
    _, nb = NN(n_neighbors=K+1).fit(xy).kneighbors(xy); nb = nb[:,1:]
    colx = col1a1[idx][nb].mean(1)
    tum = sub.is_tum_orig.values
    inter = colx[tum] * ddr1[idx][tum]; thr = np.median(inter)
    hi = inter > thr
    cd8xy = xy[sub.is_cd8.values]; tumxy = xy[tum]
    def dist(sel):
        d, _ = NN(n_neighbors=1).fit(tumxy[sel]).kneighbors(cd8xy); return np.median(d)
    obs_hi, obs_lo = dist(hi), dist(~hi)
    # A: 标签置换, 保持 HI 数量, 全切片与 FOV 内两种
    fovs = sub.fov.values[tum]; perm = []
    for b in range(300):
        p = rng.permutation(hi); perm.append((dist(p), dist(~p)))
    perm = np.array(perm); permf = []
    for b in range(300):
        p = hi.copy()
        for f in np.unique(fovs):
            w = np.where(fovs == f)[0]; p[w] = rng.permutation(hi[w])
        permf.append((dist(p), dist(~p)))
    permf = np.array(permf)
    robs = np.log(obs_hi / obs_lo); r1 = np.log(perm[:,0] / perm[:,1]); r2 = np.log(permf[:,0] / permf[:,1])
    rows_a.append(dict(tissue=tissue, n_tum=int(tum.sum()), interface_median=float(thr), frac_hi=float(hi.mean()), frac_inter_zero=float((inter==0).mean()), n_cd8=len(cd8xy),
        obs_hi_px=obs_hi, obs_lo_px=obs_lo, obs_hi_um=obs_hi*PX, obs_lo_um=obs_lo*PX, null_hi_px=np.median(perm[:,0]), null_lo_px=np.median(perm[:,1]), ratio_obs=obs_hi/obs_lo, ratio_null=np.exp(np.median(r1)), ratio_null_fov=np.exp(np.median(r2)),
        p_vs_null=float((1+(r1>=robs).sum())/(1+len(r1))), p_vs_null_fov=float((1+(r2>=robs).sum())/(1+len(r2)))))
    # B: 肿瘤细胞周围 50 um 内 CD8 数; FOV×簇 内比较, FOV 聚类 SE
    tree = BallTree(cd8xy); ncd8 = tree.query_radius(tumxy, R_px, count_only=True).astype(float)
    # 排除距切片/FOV 太近的边缘不易; 以 FOV FE 吸收
    d0 = pd.DataFrame(dict(y=ncd8, ddr1=ddr1[idx][tum], col=colx[tum], fov=fovs, cl=sub.cell_type.values[tum], lib=np.log1p(libsize[idx][tum])))
    def fe_fit(d, xs, block):
        D = d[["y"]+xs].astype(float).copy(); 
        Dm = D - D.groupby(block.values).transform("mean")
        Xm = Dm[xs].values; ym = Dm.y.values
        XtX = Xm.T @ Xm; bt = np.linalg.solve(XtX, Xm.T @ ym); res = ym - Xm @ bt
        cl = d.fov.values; meat = np.zeros_like(XtX); U = Xm * res[:,None]
        for f in np.unique(cl): s = U[cl==f].sum(0); meat += np.outer(s, s)
        I = np.linalg.inv(XtX); V = I @ meat @ I; G = len(np.unique(cl)); V *= G/(G-1)
        return bt, np.sqrt(np.diag(V))
    sdz = lambda v: (v - v.mean()) / v.std()
    for name, block in [("fov", d0.fov.astype(str)), ("fov_x_cluster", d0.fov.astype(str)+"|"+d0.cl)]:
        dd = d0.assign(ddr1=sdz(d0.ddr1), col=sdz(d0.col), lib=sdz(d0.lib)); dd["cxd"] = dd.ddr1*dd.col
        for xs in (["col","ddr1","cxd"], ["col","ddr1","cxd","lib"]):
            bt, se = fe_fit(dd, xs, block)
            for x, b_, s_ in zip(xs, bt, se): rows_b.append(dict(tissue=tissue, block=name, model="+".join(xs), term=x, beta=b_, se=s_, t=b_/s_))
    # 主簇内分层
    for cl, dc in d0.groupby("cl"):
        if len(dc) < 500: continue
        dd = dc.assign(ddr1=sdz(dc.ddr1), col=sdz(dc.col), lib=sdz(dc.lib)); dd["cxd"]=dd.ddr1*dd.col
        if dd.ddr1.std()==0 or not np.isfinite(dd.ddr1.std()): continue
        bt, se = fe_fit(dd, ["col","ddr1","cxd","lib"], dd.fov.astype(str))
        for x,b_,s_ in zip(["col","ddr1","cxd","lib"],bt,se): rows_b.append(dict(tissue=tissue, block="fov_within_"+cl, model="col+ddr1+cxd+lib", term=x, beta=b_, se=s_, t=b_/s_))
    # C: 全基因阴性对照。基线 = FOV×簇 + 文库大小
    dd = d0.copy(); blk = dd.fov.astype(str)+"|"+dd.cl; codes = pd.factorize(blk)[0]; nb_ = codes.max()+1
    def dm(v):
        s = np.bincount(codes, v, nb_)/np.bincount(codes, minlength=nb_); return v - s[codes]
    y = dm(dd.y.values); lib = dm(sdz(dd.lib).values)
    ry = y - lib*(lib@y)/(lib@lib)
    Lt = L[:, idx].tocsr()
    # 肿瘤细胞自身表达的基因 vs CD8 邻居
    own = {}; nbm = {}
    tumi = np.where(tum)[0]
    from scipy.sparse import csr_matrix
    n_nb = len(tumi)*K; rowsn = np.repeat(np.arange(len(tumi)), K); colsn = nb[tumi].ravel()
    A = csr_matrix((np.ones(n_nb)/K, (rowsn, colsn)), shape=(len(tumi), len(idx)))
    NBX = (A @ Lt.T).tocsc()  # tum x genes: 邻居平均
    OWN = Lt[:, tumi].T.tocsc()  # tum x genes
    dt = np.ones(len(genes))
    for gidx in gene_keep:
        o = np.asarray(OWN[:,gidx].todense()).ravel(); n_ = np.asarray(NBX[:,gidx].todense()).ravel()
        det = (o>0).mean(); 
        rec = dict(tissue=tissue, gene=genes[gidx], det_own=det, det_nb=(n_>0).mean(), mean_own=o.mean(), mean_nb=n_.mean())
        for lab, v in (("own",o),("nb",n_)):
            if v.std()==0: rec["b_"+lab]=np.nan; continue
            z = dm(sdz(v)); z = z - lib*(lib@z)/(lib@lib); rec["b_"+lab] = float(z@ry/(z@z))
        rows_c.append(rec)
    print(tissue, "done", flush=True)
pd.DataFrame(rows_a).to_csv(OUT/"A_T3b_density_null.tsv",sep="\t",index=False)
pd.DataFrame(rows_b).to_csv(OUT/"B_fov_within_models.tsv",sep="\t",index=False)
pd.DataFrame(rows_c).to_csv(OUT/"C_allgene_control.tsv",sep="\t",index=False)
