#!/usr/bin/env python3
"""E2 inputs. No survival model is fitted here. Reuses the frozen loaders so rows match the published run."""
import sys, json, re, hashlib
from pathlib import Path
import numpy as np, pandas as pd
RUN = Path(sys.argv[1]); P = RUN.parent; BC = P / "bulk_cohorts_20260904"; IN = RUN / "inputs"; RES = RUN / "resources"
sys.path.insert(0, str(BC)); sys.dont_write_bytecode = True
import tcga_cox as tc
expr = tc.load_expression(); surv = tc.load_survival()
common = [p for p in surv.index if p in expr.columns]; surv = surv.loc[common]
vc = surv.cancer.value_counts(); keep = surv.cancer.isin(vc[vc >= 30].index)
pts = []
for c, g in surv[keep].groupby("cancer"): pts += list(g.index)
E = expr[pts]; surv = surv.loc[pts]
print("matrix", E.shape, "cancers", surv.cancer.nunique(), flush=True)
np.asarray(E.values, dtype=np.float32).tofile(IN / "expr.bin")
(IN / "genes.txt").write_text("\n".join(E.index)); (IN / "samples.txt").write_text("\n".join(E.columns))
clin = surv[["cancer", "age", "OS", "OS.time", "PFI", "PFI.time"]].copy(); clin.insert(0, "patient", clin.index); clin.insert(0, "sample", clin.index)
clin.to_csv(IN / "clin.tsv", sep="\t", index=False)
# 父程序: 读取已冻结的逐患者打分; 原脚本建模前做全局 z 标准化, 这里同样处理
sc = pd.read_csv(P / "communication_evidence_chain_20260906/tables/tcga_proxy_scores_per_patient.tsv.gz", sep="\t", index_col=0)
assert set(pts) == set(sc.index), "patient set differs from frozen score table"
sc = sc.loc[pts]; par = sc[[c for c in sc.columns if c.startswith("sig_")]].copy(); par.columns = [c[4:] for c in par.columns]
par = (par - par.mean()) / par.std(); par.insert(0, "sample", par.index); par.to_csv(IN / "parents.tsv", sep="\t", index=False)
sl = pd.read_csv(BC / "signatures/signatures_long.tsv", sep="\t")[["signature", "gene"]]; sl.to_csv(IN / "parent_genes.tsv", sep="\t", index=False)
# 轴表: 基因定义逐字取自 03_bulk_proxy_models.py 的 AXES
g = lambda s: ",".join(s.split("_"))
rows = [("SPP1_macrophage","primary","CM01","SPP1","CD44,ITGAV,ITGB1","adverse"),
        ("LGALS1_myCAF_CD8","primary","CM08","LGALS1","PTPRC,CD69","adverse"),
        ("collagen_DDR1","control","desmoplasia_myCAF","COL1A1,COL1A2,COL3A1","DDR1","none"),
        ("COL1A1_only_DDR1","control","desmoplasia_myCAF","COL1A1","DDR1","none")]
fa = pd.read_csv(P / "communication_evidence_chain_20260906/tables/focus_axes_definition.tsv", sep="\t")
seen = set()
for r in fa.itertuples():
    key = (r.cm, r.ligand, r.receptor)
    if key in seen: continue
    seen.add(key); rows.append((f"{r.cm}|{r.ligand}|{r.receptor}", "focus", r.cm, g(r.ligand), g(r.receptor), "none"))
pr = pd.read_csv(P / "frozen/cm_core/communication_research_priorities.tsv", sep="\t")
n_pr = len(pr)
for r in pr.itertuples():
    key = (r.cm, r.ligand_complex, r.receptor_complex)
    if key in seen: continue
    seen.add(key); rows.append((f"{r.cm}|{r.ligand_complex}|{r.receptor_complex}", "priority", r.cm, g(r.ligand_complex), g(r.receptor_complex), "none"))
ax = pd.DataFrame(rows, columns=["axis_id","family","parent","ligand_genes","receptor_genes","direction"])
ax.to_csv(IN / "axes.tsv", sep="\t", index=False)
allg = set(E.index); miss = {a.axis_id: [x for x in (a.ligand_genes + "," + a.receptor_genes).split(",") if x not in allg] for a in ax.itertuples()}
miss = {k: v for k, v in miss.items() if v}
print("axes", ax.family.value_counts().to_dict(), "priority rows", n_pr, "SPP1 in priority", pr.ligand_complex.str.contains("SPP1").sum(), "missing genes", miss, flush=True)
import liana
res = liana.resource.select_resource("consensus").drop_duplicates()
lr = pd.DataFrame({"ligand_genes": res.ligand.map(g), "receptor_genes": res.receptor.map(g)}).drop_duplicates()
lr.to_csv(IN / "lr_resource.tsv", sep="\t", index=False); print("lr pairs", len(lr), "liana", liana.__version__, flush=True)
# 协变量: 只取 01 原发肿瘤
s = pd.read_csv(BC / "tcga/TCGA_survival_S1.tsv", sep="\t"); s["patient"] = s["_PATIENT"].str[:12]; s = s.drop_duplicates("patient").set_index("patient")
def stage(v):
    m = re.fullmatch(r"(?:Stage )?(IV|III|II|I)[A-C]?\d?", str(v).strip()); return m.group(1) if m else None
st = s["ajcc_pathologic_tumor_stage"].map(stage); st = st.where(st.notna(), s["clinical_stage"].map(stage))
cv = pd.DataFrame(index=pts); cv["sex"] = s["gender"].reindex(pts); cv["stage"] = st.reindex(pts)
lk = pd.read_csv(RES / "TCGA_all_leuk_estimate.masked.20170107.tsv", sep="\t", header=None, names=["c","id","v"])
lk = lk[lk.id.str[13:15] == "01"]; cv["leukocyte_fraction"] = lk.groupby(lk.id.str[:12]).v.mean().reindex(pts)
cb = pd.read_csv(RES / "TCGA.Kallisto.fullIDs.cibersort.relative.tsv", sep="\t"); cb["id"] = cb.SampleID.str.replace(".", "-", regex=False)
cb = cb[cb.id.str[13:15] == "01"]; cb["mac"] = cb[[c for c in cb.columns if c.startswith("Macrophages")]].sum(axis=1)
by = cb.groupby(cb.id.str[:12])[["mac", "T.cells.CD8"]].mean()
cv["macrophage_fraction"] = by.mac.reindex(pts) * cv.leukocyte_fraction; cv["cd8_fraction"] = by["T.cells.CD8"].reindex(pts) * cv.leukocyte_fraction
ab = pd.read_csv(RES / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t"); ab = ab[(ab["call status"] == "called") & (ab.array.str[13:15] == "01")]
cv["purity_absolute"] = ab.groupby(ab.array.str[:12]).purity.mean().reindex(pts)
cv.insert(0, "sample", cv.index); cv.to_csv(IN / "covariates_part1.tsv", sep="\t", index=False)
print("covariate coverage", cv.notna().mean().round(3).to_dict(), "stage", cv.stage.value_counts(dropna=False).to_dict(), flush=True)
