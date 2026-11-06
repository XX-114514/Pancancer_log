import sys; from pathlib import Path; import pandas as pd
IN = Path(sys.argv[1]) / "inputs"
cv = pd.read_csv(IN / "covariates_part1.tsv", sep="\t").set_index("sample")
es = pd.read_csv(IN / "estimate.tsv", sep="\t"); print(es.head(2)); 
col = "sample" if "sample" in es.columns else es.columns[0]
es = es.set_index(col); cv["estimate_stromal"] = es["stromal"].reindex(cv.index); cv["estimate_immune"] = es["immune"].reindex(cv.index)
cv["caf_fraction"] = pd.read_csv(IN / "epic_caf.tsv", sep="\t").set_index("sample").caf_fraction.reindex(cv.index)
cv.reset_index().to_csv(IN / "covariates.tsv", sep="\t", index=False); print(cv.notna().mean().round(3).to_dict())
