#!/usr/bin/env python3
"""Inventory of downloaded non-Xenium sources: per-section cell/spot and gene counts, coordinates,
author annotations, and coverage of the same gene lists used for Xenium panel qualification.

Inspection only: files are read, nothing is converted, filtered or normalised.
Usage: python -I inspect_sources.py RAW_ROOT FROZEN_PLAN_JSON OUT_DIR
Writes OUT_DIR/source_sections.tsv, OUT_DIR/source_gene_coverage.tsv and OUT_DIR/SOURCES.md.
"""
import gzip
import io
import json
import re
import sys
import tarfile
from pathlib import Path

import h5py
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from panel_coverage import CD3, gene_lists  # noqa: E402

# Probe/control features that are not genes on each platform.
CONTROL = re.compile(r"^(Neg|NegPrb|SystemControl|Blank|BLANK|NegControl|Unassigned|fov$|cell_ID$|cell$|Unnamed)")


def csv_header(path):
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split(",")


def count_lines(path):
    with gzip.open(path, "rb") as fh:
        return sum(1 for _ in fh) - 1


def tenx_h5_shape(path):
    """(n_barcodes, gene names) from a 10x feature-barcode h5, gzipped or not."""
    raw = Path(path)
    data = gzip.open(raw).read() if raw.suffix == ".gz" else None
    with h5py.File(io.BytesIO(data) if data else raw, "r") as f:
        m = f["matrix"]
        names = [x.decode() for x in m["features"]["name"][:]]
        types = [x.decode() for x in m["features"]["feature_type"][:]]
        genes = [n for n, t in zip(names, types) if t == "Gene Expression"]
        return int(m["shape"][1]), genes


def cosmx(d):
    rows = []
    for expr in sorted(d.rglob("*exprMat_file.csv.gz")):
        header = csv_header(expr)
        meta = next(iter(sorted(expr.parent.glob(expr.name.replace("exprMat_file", "metadata_file")))), None)
        meta_cols = csv_header(meta) if meta else []
        genes = [g for g in header if not CONTROL.match(g)]
        rows.append(dict(section=expr.name.split("_exprMat")[0], unit="cell", n_obs=count_lines(expr), genes=genes,
                         n_controls=sum(1 for g in header if CONTROL.match(g)),
                         coordinates="CenterX_global_px" in meta_cols,
                         author_annotation=any(expr.parent.glob("*.RDS*")),
                         notes="leading unnamed index column" if header and header[0] in ("", "Unnamed: 0") else ""))
    return rows


def merfish(d):
    rows = []
    for cbg in sorted(d.rglob("*cell_by_gene.csv.gz")):
        header = csv_header(cbg)
        meta = next(iter(sorted(cbg.parent.glob("*cell_metadata.csv.gz"))), None)
        meta_cols = csv_header(meta) if meta else []
        genes = [g for g in header[1:] if not CONTROL.match(g)]
        rows.append(dict(section=cbg.name.split("_cell_by_gene")[0], unit="cell", n_obs=count_lines(cbg), genes=genes,
                         n_controls=sum(1 for g in header if CONTROL.match(g)),
                         coordinates=any(c in meta_cols for c in ("center_x", "x")),
                         author_annotation=False, notes="cellpose segmentation" if "cellpose" in cbg.name else ""))
    return rows


def h5ad_sections(d):
    import anndata as ad
    rows = []
    for p in sorted(d.rglob("*.h5ad")):
        a = ad.read_h5ad(p, backed="r")
        x = a.X if a.X is not None else None
        sample = x[: min(200, a.n_obs)]
        sample = sample.toarray() if hasattr(sample, "toarray") else sample
        integer = bool((sample == sample.round()).all()) if sample is not None else None
        label_cols = [c for c in a.obs.columns if re.search(r"(type|annot|label|cluster|niche)", c, re.I)]
        rows.append(dict(section=p.stem, unit="cell_or_spot", n_obs=a.n_obs, genes=list(map(str, a.var_names)), n_controls=0,
                         coordinates="spatial" in a.obsm, author_annotation=bool(label_cols),
                         notes=f"X integer-valued in first rows={integer}; layers={list(a.layers.keys())}; "
                               f"label cols={label_cols[:6]}; obsm={list(a.obsm.keys())}"))
        a.file.close()
    return rows


def visium(d):
    rows = []
    for h5 in sorted(d.rglob("*filtered_feature_bc_matrix.h5")):
        n, genes = tenx_h5_shape(h5)
        g = h5.parent
        rows.append(dict(section=h5.name.split("_filtered")[0], unit="spot", n_obs=n, genes=genes, n_controls=0,
                         coordinates=any(g.glob(h5.name.split("_filtered")[0] + "_tissue_positions*")),
                         author_annotation=any(g.glob(h5.name.split("_filtered")[0] + "_*annotation*")), notes=""))
    return rows


def tenx_h5_gz(d):
    rows = []
    for h5 in sorted(d.rglob("*filtered_feature_bc_matrix.h5.gz")):
        n, genes = tenx_h5_shape(h5)
        rows.append(dict(section=h5.name.split("_filtered")[0], unit="nucleus_or_cell", n_obs=n, genes=genes,
                         n_controls=0, coordinates=False, author_annotation=False, notes=""))
    return rows


def tar_samples(d):
    """GSE176021-style per-library tar.gz: list members only (no extraction)."""
    rows = []
    for t in sorted(d.rglob("*.tar.gz")):
        with tarfile.open(t) as tf:
            names = tf.getnames()
        kind = "vdj" if ".vdj." in t.name else "gex"
        rows.append(dict(section=t.name.replace(".tar.gz", ""), unit=f"library_{kind}", n_obs=None, genes=[],
                         n_controls=0, coordinates=False, author_annotation=False,
                         notes=";".join(sorted({Path(n).name for n in names})[:5])))
    return rows


def mtx_series(d):
    rows = []
    for m in sorted(d.rglob("*counts.mtx.gz")):
        with gzip.open(m, "rt") as fh:
            line = fh.readline()
            while line.startswith("%"):
                line = fh.readline()
        n_feat, n_cells, nnz = map(int, line.split())
        feat = m.with_name(m.name.replace("counts.mtx.gz", "feature.tsv.gz"))
        genes = [l.rstrip("\n").split("\t")[-1] for l in gzip.open(feat, "rt")] if feat.exists() else []
        rows.append(dict(section=m.name.replace(".counts.mtx.gz", ""), unit="cell", n_obs=n_cells, genes=genes,
                         n_controls=0, coordinates=False, author_annotation=False,
                         notes=f"features={n_feat}; nnz={nnz}"))
    return rows


def imc(d):
    rows = []
    for z in sorted(d.glob("*.zip")):
        import zipfile
        try:
            with zipfile.ZipFile(z) as zf:
                names = zf.namelist()
        except zipfile.BadZipFile:
            rows.append(dict(section=z.name, unit="archive", n_obs=None, genes=[], n_controls=0, coordinates=False,
                             author_annotation=False, notes="incomplete or unreadable zip"))
            continue
        tops = sorted({n.split("/")[1] if n.count("/") > 1 else n.split("/")[0] for n in names})
        rows.append(dict(section=z.name, unit="archive", n_obs=None, genes=[], n_controls=0, coordinates=False,
                         author_annotation=False, notes=f"{len(names)} members; top entries: {tops[:12]}"))
    return rows


READERS = {"CosMx": cosmx, "MERFISH": merfish, "LUAD_spatial": h5ad_sections, "Visium": visium,
           "snRNA": tenx_h5_gz, "scRNA_NSCLC": tar_samples, "scRNA_TNBC": mtx_series, "IMC": imc}


def reader_for(name):
    for key, fn in READERS.items():
        if key in name:
            return fn
    return None


def main(raw_root, frozen_plan, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lists = gene_lists(frozen_plan)
    sec_rows, cov_rows = [], []
    for src in sorted(p for p in Path(raw_root).iterdir() if p.is_dir()):
        fn = reader_for(src.name)
        if fn is None:
            continue
        try:
            rows = fn(src)
        except Exception as e:  # keep inventory going; the failure is reported, not hidden
            rows = [dict(section="(source)", unit="", n_obs=None, genes=[], n_controls=0, coordinates=False,
                         author_annotation=False, notes=f"READ FAILED: {type(e).__name__}: {e}")]
        for r in rows:
            genes = set(r.pop("genes"))
            sec_rows.append(dict(source=src.name, n_genes=len(genes) or None, **r))
            if genes:
                for name, (glist, _) in lists.items():
                    cov_rows.append(dict(source=src.name, section=r["section"], gene_list=name,
                                         present=sum(g in genes for g in glist), total=len(glist),
                                         missing=",".join(g for g in glist if g not in genes)))
                cov_rows.append(dict(source=src.name, section=r["section"], gene_list="CD3", total=len(CD3),
                                     present=sum(g in genes for g in CD3),
                                     missing=",".join(g for g in CD3 if g not in genes)))
    sec = pd.DataFrame(sec_rows)
    cov = pd.DataFrame(cov_rows)
    sec.to_csv(out_dir / "source_sections.tsv", sep="\t", index=False)
    cov.to_csv(out_dir / "source_gene_coverage.tsv", sep="\t", index=False)

    agg = sec.groupby("source").agg(sections=("section", "size"), units=("unit", "first"),
                                    n_obs=("n_obs", lambda s: int(s.sum()) if s.notna().any() else None),
                                    genes_min=("n_genes", "min"), genes_max=("n_genes", "max"),
                                    coordinates=("coordinates", "all"), author_annotation=("author_annotation", "any"))
    lines = ["# Downloaded source inventory (inspection only)", "", agg.to_markdown(), ""]
    if not cov.empty:
        piv = (cov.assign(cell=cov["present"].astype(str) + "/" + cov["total"].astype(str))
               .drop_duplicates(["source", "gene_list"]).pivot(index="source", columns="gene_list", values="cell"))
        lines += ["## Gene-list coverage (first section per source)", "", piv.to_markdown(), ""]
    notes = sec[sec["notes"].astype(str).str.len() > 0][["source", "section", "notes"]]
    lines += ["## Notes", "", notes.to_markdown(index=False), ""]
    (out_dir / "SOURCES.md").write_text("\n".join(lines))
    print(f"{len(sec)} sections from {sec['source'].nunique()} sources -> {out_dir}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
