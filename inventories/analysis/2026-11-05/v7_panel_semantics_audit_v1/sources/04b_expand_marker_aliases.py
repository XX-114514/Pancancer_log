#!/usr/bin/env python3
"""Extend stage 04's single hardcoded KIM1->HAVCR1 patch into a small,
per-GSE-fallback alias table covering four more legacy symbols found by an
independent mygene.info (HGNC-backed) audit of all 359 marker-panel genes
against the 54 per-GSE source files:

    SELENOP -> SEPP1   (recoverable in 18/20 GSEs currently missing SELENOP)
    ACKR1   -> DARC    (recoverable in 5/8 GSEs currently missing ACKR1)
    JCHAIN  -> IGJ      (recoverable in 5/6 GSEs currently missing JCHAIN)
    CXCL8   -> IL8      (recoverable in 5/6 GSEs currently missing CXCL8)

Unlike KIM1 (never a valid symbol in any of the 54 sources, already
resolved by stage 04 renaming that column to HAVCR1 outright), these four
genes are already present under their primary symbol in most GSEs. This
stage therefore resolves per (gene, GSE): use the primary symbol's counts
where the source file has it, fall back to the alias's counts only where
the primary symbol is absent but the alias is present, and leaves the
column unavailable (masked, structural zero) otherwise. Column names are
NOT renamed (unlike KIM1->HAVCR1) since the primary symbol is already the
correct canonical name used by pancancer_evidence_marker_specs.py.

Reads stage 04's output (v1); does not modify or delete it (append-only
contract). Writes a v2 object with only these 5 columns' counts changed;
all other 354 marker columns are copied through unchanged from v1.
"""

from __future__ import annotations

import os
import socket
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import zarr

from common import ROOT, append_status, ensure_layout, load_config, now_iso, refuse_existing, resolve_project_path, write_json_new


ALIASES: dict[str, str] = {
    "SELENOP": "SEPP1",
    "ACKR1": "DARC",
    "JCHAIN": "IGJ",
    "CXCL8": "IL8",
}
# KIM1 is not handled here: stage 04 already resolved it by renaming the
# column to HAVCR1 outright (KIM1 was never a valid symbol in any of the 54
# sources, so there is no "primary vs. fallback" choice to make), and
# 06_finalize_evidence_annotation.py's own ALIASES dict already redirects any
# downstream "KIM1" marker-panel lookup to that renamed HAVCR1 column.


def source_symbols(adata: ad.AnnData) -> pd.Index:
    names = pd.Index(adata.var_names.astype(str))
    looks_ensembl = names.str.match(r"^ENSG\d+(\.\d+)?$").mean() > 0.5
    symbol_col = next((c for c in ("gene_symbol", "gene_symbols", "gene_name", "symbol", "SYMBOL") if c in adata.var.columns), None)
    if looks_ensembl and symbol_col is not None:
        names = pd.Index(adata.var[symbol_col].astype(str))
    return names.str.strip().str.upper()


def gene_column_values(source: ad.AnnData, symbols: pd.Index, gene: str, wanted_obs_names: pd.Series) -> tuple[np.ndarray, bool]:
    cols = np.flatnonzero(symbols == gene)
    if len(cols) == 0:
        return np.zeros(len(wanted_obs_names), dtype=np.float32), False
    view = source[:, cols]
    raw = view.layers["counts"] if "counts" in source.layers else view.X
    raw = raw.tocsr() if sp.issparse(raw) else sp.csr_matrix(raw)
    if raw.shape[1] > 1:
        raw = sp.csr_matrix(np.asarray(raw.sum(axis=1)))
    source_index = pd.Index(source.obs_names.astype(str))
    positions = source_index.get_indexer(wanted_obs_names.to_numpy())
    if (positions < 0).any():
        raise ValueError(f"missing {int((positions < 0).sum())} source rows for gene {gene}")
    values = np.asarray(raw[positions, :].toarray()).ravel().astype(np.float32)
    return values, True


def lognormalize(counts: sp.csr_matrix, totals: np.ndarray) -> sp.csr_matrix:
    matrix = counts.astype(np.float32, copy=True).tocsr()
    scales = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    matrix.data *= np.repeat(scales, np.diff(matrix.indptr))
    np.log1p(matrix.data, out=matrix.data)
    matrix.eliminate_zeros()
    return matrix


def main() -> None:
    started = time.time()
    started_iso = now_iso()
    ensure_layout()
    config = load_config()
    v1_path = ROOT / "objects" / "wide_marker_panel_full_library_masked_v1.h5ad"
    output = ROOT / "objects" / "wide_marker_panel_full_library_masked_v2.h5ad"
    alias_path = ROOT / "tables" / "marker_alias_resolution_v2.tsv"
    summary_path = ROOT / "reports" / "04b_expand_marker_aliases_summary.json"
    refuse_existing(output, alias_path, summary_path)
    if not v1_path.exists():
        raise FileNotFoundError(f"Stage 04 output missing: {v1_path}")
    append_status("04b", "expand_marker_aliases", "running", started_iso,
                  input_ref=str(v1_path), job_id=str(os.getpid()))

    wide = ad.read_h5ad(v1_path)
    marker_names = pd.Index(wide.var_names.astype(str).str.upper())
    for gene in ALIASES:
        if gene not in marker_names:
            raise ValueError(f"Expected marker column {gene} in stage 04 output")
    counts = wide.layers["counts"].tocsr() if sp.issparse(wide.layers["counts"]) else sp.csr_matrix(wide.layers["counts"])
    availability = wide.varm["availability_by_gse"].copy()
    gse_order = list(wide.uns["availability_gse_order"])

    obs = wide.obs
    gse_values = obs["gse_id"].astype(str).to_numpy()
    merge_root_zarr = zarr.open_group(str(resolve_project_path(config["merge_zarr"])), mode="r")
    merge_obs = ad.io.read_elem(merge_root_zarr["obs"])
    if not merge_obs.index.equals(wide.obs_names):
        merge_obs = merge_obs.reindex(wide.obs_names)

    n_before = {gene: int(availability.loc[gene].sum()) for gene in ALIASES}
    col_values = {gene: np.zeros(len(wide), dtype=np.float32) for gene in ALIASES}
    col_available = {gene: availability.loc[gene].copy() for gene in ALIASES}

    for gse_id in gse_order:
        global_pos = np.flatnonzero(gse_values == gse_id)
        if len(global_pos) == 0:
            continue
        group = merge_obs.iloc[global_pos]
        source_file = Path(str(group["phase07_source_h5ad"].astype(str).unique()[0]))
        source = sc.read_h5ad(source_file, backed="r")
        symbols = source_symbols(source)
        symbol_set = set(symbols)
        wanted = group["phase07_source_obs_name"].astype(str)
        for gene, alias in ALIASES.items():
            if gene in symbol_set:
                values, found = gene_column_values(source, symbols, gene, wanted)
            elif alias in symbol_set:
                values, found = gene_column_values(source, symbols, alias, wanted)
            else:
                values, found = np.zeros(len(wanted), dtype=np.float32), False
            col_values[gene][global_pos] = values
            col_available[gene][gse_id] = found
        source.file.close()
        print(f"scanned {gse_id} for {len(ALIASES)} alias-fallback genes", flush=True)

    alias_rows: list[dict[str, object]] = []
    for gene, alias in ALIASES.items():
        gene_col = int(marker_names.get_loc(gene))
        counts_col = sp.csr_matrix(col_values[gene].reshape(-1, 1))
        counts = sp.hstack([counts[:, :gene_col], counts_col, counts[:, gene_col + 1:]], format="csr")
        availability.loc[gene] = col_available[gene]
        n_after = int(col_available[gene].sum())
        alias_rows.append({
            "requested_symbol": gene, "resolved_alias": alias, "authority": "HGNC_current_symbol_or_validated_alias",
            "n_gse_available_before": n_before[gene], "n_gse_available_after": n_after,
            "n_cells_expressing_after": int((col_values[gene] != 0).sum()),
        })
        print(f"{gene} (-> {alias} fallback): GSE availability {n_before[gene]} -> {n_after}", flush=True)

    totals = pd.to_numeric(wide.obs["total_counts_full_library"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)
    corrected_x = lognormalize(counts, totals)

    var = wide.var.copy()
    for gene, alias in ALIASES.items():
        var.loc[gene, "marker_alias_fallback"] = alias
    if "marker_alias_fallback" not in var.columns:
        var["marker_alias_fallback"] = pd.NA
    var["n_gse_available"] = availability.sum(axis=1).astype(np.int16).to_numpy()
    var["fraction_gse_available"] = availability.mean(axis=1).astype(np.float32).to_numpy()

    corrected = ad.AnnData(X=corrected_x, obs=wide.obs.copy(), var=var)
    corrected.layers["counts"] = counts.astype(np.float32)
    corrected.varm["availability_by_gse"] = availability
    corrected.uns["availability_gse_order"] = gse_order
    corrected.uns["normalization"] = dict(wide.uns["normalization"])
    corrected.uns["marker_aliases"] = dict(ALIASES)
    corrected.uns["marker_alias_provenance"] = (
        "per-(gene,GSE) fallback: primary symbol used where the source file has it; "
        "alias used only where primary symbol is absent but alias is present; "
        "aliases validated by reverse mygene.info lookup (unambiguous only)"
    )
    ad.settings.allow_write_nullable_strings = True
    corrected.write_h5ad(output, compression="gzip", compression_opts=4)

    pd.DataFrame(alias_rows).to_csv(alias_path, sep="\t", index=False)

    ended_iso = now_iso()
    duration = time.time() - started
    payload = {
        "run_id": config["run_id"], "stage": "04b_expand_marker_aliases", "status": "success",
        "start_time": started_iso, "end_time": ended_iso, "duration_sec": duration,
        "params": {"aliases": ALIASES},
        "inputs": [str(v1_path)], "outputs": [str(output), str(alias_path)],
        "warnings": ["v1 object retained unchanged; only these 5 marker columns differ in v2"],
        "error": None, "source_scripts": [str(Path(__file__).resolve())], "logic_version": config["logic_version"],
        "hostname": socket.gethostname(), "pid": os.getpid(), "shape": [int(corrected.n_obs), int(corrected.n_vars)],
        "alias_recovery": alias_rows,
    }
    write_json_new(summary_path, payload)
    append_status("04b", "expand_marker_aliases", "success", started_iso, ended_iso, round(duration, 3),
                  str(v1_path), str(output), str(os.getpid()), "")
    print(payload, flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        import logging
        logging.exception("expand marker aliases failed")
        append_status("04b", "expand_marker_aliases", "failed", error_log=repr(exc), job_id=str(os.getpid()))
        raise
