#!/usr/bin/env python3
"""Upper bound for panel annotation: the same two-step CT75-scheme annotation with all 12,783 genes.

Held-out atlas cells, raw counts, CP10K over the gene set used. Diagnostic only:
it shows how separable the CT75 nodes are from expression at all, independent
of any spatial panel. Run inside wehi_oldpath_sandbox.sh.
Usage: diag_node_separability_v1.py --v5 OLD_V5 --reference REF_DIR --ct-model PKL --model-dir DIR --output DIR
"""
import argparse, json, sys
from pathlib import Path
import anndata as ad
import celltypist

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotate_spatial_cells_v1 import Scheme, panel_model, validate
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--v5', required=True); ap.add_argument('--reference', required=True)
    ap.add_argument('--ct-model', required=True); ap.add_argument('--model-dir', required=True)
    ap.add_argument('--output', required=True); ap.add_argument('--seed', type=int, default=20261102)
    a = ap.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    scheme = Scheme(a.v5); ct = celltypist.models.Model.load(a.ct_model)
    ref = ad.read_h5ad(Path(a.reference) / 'atlas_reference_v1.h5ad')
    ct75 = pd.read_csv(Path(a.v5) / 'runs/20261101_TME_immune_CM_main_CT_identity_own_v1/outputs/seed_20261030/private_W_v1.tsv', sep='\t', index_col=0).index.tolist()
    genes = list(ref.var_names)
    model, key = panel_model(ref, genes, a.model_dir, a.seed)
    t, s = validate(ref, genes, model, ct, scheme, ct75)
    t.to_csv(out / 'fullgene_heldout_node_F1.tsv', sep='\t', index=False)
    s.update(n_genes=len(genes), panel=key); json.dump(s, open(out / 'fullgene_heldout_summary.json', 'w'), indent=2)
    print(s); print((t.F1 >= .5).sum(), 'of', len(t), 'nodes F1>=0.5')


if __name__ == '__main__':
    main()
