#!/usr/bin/env python3
"""Annotate single-cell spatial sections with the CT75 annotation scheme, and validate on held-out atlas cells.

The atlas node of a cell is (parent compartment, numeric Leiden cluster, identity):
- identity: official CellTypist Immune_All_Low, best match, no majority voting,
  mapped to canonical identity through the frozen reference-semantics ontology;
  kept only when its coarse reference parent belongs to the cell's parent
  compartment, otherwise the cell becomes the parent x cluster
  'reference_unresolved_or_conflict' node (same rule as the atlas). Non-immune
  parents use parent x cluster only.
- parent x cluster: cannot be recomputed off-atlas, so it is transferred with a
  CellTypist model trained (same package, SGD logistic regression) on the atlas
  reference cells restricted to the section's panel genes.
The identical function runs on held-out atlas cells restricted to the same panel
(raw counts, CP10K over panel genes, as for the spatial cells); per-node F1 on
those cells decides which CM members a panel can measure (pre-registered: F1>=0.5).

Run inside wehi_oldpath_sandbox.sh (reads frozen JSONs).
Usage: annotate_spatial_cells_v1.py --reference REF_DIR --model-dir MODELS --sections H5AD [H5AD ...] --output DIR
"""
import argparse, hashlib, json, os
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import celltypist

IMMUNE = ['T_NK', 'B_Plasma', 'Myeloid', 'Mast']


def load(p):
    return json.loads(Path(p).read_text())


def cp10k_log(X):
    X = sp.csr_matrix(X, dtype=np.float64)
    tot = np.asarray(X.sum(axis=1)).ravel(); s = np.zeros_like(tot); s[tot > 0] = 1e4 / tot[tot > 0]
    Y = sp.diags(s) @ X; Y = sp.csr_matrix(Y); Y.data = np.log1p(Y.data)
    return Y, tot


class Scheme:
    def __init__(self, v5):
        runs = Path(v5) / 'runs'
        fit = runs / '20261101_TME_immune_CM_main_CT_identity_own_v1'
        acc = load(load(fit / 'config/spec_v1.json')['input_acceptance']['path'])
        m = load(acc['annotation_acceptance']['path'])
        sem = load(m['reference_semantics_acceptance']['path'])
        ont = load(sem['exact_ontology']['path'])['models']['immune_low']
        pol = load(runs / '20261030_TME_shared_compartment_scope_v3/config/frozen_compartment_policy_v1.json')
        family = {x: k for k, v in pol['compartments'].items() for x in v['identity_parents']}
        self.label_identity = {r['label']: r['canonical_identity'] for r in ont['rows']}
        self.label_family = {r['label']: family.get(r['coarse_reference_parent']) for r in ont['rows']}
        parents = {r['code']: r['name'] for r in m['parents']}
        self.node = {(parents[n['parent_code']], int(n['numeric_cluster']), n['reference_label']): n['name'] for n in m['nodes']}
        model_sha = [x for x in sem['models'] if x['model'] == 'immune_low'][0]['reference_model_SHA256']
        self.model_sha = model_sha

    def assign(self, parent_cluster, ct_label):
        out = []
        for pc, lab in zip(parent_cluster, ct_label):
            p, c = pc.split('|C'); c = int(c)
            if p in IMMUNE:
                fid = self.label_identity[lab] if self.label_family.get(lab) == p else 'reference_unresolved_or_conflict'
            else:
                fid = 'numeric_cluster'
            out.append(self.node.get((p, c, fid), 'unseen_combination'))
        return np.array(out)


def panel_model(ref, genes, model_dir, seed):
    key = hashlib.sha256('\n'.join(genes).encode()).hexdigest()[:16]
    path = Path(model_dir) / f'parent_cluster_{key}.pkl'
    if not path.exists():
        r = ref[ref.obs.split == 'reference', genes]
        Y, tot = cp10k_log(r.X); keep = tot > 0
        a = ad.AnnData(Y[keep], obs=r.obs.loc[keep, ['parent_cluster']].copy(), var=pd.DataFrame(index=genes))
        np.random.seed(seed)
        mdl = celltypist.train(a, labels='parent_cluster', check_expression=True, use_SGD=True, mini_batch=True,
                               feature_selection=False, n_jobs=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
        Path(model_dir).mkdir(parents=True, exist_ok=True); mdl.write(str(path))
    return celltypist.models.Model.load(str(path)), key


def annotate(X, genes, pc_model, ct_model, scheme):
    Y, tot = cp10k_log(X); keep = tot > 0
    a = ad.AnnData(Y[keep], var=pd.DataFrame(index=genes))
    p1 = celltypist.annotate(a, model=pc_model, mode='best match', majority_voting=False).predicted_labels
    p2 = celltypist.annotate(a, model=ct_model, mode='best match', majority_voting=False).predicted_labels
    res = pd.DataFrame({'parent_cluster': p1.predicted_labels.astype(str).to_numpy(),
                        'celltypist_label': p2.predicted_labels.astype(str).to_numpy()})
    res['node'] = scheme.assign(res.parent_cluster, res.celltypist_label)
    res['parent'] = res.parent_cluster.str.split('|').str[0]
    return res, keep


def validate(ref, genes, pc_model, ct_model, scheme, ct75):
    h = ref[ref.obs.split == 'heldout', genes]
    res, keep = annotate(h.X, genes, pc_model, ct_model, scheme)
    truth = h.obs.loc[keep]
    res.index = truth.index
    rows = []
    for n in ct75:
        t = truth.node.to_numpy() == n; p = res.node.to_numpy() == n
        tp = int((t & p).sum()); prec = tp / p.sum() if p.sum() else np.nan; rec = tp / t.sum() if t.sum() else np.nan
        f1 = 2 * prec * rec / (prec + rec) if prec and rec and np.isfinite(prec) and np.isfinite(rec) else 0.0
        rows.append(dict(node=n, n_true=int(t.sum()), n_pred=int(p.sum()), precision=prec, recall=rec, F1=f1))
    pcacc = float((res.parent_cluster.to_numpy() == truth.parent_cluster.to_numpy()).mean())
    pacc = float((res.parent.to_numpy() == truth.parent.to_numpy()).mean())
    return pd.DataFrame(rows), dict(n_heldout=int(keep.sum()), parent_accuracy=pacc, parent_cluster_accuracy=pcacc,
                                    node_accuracy_all=float((res.node.to_numpy() == truth.node.to_numpy()).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--v5', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--ct-model', required=True)
    ap.add_argument('--model-dir', required=True)
    ap.add_argument('--sections', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--seed', type=int, default=20261102)
    a = ap.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    scheme = Scheme(a.v5)
    assert hashlib.sha256(Path(a.ct_model).read_bytes()).hexdigest() == scheme.model_sha, 'CellTypist model differs from atlas'
    ct_model = celltypist.models.Model.load(a.ct_model)
    ref = ad.read_h5ad(Path(a.reference) / 'atlas_reference_v1.h5ad')
    ct75 = pd.read_csv(Path(a.v5) / 'runs/20261101_TME_immune_CM_main_CT_identity_own_v1/outputs/seed_20261030/private_W_v1.tsv', sep='\t', index_col=0).index.tolist()
    validated = {}
    for f in a.sections:
        name = Path(f).name.replace('.counts.h5ad', '')
        if (out / f'{name}.cells.tsv.gz').exists():
            continue
        s = ad.read_h5ad(f)
        s = s[s.obs.qc_pass.astype(bool).to_numpy()]
        genes = [g for g in s.var_names if g in set(ref.var_names)]
        pc_model, key = panel_model(ref, genes, a.model_dir, a.seed)
        if key not in validated:
            vt, vs = validate(ref, genes, pc_model, ct_model, scheme, ct75)
            vt.to_csv(out / f'panel_{key}_heldout_node_F1.tsv', sep='\t', index=False)
            vs.update(panel=key, n_panel_genes=len(genes), first_section=name)
            json.dump(vs, open(out / f'panel_{key}_heldout_summary.json', 'w'), indent=2)
            validated[key] = vs
        res, keep = annotate(s[:, genes].X, genes, pc_model, ct_model, scheme)
        cells = s.obs.loc[keep, ['x_centroid', 'y_centroid']].copy()
        for c in ['patient', 'run', 'sample_id']:
            if c in s.obs:
                cells[c] = s.obs.loc[keep, c].astype(str).to_numpy()
        for c in res.columns:
            cells[c] = res[c].to_numpy()
        cells['panel'] = key
        cells.to_csv(out / f'{name}.cells.tsv.gz', sep='\t')
        print(name, len(cells), 'panel', key, cells.node.isin(ct75).mean().round(3), flush=True)


if __name__ == '__main__':
    main()
