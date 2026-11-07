# Annotation v2 (2026-11-07, user decision): label transfer for low-depth cells within one lineage object. These cells stay in
# the lineage total; their subtype is marked as inferred and is left out of the main subtype proportions (sensitivity analysis only).
#  query     : cells of the pool flagged low depth in round 1 (r1_action) or sitting in a keep_flag_low_depth cluster
#  reference : other cells in identity clusters (action keep or keep_tissue); state, source-restricted and sample-specific clusters are not used
#  method    : 30 nearest reference cells in X_emb (query cells are never neighbours), inverse-distance vote; confidence = top label's
#              weight share; below MIN_CONF the cell is <pool>_unresolved_low_depth
#  validation: per label, the 10% of reference cells with fewest genes are held out and predicted from the rest. They are the closest
#              available stand-in for low depth but are still deeper than the query, so the accuracy is an upper bound.
# Usage: transfer_low_depth_v1.py OBJ_H5AD REVIEW_TSV POOL CELL_LABELS OUT_PREFIX SEED   (env COLLAPSE="old=new,old=new", EXTRA_REF_ACTIONS)
import sys, os, re
import numpy as np, pandas as pd, h5py, pynndescent
from anndata.io import read_elem

obj, rev, pool, cl, out, seed = sys.argv[1:7]; seed = int(seed)
K, MIN_CONF = 30, 0.6
R = pd.read_csv(rev, sep='\t', dtype={'cluster': str}); R = R[R.pool == pool].set_index('cluster')
with h5py.File(obj, 'r') as f:
    o = read_elem(f['obs'])[['leiden_r1', 'n_genes']] if 'n_genes' in f['obs'] else read_elem(f['obs'])
    Z = f['obsm/X_emb'][:]
genes_col = next(c for c in ('n_genes', 'n_genes_by_counts', 'nFeature_RNA', 'n_genes_qc') if c in o.columns)
L = pd.read_csv(cl, sep='\t', usecols=['cell_id', 'r1_action', 'major_label'], low_memory=False).set_index('cell_id').reindex(o.index.astype(str))
c = o['leiden_r1'].astype(str).values; assert set(c) == set(R.index), pool
act = pd.Series(c).map(R['action']).values; lab = pd.Series(c).map(R['r2_label']).values
for kv in filter(None, os.environ.get('COLLAPSE', '').split(',')):
    pat, new = kv.split('='); lab = np.array([new if re.fullmatch(pat, x) else x for x in lab])
low = (L['r1_action'].values == 'keep_flag_low_depth') | (act == 'keep_flag_low_depth')
refact = {'keep', 'keep_tissue'} | set(filter(None, os.environ.get('EXTRA_REF_ACTIONS', '').split(',')))
ref = np.isin(act, list(refact)) & ~low
q = low & np.char.startswith(L['major_label'].values.astype(str), 'r3_') & ~np.isin(act, ['move_to_Myeloid', 'move_to_Endothelial', 'move_to_B_Plasma'])
def vote(zr, yr, zq):
    I = pynndescent.NNDescent(zr, n_neighbors=K + 1, random_state=seed, n_jobs=-1); I.prepare()
    i, d = I.query(zq, k=K); w = 1 / (d + 1e-6); w /= w.sum(1, keepdims=True)
    codes, names = pd.factorize(yr); s = np.zeros((len(zq), len(names)), np.float32)
    np.add.at(s, (np.arange(len(zq))[:, None].repeat(K, 1), codes[i]), w)
    return names[s.argmax(1)], s.max(1)
g = pd.to_numeric(o[genes_col]).values
vd = pd.DataFrame({'lab': lab[ref], 'g': g[ref], 'i': np.flatnonzero(ref)}); hold = vd.loc[vd.groupby('lab')['g'].rank(method='first', pct=True) <= .10, 'i'].values
keep = np.setdiff1d(np.flatnonzero(ref), hold)
p, cf = vote(Z[keep].astype(np.float32), lab[keep], Z[hold].astype(np.float32))
V = pd.DataFrame({'true': lab[hold], 'pred': p, 'conf': cf, 'genes': g[hold]}); V['confident'] = V.conf >= MIN_CONF; V['correct'] = V.true == V.pred
val = V.groupby('true').agg(cells=('pred', 'size'), median_genes=('genes', 'median'), accuracy=('correct', 'mean'), confident_share=('confident', 'mean')).round(3)
val['accuracy_if_confident'] = V[V.confident].groupby('true')['correct'].mean().round(3)
val.loc['ALL (macro)'] = [len(V), V.genes.median(), val.accuracy.mean().round(3), V.confident.mean().round(3), val.accuracy_if_confident.mean().round(3)]
val.to_csv(f'{out}.validation.tsv', sep='\t')
pq, cq = vote(Z[ref].astype(np.float32), lab[ref], Z[q].astype(np.float32))
T = pd.DataFrame({'cell_id': o.index[q].astype(str), 'pool': pool, 'r2_cluster': c[q], 'cluster_label': lab[q], 'cluster_action': act[q], 'transfer_label': pq, 'transfer_conf': cq.round(3)})
T['transfer_status'] = np.where(T.transfer_conf >= MIN_CONF, 'transferred', 'unresolved')
T.loc[T.transfer_status == 'unresolved', 'transfer_label'] = f'{pool}_unresolved_low_depth'
T.to_csv(f'{out}.cells.tsv.gz', sep='\t', index=False, compression='gzip')
S = T.groupby(['transfer_label', 'transfer_status']).size().rename('cells').reset_index().sort_values('cells', ascending=False); S.to_csv(f'{out}.summary.tsv', sep='\t', index=False)
idc = T[T.cluster_action.isin(refact)]
agr = (idc.cluster_label == idc.transfer_label).mean() if len(idc) else float('nan')
print(pool, 'ref', int(ref.sum()), 'query', int(q.sum()), 'unresolved', round((T.transfer_status == 'unresolved').mean(), 3), 'agreement with own identity cluster', round(agr, 3), 'n', len(idc)); print(val.to_string()); print(S.to_string(index=False))
