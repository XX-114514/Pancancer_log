# Annotation v2 (2026-11-07): apply the reviewed round-2 decisions (inventories/annotation_v2/lineage_r2_reviewed.tsv).
# Writes major_label = r3_<pool> for cells that stay or move into a pool, r2_<action> for cells that leave. Round-3 reruns
# (Myeloid, Endothelial) select their cells from this file. origin_pool, r2_cluster and r2_label keep where each cell came from.
# Usage: apply_lineage_r2_review.py R2_DIR REVIEW_TSV CELL_LABELS_R2 OUT_CELLS OUT_SUMMARY
import sys
import pandas as pd, h5py
from anndata.io import read_elem

r2, rev, cl, outc, outs = sys.argv[1:6]
R = pd.read_csv(rev, sep='\t', dtype={'cluster': str})
C = pd.read_csv(cl, sep='\t').set_index('cell_id')
for k in ('origin_pool', 'r2_cluster', 'r2_label', 'r2_action'):
    C[k] = pd.NA
rows = []
for pool, d in R.groupby('pool'):
    with h5py.File(f'{r2}/{pool}/object.h5ad', 'r') as f:
        c = read_elem(f['obs'])['leiden_r1'].astype(str)
    assert set(c) == set(d.cluster), pool
    D = d.set_index('cluster'); act = c.map(D['action']); lab = c.map(D['r2_label'])
    new = act.map(lambda a: f'r3_{a[8:]}' if a.startswith('move_to_') else (f'r3_{pool}' if a.startswith('keep') else f'r2_{a}'))
    i = c.index.astype(str)
    assert (C.loc[i, 'major_label'] == f'r2_{pool}').all(), pool
    C.loc[i, 'major_label'] = new.values; C.loc[i, 'origin_pool'] = pool; C.loc[i, 'r2_cluster'] = c.values
    C.loc[i, 'r2_label'] = lab.values; C.loc[i, 'r2_action'] = act.values
    rows.append(pd.DataFrame({'pool': pool, 'action': act.values, 'major_label': new.values}).value_counts().rename('cells').reset_index())
assert not C['major_label'].astype(str).str.match(r'^r2_(T_NK|B_Plasma|Myeloid|Stromal|Endothelial|Kupffer|Microglia|Alveolar_macrophage|Sinusoidal_EC)$').any()
C.to_csv(outc, sep='\t', compression='gzip')
S = pd.concat(rows).sort_values(['pool', 'cells'], ascending=[True, False]); S.to_csv(outs, sep='\t', index=False)
print(S.to_string(index=False)); print(C['major_label'].value_counts().to_string())
