# Annotation v2 (2026-11-07): apply the reviewed QC-round decisions (inventories/annotation_v2/lineage_r1_reviewed.tsv) to
# the cell-level major labels. Writes major_label = r2_<pool> for cells going to round 2, r1_<action> for cells that do not;
# labels outside the five pools are unchanged, except the tissue populations which become their own round-2 pools.
# Usage: apply_lineage_r1_review.py R1_DIR REVIEW_TSV CELL_MAJOR_LABELS OUT_CELLS OUT_SUMMARY
import sys
import pandas as pd, anndata as ad

r1, rev, cl, outc, outs = sys.argv[1:6]
R = pd.read_csv(rev, sep='\t', dtype={'cluster': str})
C = pd.read_csv(cl, sep='\t').set_index('cell_id')
C['label_major'] = C['major_label']
rows = []
for pool, d in R.groupby('pool'):
    o = ad.read_h5ad(f'{r1}/{pool}/object.h5ad', backed='r').obs
    c = o['leiden_r1'].astype(str)
    assert set(c) == set(d.cluster), pool
    act = c.map(d.set_index('cluster')['action']); lab = c.map(d.set_index('cluster')['r1_label'])
    new = act.map(lambda a: f'r2_{pool}' if a.startswith('keep') else ('r2_Myeloid' if a == 'move_to_Myeloid' else f'r1_{a}'))
    C.loc[o.index.astype(str), 'major_label'] = new.astype(str).values
    C.loc[o.index.astype(str), 'r1_label'] = lab.astype(str).values
    C.loc[o.index.astype(str), 'r1_action'] = act.astype(str).values
    rows.append(pd.DataFrame({'pool': pool, 'action': act.astype(str).values}).value_counts().rename('cells').reset_index())
for t in ('Kupffer', 'Microglia', 'Alveolar_macrophage', 'Sinusoidal_EC'):
    C.loc[C['major_label'] == t, 'major_label'] = f'r2_{t}'
C.to_csv(outc, sep='\t', compression='gzip')
S = pd.concat(rows).sort_values(['pool', 'cells'], ascending=[True, False]); S.to_csv(outs, sep='\t', index=False)
print(S.groupby('pool').cells.sum().to_dict()); print(S.to_string(index=False)); print(C['major_label'].value_counts().to_string())
