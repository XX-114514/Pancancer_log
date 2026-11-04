# Apply inventories/tme_v1/scope_exclusions_v1.tsv to tme_v1.zarr (2026-11-05): cells of excluded sources /
# units get obs['scope_excluded'] = reason and covarnet_v1_include = False. Cells are not removed from the
# zarr (it was integrated before the decision); the next TME integration leaves them out entirely.
# Usage: patch_tme_scope_v1.py TME_ZARR EXCLUSIONS_TSV
import sys, json
import numpy as np, pandas as pd, zarr
from anndata.io import read_elem, write_elem

zpath, exf = sys.argv[1:3]
ex = pd.read_csv(exf, sep='\t', dtype=str).fillna('')
g = zarr.open(zpath, 'r+')
obs = read_elem(g['obs'])
src, unit = obs['source'].astype(str).values, obs['unit_id'].astype(str).values
reason = np.array([''] * len(obs), dtype=object)
counts = {}
for _, r in ex.iterrows():
    m = (src == r['source']) if r['level'] == 'source' else (unit == r['unit_id'])
    reason[m] = r['reason'].split(' ')[0]
    counts[r['source'] + ('' if r['level'] == 'source' else '|' + r['unit_id'].split('_')[-1])] = int(m.sum())
obs['scope_excluded'] = pd.Categorical(reason)
obs['covarnet_v1_include'] = (obs['tme_flag'].astype(str) == 'non_malignant').values & (reason == '')
del g['obs']; write_elem(g, 'obs', obs)
uns = read_elem(g['uns'])
uns['scope_exclusions'] = ex.to_json(orient='records')
del g['uns']; write_elem(g, 'uns', uns)
zarr.consolidate_metadata(zpath)
print(json.dumps({'cells': int(len(obs)), 'excluded_cells_in_zarr': counts, 'covarnet_v1_include': int(obs['covarnet_v1_include'].sum())}))
