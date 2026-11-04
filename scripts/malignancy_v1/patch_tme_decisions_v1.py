# Apply the 2026-11-04 decisions (inventories/tme_v1/source_decisions_v1.tsv) to tme_v1.zarr obs/uns:
#  - covarnet_v1_include: True only for tme_flag == non_malignant (uncertain cells set aside from CoVarNet v1)
#  - fresh_frozen of GSE189903 -> cryopreserved_or_frozen (paper Methods: suspension cryopreserved, thawed
#    before 10x 5'); original value kept in fresh_frozen_v2
#  - uns excluded_sources: GSE179994 reason corrected to author_T_cell_subset_only
# Usage: patch_tme_decisions_v1.py TME_ZARR
import sys, json
import numpy as np, pandas as pd, zarr
from anndata.io import read_elem, write_elem

zpath = sys.argv[1]
g = zarr.open(zpath, 'r+')
obs = read_elem(g['obs'])
obs['covarnet_v1_include'] = (obs['tme_flag'].astype(str) == 'non_malignant').values
obs['fresh_frozen_v2'] = obs['fresh_frozen'].astype(str).values
ff = obs['fresh_frozen'].astype(str).values.copy()
ff[obs['source'].astype(str).values == 'GSE189903'] = 'cryopreserved_or_frozen'
obs['fresh_frozen'] = pd.Categorical(ff)
del g['obs']; write_elem(g, 'obs', obs)
uns = read_elem(g['uns'])
uns['excluded_sources'] = json.dumps({'GSE179994': 'author_T_cell_subset_only'})
uns['decisions'] = 'inventories/tme_v1/source_decisions_v1.tsv (2026-11-04): uncertain cells kept but covarnet_v1_include=False; GSE189903 fresh_frozen corrected'
del g['uns']; write_elem(g, 'uns', uns)
zarr.consolidate_metadata(zpath)
print(json.dumps({'cells': int(len(obs)), 'covarnet_v1_include': int(obs['covarnet_v1_include'].sum()),
                  'GSE189903_cells_fixed': int((obs['source'].astype(str) == 'GSE189903').sum())}))
