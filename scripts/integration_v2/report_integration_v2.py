# Composition tables and read-back checks for integration_v2.zarr.
# Usage: report_integration_v2.py RUN_DIR REPORT_OUT
import sys, json
import numpy as np, pandas as pd, zarr
from anndata.io import read_elem

run, rep = sys.argv[1:3]
g = zarr.open(f'{run}/integration_v2.zarr', 'r')
obs = read_elem(g['obs']).reset_index(drop=True)
chk = {'n_obs': len(obs), 'cell_id_unique': bool(obs['cell_id'].is_unique),
       'X_shape': g['X'].attrs['shape'], 'layers': list(g['layers'].keys()), 'obsm': {k: list(g['obsm'][k].shape) for k in g['obsm']},
       'umap_finite': bool(np.isfinite(g['obsm']['X_umap'][:]).all()),
       'counts_integer_sample': bool(np.all(np.mod(g['layers']['counts']['data'][:2_000_000], 1) == 0)),
       'n_selected_HVG': int(read_elem(g['var'])['selected_HVG'].sum())}
json.dump(chk, open(f'{rep}/readback_v2.json', 'w'), indent=1, default=str)
def comp(by):
    return obs.groupby(by, observed=True).agg(cells=('cell_id', 'size'), units=('unit_id', 'nunique'), patients=('patient_id', 'nunique'), sources=('source', 'nunique')).reset_index()
comp(['data_batch']).to_csv(f'{rep}/composition_by_batch.tsv', sep='\t', index=False)
comp(['analysis_main']).sort_values('cells', ascending=False).to_csv(f'{rep}/composition_by_cancer.tsv', sep='\t', index=False)
comp(['source', 'data_batch', 'analysis_main']).to_csv(f'{rep}/composition_by_source.tsv', sep='\t', index=False)
for c in ['sampling_site', 'treatment_status', 'fresh_frozen', 'technology', 'capture_scope', 'qc_flag', 'oncotree_code']:
    comp([c]).to_csv(f'{rep}/composition_by_{c}.tsv', sep='\t', index=False)
c = obs['decontX_contamination'].astype(float)
summ = {'cells': len(obs), 'sources': int(obs.source.nunique()), 'units': int(obs.unit_id.nunique()), 'patients': int(obs.patient_id.nunique()),
        'cancers': int(obs.analysis_main.nunique()), 'decontX_contamination_median': float(c.median()), 'decontX_contamination_p90': float(c.quantile(0.9))}
json.dump(summ, open(f'{rep}/summary.json', 'w'), indent=1)
print(json.dumps(chk, default=str)); print(json.dumps(summ))
