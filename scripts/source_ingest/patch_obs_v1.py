# Patch obs of pancancer_analysis_ready_v1.zarr in place (X untouched) and refresh composition tables.
# 2026-11-01 patch: GSE125449 patient_id = sample (P## numbers repeat across Set1/Set2 with different
# cancer types, so they are not patient identifiers; the paper reports 19 patients for 19 samples).
# Usage: patch_obs_v1.py ZARR REPORT_DIR
import sys, json
import pandas as pd, zarr
from anndata.io import read_elem, write_elem

z, rep = sys.argv[1:3]
g = zarr.open(z, 'r+')
obs = read_elem(g['obs'])
m = obs['source'].astype(str) == 'GSE125449'
pid = obs['patient_id'].astype(str)
pid[m] = 'GSE125449:' + obs.loc[m, 'sample_id'].astype(str)
obs['patient_id'] = pd.Categorical(pid)
del g['obs']; write_elem(g, 'obs', obs)
zarr.consolidate_metadata(z)
o = obs.reset_index(drop=True)
comp = lambda by: o.groupby(by, observed=True).agg(cells=('cell_id', 'size'), units=('unit_id', 'nunique'),
                                                   patients=('patient_id', 'nunique'), sources=('source', 'nunique')).reset_index()
comp(['data_batch']).to_csv(f'{rep}/composition_by_batch.tsv', sep='\t', index=False)
comp(['analysis_main']).sort_values('cells', ascending=False).to_csv(f'{rep}/composition_by_cancer.tsv', sep='\t', index=False)
comp(['source', 'data_batch', 'analysis_main']).to_csv(f'{rep}/composition_by_source.tsv', sep='\t', index=False)
for col in ['sampling_site', 'treatment_status', 'fresh_frozen', 'technology', 'capture_scope', 'qc_flag']:
    comp([col]).to_csv(f'{rep}/composition_by_{col}.tsv', sep='\t', index=False)
s = json.load(open(f'{rep}/summary.json')); s['patients'] = int(o.patient_id.nunique()); s['patch'] = 'GSE125449 patient_id = sample (2026-11-01)'
json.dump(s, open(f'{rep}/summary.json', 'w'), indent=1)
print(json.dumps(s), m.sum(), o.loc[o.source == 'GSE125449', 'patient_id'].nunique())
