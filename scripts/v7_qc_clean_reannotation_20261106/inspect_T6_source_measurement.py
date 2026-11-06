import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import json
import numpy as np,pandas as pd,h5py
from anndata.io import read_elem
r=PROJECT_ROOT;run=r/'runs/20261106_v7_QC_clean_reannotation_v1';q=pd.read_pickle(run/'stage_c_clean_views/T6_cleaning_decisions_PRIVATE.pkl.gz');rows=[]
for gs,t in q.groupby('gse_id',observed=True):
 p=r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(str(gs)+'_symbol_counts.h5ad')
 with h5py.File(p,'r') as f:vv=read_elem(f['var']);var=pd.Index(vv.index.astype(str))
 cd3=[g for g in ['CD3D','CD3E','CD3G'] if g in var];tcr=[g for g in ['TRAC','TRBC1','TRBC2'] if g in var]
 if not cd3 or not tcr:
  rows.append({'gse_id':str(gs),'n_T6':len(t),'n_kept':int(t.primary_QC_keep.sum()),'n_source_features':len(var),'CD3_available':','.join(cd3) or 'none','TCR_available':','.join(tcr) or 'none','FOXP3_available':'FOXP3' in var,'n_ENSG_like_features':int(var.str.match('^ENSG').sum()),'var_columns':list(vv.columns),'CD3_alias_candidates':[g for g in var if g.upper().startswith(('CD3D','CD3E','CD3G','TRAC','TRBC'))][:20]})
out={'scope':'Source feature-axis availability; not a verified absence of RNA or immune identity.','records':rows,'n_T6_marker_component_missing':int((~q.T_parent_panel_evaluable).sum())};(run/'stage_c_clean_views/T6_source_measurement_diagnostic.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
