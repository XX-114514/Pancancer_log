import json
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');t=pd.read_pickle(r/'stage_c4_harmonized_identity_freeze/annotation_all_cells_with_release_status_PRIVATE.pkl.gz');assert t.index.equals(q.index)
take=q.final_annotation_v7.astype(str).isin(['Bcell_coarse_provisional','Bcell_coarse_subtype_pending'])&q.current_QC_keep;rows=[]
for label,g in q.loc[take].groupby('final_annotation_v7',observed=True):
 z=t.loc[g.index];d=z.groupby(['annotation_frozen','final_parent','final_annotation','freeze_status'],observed=True,dropna=False).size().rename('n_cells').reset_index();d.insert(0,'historical_label',str(label));rows.append(d)
pd.concat(rows,ignore_index=True).to_csv(r/'stage_c4_harmonized_identity_freeze/B_previous_coarse_to_frozen_identity.tsv',sep='\t',index=False,na_rep='')
print(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'n_B_previous_coarse_QC_keep':int(take.sum()),'n_new_frozen':int(t.loc[take,'annotation_frozen'].sum()),'n_outside_new_freeze':int((~t.loc[take,'annotation_frozen']).sum())}))
