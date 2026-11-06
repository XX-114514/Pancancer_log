"""Retain explicit grouped reasons for QC-kept cells outside the formal identity release."""
import json
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';c=r/'stage_c4_harmonized_identity_freeze';o=r/'stage_g_outside_identity_diagnostics';o.mkdir(exist_ok=False)
assert json.loads((c/'independent_freeze_audit.json').read_text())['status']=='PASS'
t=pd.read_pickle(c/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');assert t.index.equals(q.index);policy=json.loads((r/'identity_policy_v2.json').read_text());mask=t.QC_keep&~t.annotation_frozen;z=q.loc[mask].copy();z['freeze_status']=t.loc[mask,'freeze_status'];rows=[]
for (group,status),v in z.groupby(['audit_group','freeze_status'],observed=True):
 record={'audit_group':group,'outside_reason':status,'n_cells':len(v),'n_sources':v.gse_id.nunique(),'historical_parent_mode':v.major_celltype_v7.mode().iloc[0],'previous_working_parent_mode':v.current_parent.mode().iloc[0]}
 for j,p in enumerate(policy['rules']):record['positive_rule_'+str(j)+'_'+p['parent']]=int(v['positive_'+str(j)].sum())
 rows.append(record)
df=pd.DataFrame(rows).sort_values('n_cells',ascending=False);assert df.n_cells.sum()==749597;df.to_csv(o/'outside_identity_group_diagnostics.tsv',sep='\t',index=False)
out={'status':'PASS_outside_identity_accounting','UTC':datetime.now(timezone.utc).isoformat(),'n_QC_keep_outside_identity_release':int(mask.sum()),'n_group_reason_records':len(df),'counts_conserved':True,'positive_programmes_are_overlapping_not_cell_type_proportions':True,'not_new_QC_or_doublet_exclusions':True};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
