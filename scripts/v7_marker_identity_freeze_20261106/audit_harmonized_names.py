"""Exact comparison to the independently audited full-axis membership, plus all renamed values."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';c=r/'stage_c3_exact_identity_freeze';o=r/'stage_c4_harmonized_identity_freeze';assert json.loads((c/'independent_freeze_audit.json').read_text())['status']=='PASS'
a=pd.read_pickle(c/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');b=pd.read_pickle(o/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');assert a.index.equals(b.index) and a.columns.equals(b.columns)
for col in a:
 if col in ['final_annotation','representative_markers','release_version']:continue
 x=a[col].astype(object).to_numpy();y=b[col].astype(object).to_numpy();na=pd.isna(x);assert np.array_equal(na,pd.isna(y)) and np.array_equal(x[~na],y[~na]),col
reg=pd.read_csv(o/'harmonized_name_registry.tsv',sep='\t').set_index('accepted_rule_id');take=b.annotation_frozen.to_numpy()
for j,row in reg.iterrows():
 mask=take&b.accepted_rule_id.eq(j).fillna(False).to_numpy(bool)
 for col in ['final_parent','final_annotation','representative_markers']:assert b.loc[mask,col].astype(str).eq(row[col]).all()
assert b.final_annotation.notna().to_numpy().tolist()==take.tolist();assert not b.loc[take,'final_annotation'].astype(str).str.contains('unresolved|coarse|provisional|pending|review|candidate|boundary|unknown|ambiguous',case=False).any()
f=pd.read_pickle(o/'marker_identity_frozen_cells_PRIVATE.pkl.gz');expected=b.loc[take];assert f.index.equals(expected.index)
for col in f:
 x=f[col].astype(object).to_numpy();y=expected[col].astype(object).to_numpy();na=pd.isna(x);assert np.array_equal(na,pd.isna(y)) and np.array_equal(x[~na],y[~na]),col
for name,total,count in [('frozen_identity_counts.tsv',len(f),'n_cells'),('frozen_group_counts.tsv',len(f),'n_frozen_cells'),('name_consolidation_audit.tsv',len(f),'n_cells'),('all_cell_release_status_counts.tsv',len(b),'n_cells'),('historical_to_new_identity_counts.tsv',len(b),'n_cells')]:assert pd.read_csv(o/name,sep='\t')[count].sum()==total
position=0
for chunk in pd.read_csv(o/'marker_identity_frozen_cells_PRIVATE.tsv.gz',sep='\t',usecols=['cell_id','final_annotation','final_parent'],chunksize=100000):
 n=len(chunk);v=f.iloc[position:position+n];assert np.array_equal(chunk.cell_id.to_numpy(),v.index.to_numpy())
 for col in ['final_annotation','final_parent']:assert np.array_equal(chunk[col].to_numpy(),v[col].astype(str).to_numpy())
 position+=n
assert position==len(f)
for z in json.loads((o/'complete.json').read_text())['files']:
 p=o/z['name'];assert p.stat().st_size==z['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256']
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_all_axis':len(b),'n_frozen':len(f),'n_labels':f.final_annotation.nunique(),'n_parents':f.final_parent.nunique(),'all_membership_QC_CNV_state_parent_rule_and_provenance_fields_exactly_unchanged':True,'all_names_match_harmonized_registry':True,'all_private_TSV_rows_and_label_values_verified':True,'all_artifact_hashes_verified':True,'all_atlas_publication_acceptance':False};(o/'independent_freeze_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
