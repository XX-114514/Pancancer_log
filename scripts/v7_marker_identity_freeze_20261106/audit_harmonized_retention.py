import json
from datetime import datetime,timezone
import numpy as np,pandas as pd,zarr
from anndata.io import read_elem
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';c=r/'stage_c4_harmonized_identity_freeze';o=r/'stage_e4_release_retention';o.mkdir(exist_ok=False);assert json.loads((c/'independent_freeze_audit.json').read_text())['status']=='PASS'
q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');t=pd.read_pickle(c/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');assert q.index.equals(t.index);mask=q.current_QC_keep.to_numpy();z=pd.DataFrame({'gse_id':q.gse_id.astype(str),'historical_parent':q.major_celltype_v7.astype(str),'included':t.annotation_frozen.to_numpy(),'n_cells':np.ones(len(q),np.int8)},index=q.index)
g=zarr.open_group(str(PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');assert pd.Index(read_elem(g['obs/_index'])).equals(q.index)
if 'cancer_type' in g['obs']:z['historical_cancer_type']=np.asarray(read_elem(g['obs/cancer_type'])).astype(str)
files=[]
for keys,name in [(['gse_id'],'source_retention.tsv'),(['historical_parent'],'historical_parent_retention.tsv'),(['gse_id','historical_parent'],'source_by_historical_parent_retention.tsv')]+([(['historical_cancer_type'],'historical_cancer_context_retention.tsv')] if 'historical_cancer_type' in z else []):
 b=z.loc[mask].groupby(keys,observed=True,dropna=False).agg(n_QC_keep_before_identity_release=('n_cells','sum'),n_identity_frozen=('included','sum')).reset_index();b['n_not_in_identity_freeze']=b.n_QC_keep_before_identity_release-b.n_identity_frozen;b['fraction_identity_frozen']=b.n_identity_frozen/b.n_QC_keep_before_identity_release
 assert b.n_QC_keep_before_identity_release.sum()==4668312 and b.n_identity_frozen.sum()==t.annotation_frozen.sum();b.to_csv(o/name,sep='\t',index=False);files.append(name)
out={'status':'PASS_retention_accounting','UTC':datetime.now(timezone.utc).isoformat(),'n_QC_keep':int(mask.sum()),'n_frozen':int(t.annotation_frozen.sum()),'n_sources':int(z.loc[mask,'gse_id'].nunique()),'files':files,'meaning':'Selection retention only; not unbiased population abundances or source/identity final publication acceptance. No sample or patient rows exported.'};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
