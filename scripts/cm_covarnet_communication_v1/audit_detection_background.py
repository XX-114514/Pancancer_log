"""Export aggregate background distributions; sample rows remain private."""
from pathlib import Path
import pandas as pd,numpy as np,zarr,json,tarfile,hashlib
from anndata.io import read_elem
r=Path(__file__).resolve().parents[1];q=pd.read_pickle(r/'private/cell_metadata.pkl.gz');z=zarr.open(str(r.parent/'20260818_pancancer_5M_merge_prep/02_merge_new_cohort/objects/5M_merge_concat.zarr'),mode='r');obs=read_elem(z['obs']);assert obs.index.equals(q.index)
valid=q.CMT.notna();a=pd.DataFrame({'sampleID':q.sampleID,'CMT':q.CMT,'identity':q.final_annotation_with_malignancy_v6.astype(str),'source_QC_total_counts':obs.total_counts,'source_QC_n_genes':obs.n_genes_by_counts});a=a[valid]
s=a.groupby('sampleID',observed=True).agg(CMT=('CMT','first'),n_cells=('CMT','size'),median_source_QC_counts=('source_QC_total_counts','median'),median_source_QC_genes=('source_QC_n_genes','median'));s.to_csv(r/'private/sample_detection_eligibility.tsv',sep='\t');out=[]
for cm in sorted(s.CMT.unique()):
 for group,mask in [('High',s.CMT==cm),('Low',s.CMT!=cm)]:
  for metric in ['n_cells','median_source_QC_counts','median_source_QC_genes']:
   v=s.loc[mask,metric].dropna();out.append(dict(cm=cm,group=group,metric=metric,n_samples=len(v),q25=v.quantile(.25),median=v.median(),q75=v.quantile(.75),note='source_QC_fields_may_use_original_gene_panels; descriptive_not_adjusted_effect'))
pd.DataFrame(out).to_csv(r/'audit/HighLow_detection_background.tsv',sep='\t',index=False)
counts=a.groupby(['CMT','sampleID','identity'],observed=True).size().rename('n_cells').reset_index();counts.to_csv(r/'private/sample_identity_cell_counts.tsv.gz',sep='\t',index=False);members=pd.read_csv(r/'audit/network_members.tsv',sep='\t');rows=[]
for row in members.itertuples():
 for group,mask in [('High',s.CMT==row.cm),('Low',s.CMT!=row.cm)]:
  ids=s.index[mask];v=counts[counts.identity==row.subCluster].set_index('sampleID').n_cells.reindex(ids,fill_value=0);pos=v[v>0];rows.append(dict(cm=row.cm,identity=row.subCluster,group=group,n_samples=len(v),n_identity_absent=int((v==0).sum()),n_identity_lt20=int(((v>0)&(v<20)).sum()),n_identity_ge20=int((v>=20).sum()),median_cells_when_present=pos.median(),q25_cells_when_present=pos.quantile(.25),q75_cells_when_present=pos.quantile(.75)))
pd.DataFrame(rows).to_csv(r/'audit/HighLow_member_cell_eligibility.tsv',sep='\t',index=False)
files=[r/'audit/HighLow_detection_background.tsv',r/'audit/HighLow_member_cell_eligibility.tsv']+list((r/'audit/metadata_v2').glob('CM*.tsv'))
with tarfile.open(r/'background_review_tables.tar.gz','w:gz') as tar:
 for p in files:tar.add(p,arcname=p.name,recursive=False)
print(json.dumps({'rows':len(out),'member_group_rows':len(rows),'bytes':(r/'background_review_tables.tar.gz').stat().st_size,'sha256':hashlib.sha256((r/'background_review_tables.tar.gz').read_bytes()).hexdigest()}))
