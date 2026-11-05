"""Read-only V7 audit with new aggregate outputs; no annotation overwrite.
Environment: PROJECT_ROOT required. All cell-level metadata stays on server.
"""
import os,json,hashlib,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))))
assert ROOT.is_absolute()
RUN=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1'
OUT=RUN/'stage_a_frozen_audit'
OUT.mkdir(parents=True,exist_ok=False)
START=datetime.now(timezone.utc).isoformat()
B=ROOT/'runs/20260818_pancancer_5M_merge_prep'
V=B/'11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1'
F=B/'07_fine_annotation/20260823_corrected_evidence_v1'
U=B/'07_fine_annotation/20260905_uncertainty_resolution_v1'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8*1024**2),b''):h.update(x)
 return h.hexdigest()
def save(d,name):
 p=OUT/name;d.to_csv(p,index=False,sep='\t');q=pd.read_csv(p,sep='\t');assert len(q)==len(d)
 return {'path':str(p.relative_to(ROOT)),'rows':len(d),'sha256':sha(p)}
def counts(s):return json.dumps({str(k):int(v) for k,v in s.value_counts(dropna=False).items()},ensure_ascii=False,sort_keys=True)
p=V/'checkpoints/final_annotation_v7.pkl.gz';fingerprint=sha(p)
assert fingerprint=='e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e'
df=pd.read_pickle(p);assert len(df)==4676787 and df.index.is_unique
files=[];records=[]
for (major,label),g in df.groupby(['major_celltype_v7','final_annotation_v7'],observed=True):
 records.append({'major':major,'label_v7':label,'n_cells':len(g),'n_gse':g.gse_id.nunique(),
 'n_samples':g.groupby(['gse_id','sampleID'],observed=True).ngroups,
 'n_strict_major':int(g.strict_major_eligible_v7.sum()),'n_strict_subtype':int(g.strict_subtype_eligible_v7.sum()),
 'n_exclude_primary_major':int(g.exclude_from_primary_major_v2.sum()),
 'n_requires_recluster':int(g.requires_recluster_for_subtype_v2.sum()),
 'confidence_counts':counts(g.annotation_confidence_v7),
 'major_reason_counts':counts(g.major_resolution_reason_v2),
 'subtype_method_counts':counts(g.subtype_resolution_method_v2),
 'legacy_cluster_counts':counts(g.cluster_uid),'celltypist_major_counts':counts(g.celltypist_major),
 'celltypist_confidence_median':float(g.celltypist_confidence.median()),
 'old_parent_counts':counts(g.major_celltype_v5)})
labels=pd.DataFrame(records);assert labels.n_cells.sum()==len(df)
files.append(save(labels,'v7_label_reason_audit.tsv'))
major=labels.groupby('major',as_index=False)[['n_cells','n_strict_major','n_strict_subtype','n_exclude_primary_major','n_requires_recluster']].sum()
files.append(save(major,'v7_major_reason_audit.tsv'))
ct=df.groupby(['major_celltype_v7','final_annotation_v7','cluster_uid','subtype_resolution_method_v2','requires_recluster_for_subtype_v2'],dropna=False,observed=True).size().reset_index(name='n_cells')
assert ct.n_cells.sum()==len(df)
files.append(save(ct,'v7_label_legacy_cluster_crosswalk.tsv'))
for name,src in [('legacy_initial_cluster_gates.tsv',F/'tables/cluster_annotation_evidence_v5.tsv'),('legacy_expanded_cluster_gates.tsv',U/'tables/expanded_panel_cluster_rescue_v1.tsv'),('legacy_descriptive_cluster_evidence.tsv',U/'tables/de_novo_unresolved_labels_v1.tsv')]:
 files.append(save(pd.read_csv(src,sep='\t'),name))
# Wide legacy panel is log1p CP10K normalized against each source full library.
# It is used to diagnose shared identity vs subtype DE, never for new HVG/PCA.
a=ad.read_h5ad(F/'objects/wide_marker_panel_full_library_masked_v2.h5ad')
assert a.obs_names.equals(df.index)
genes=a.var_names.astype(str).to_numpy();assert len(genes)==359
X=a.X.tocsr(copy=False);nan_count=int(np.isnan(X.data).sum());X.data=np.nan_to_num(X.data,nan=0.0);X.eliminate_zeros()
assert np.isfinite(X.data).all() and X.data.min(initial=0)>=0
availability=np.asarray(a.varm['availability_by_gse']).astype(bool)
gses=np.asarray(a.uns['availability_gse_order']).astype(str)
assert availability.shape==(len(genes),len(gses))
gidx=pd.Index(gses).get_indexer(df.gse_id.astype(str));assert (gidx>=0).all()
# Each cell has one final-label group and one legacy-cluster group, kept separate.
meta=pd.DataFrame({'kind':'label','group':df.final_annotation_v7.astype(str).to_numpy()})
label_codes,label_names=pd.factorize(meta['group'],sort=True)
cluster_codes,cluster_names=pd.factorize(df.cluster_uid.fillna('NO_LEGACY_CLUSTER').astype(str),sort=True)
all_marker=[];coexpr=[]
programs={'T_identity':['CD3D','CD3E','TRAC'],'NK_cytotoxic_shared':['NKG7','GNLY','KLRD1'],
 'Mast_identity':['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],
 'B_identity':['MS4A1','CD79A','CD79B'],'Myeloid_shared':['LYZ','LST1','TYROBP'],
 'Erythroid_non_HB':['AHSP','GYPA','KLF1','SLC4A1'],
 'Platelet_identity':['PF4','PPBP','GP9','ITGA2B'],'HSPC_candidate':['CD34','SPINK2','GATA2','LMO2'],
 'Epithelial_identity':['EPCAM','KRT8','KRT18','KRT19']}
program_states={}
for name,panel in programs.items():
 ix=np.flatnonzero(np.isin(genes,panel));det=np.asarray((X[:,ix]>0).sum(axis=1)).ravel()
 measured=availability[ix].sum(axis=0)[gidx]>=2
 program_states[name]=(det>=2)&measured,measured
for kind,codes,names in [('label',label_codes,label_names),('legacy_cluster',cluster_codes,cluster_names)]:
 n=len(names);H=sp.csr_matrix((np.ones(len(df),dtype=np.int32),(codes,np.arange(len(df)))),shape=(n,len(df)))
 group_n=np.bincount(codes,minlength=n);sg=sp.coo_matrix((np.ones(len(df),dtype=np.int32),(codes,gidx)),shape=(n,len(gses))).tocsr()
 denom=np.asarray(sg@availability.T.astype(np.int32));sums=(H@X).toarray();Z=X.copy();Z.data[:]=1;detected=(H@Z).toarray()
 linear=X.copy();linear.data=np.expm1(linear.data);linear_sums=(H@linear).toarray()
 for i,name in enumerate(names):
  for j,gene in enumerate(genes):
   den=int(denom[i,j]);all_marker.append({'group_kind':kind,'group':name,'gene':gene,'n_cells':int(group_n[i]),'n_measured':den,
    'n_expressing':int(round(detected[i,j])),'pct_expressing_measured':float(detected[i,j]/den) if den else np.nan,
    'mean_log1p_cp10k_measured':float(sums[i,j]/den) if den else np.nan,
    'mean_cp10k_measured':float(linear_sums[i,j]/den) if den else np.nan})
 for name,(positive,evaluable) in program_states.items():
  nn=np.bincount(codes,weights=positive,minlength=n).astype(int);ne=np.bincount(codes,weights=evaluable,minlength=n).astype(int)
  for i,g in enumerate(names):coexpr.append({'group_kind':kind,'group':g,'program':name,'n_cells':int(group_n[i]),'n_evaluable':int(ne[i]),'n_ge2_markers':int(nn[i]),'fraction_ge2_among_evaluable':float(nn[i]/ne[i]) if ne[i] else np.nan})
 for pa,pb in [('T_identity','NK_cytotoxic_shared'),('T_identity','Mast_identity'),('T_identity','Myeloid_shared'),('Mast_identity','Myeloid_shared'),('Erythroid_non_HB','Myeloid_shared')]:
  positive=program_states[pa][0]&program_states[pb][0];evaluable=program_states[pa][1]&program_states[pb][1]
  nn=np.bincount(codes,weights=positive,minlength=n).astype(int);ne=np.bincount(codes,weights=evaluable,minlength=n).astype(int)
  for i,g in enumerate(names):coexpr.append({'group_kind':kind,'group':g,'program':pa+'__AND__'+pb,'n_cells':int(group_n[i]),'n_evaluable':int(ne[i]),'n_ge2_markers':int(nn[i]),'fraction_ge2_among_evaluable':float(nn[i]/ne[i]) if ne[i] else np.nan})
files.append(save(pd.DataFrame(all_marker),'legacy_panel_expression_audit.tsv'))
files.append(save(pd.DataFrame(coexpr),'legacy_panel_same_cell_program_audit.tsv'))
manifest={'stage':'A_completed_diagnostic_only','UTC_start':START,'UTC_end':datetime.now(timezone.utc).isoformat(),'input_sha256':fingerprint,
 'n_cells':len(df),'n_final_labels':len(labels),'n_major':len(major),'original_checkpoint_unchanged':sha(p)==fingerprint,
 'panel_axis_exact':True,'n_panel_genes':len(genes),'panel_nan_entries_replaced_for_aggregate_with_measured_denominator':nan_count,
 'outputs':files,'annotation_changes':0,
 'limits':['Aggregate >=2-marker coexpression is descriptive evidence, not a doublet or identity classifier.','Legacy panel does not replace full-gene lineage refitting.','Measured-gene denominator uses explicit per-source availability; absent genes are not biological negatives.','Frozen malignancy and original annotation fields are unchanged.']}
(OUT/'complete.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest))
