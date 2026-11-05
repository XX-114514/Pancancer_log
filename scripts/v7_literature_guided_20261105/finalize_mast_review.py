"""New evidence/review overlay only; original/final V7 labels are preserved."""
import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,h5py
from anndata.io import read_elem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))))
r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';src=r/'stage_c2_mast_initial_diagnostic';o=r/'stage_f_mast_review';o.mkdir(exist_ok=False)
terminal=json.loads((r/'launch/mast_initial_v2_terminal.json').read_text());assert terminal['exit_code']==0
md=pd.read_csv(src/'mast_initial_metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id')
a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');assert md.index.equals(a.obs_names)
cl=md['RNA_snn_res.0.4'].astype(int).to_numpy();assert set(cl)==set(range(27))
qc=pd.read_csv(src/'mast_initial_cluster_QC_res0.4.tsv',sep='\t');source=pd.read_csv(src/'cluster_source_concentration.tsv',sep='\t')
markers=pd.read_csv(src/'mast_initial_RNA_FindAllMarkers_res0.4.tsv',sep='\t')
# Read actual source feature headers: negatives must not use the full union denominator.
gses=sorted(a.obs.gse_id.astype(str).unique());av=[];header_manifest=[]
for gs in gses:
 p=ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad')
 assert p.exists(),gs
 with h5py.File(p,'r') as f:
  vg=read_elem(f['var']);native=pd.Index(vg.index.astype(str));assert native.is_unique
 av.append(a.var_names.isin(native));header_manifest.append({'source':gs,'n_features':len(native),'n_core_features_measured':int(av[-1].sum())})
av=np.asarray(av).T;assert av.shape==(a.n_vars,len(gses))
gidx=pd.Index(gses).get_indexer(a.obs.gse_id.astype(str));H=sp.csr_matrix((np.ones(a.n_obs,np.int32),(cl,np.arange(a.n_obs))),shape=(27,a.n_obs))
S=sp.csr_matrix((np.ones(a.n_obs,np.int32),(cl,gidx)),shape=(27,len(gses)))
panel=['TPSAB1','TPSB2','CPA3','KIT','MS4A2','HDC','LTC4S','IL1RL1','CMA1','CTSG','CD3D','CD3E','TRAC','LYZ','LST1','TYROBP','EPCAM','KRT8','KRT18','KRT19','COL1A1','COL1A2','DCN','MKI67','TOP2A','UBE2C','BIRC5','ISG15','IFIT1','FOS','JUN']
j=a.var_names.get_indexer(panel);assert (j>=0).all()
X=a.X[:,j].tocsr();det=X.copy();det.data[:]=1
nd=(H@det).toarray();den=np.asarray(S@av[j].T.astype(np.int32));tot=np.asarray(a.X.sum(axis=1)).ravel();norm=X.multiply(10000/np.maximum(tot,1)[:,None]).tocsr();linear=(H@norm).toarray();n=np.bincount(cl,minlength=27)
records=[]
for ci in range(27):
 for jj,gene in enumerate(panel):
  dd=int(den[ci,jj]);records.append({'cluster':ci,'gene':gene,'n_cells':int(n[ci]),'n_measured':dd,'coverage_fraction':float(dd/n[ci]),'pct_expressing_measured':float(nd[ci,jj]/dd) if dd else None,'mean_cp10k_measured':float(linear[ci,jj]/dd) if dd else None})
e=pd.DataFrame(records);e.to_csv(o/'mast27_measured_marker_expression.tsv',sep='\t',index=False)
# Manual diagnostic dispositions based on RNA marker programs + QC + source data.
# These are review routes, not automatic doublet calls or final cell identities.
held={7:'epithelial_program_and_source_concentration_review',10:'T_plus_Mast_program_with_elevated_doublet_score_review',13:'stromal_epithelial_program_high_decontx_review',15:'myeloid_plus_Mast_program_high_doublet_and_decontx_review',17:'other_lineage_program_high_decontx_review',25:'keratin_program_high_decontx_review',26:'prostate_epithelial_program_high_decontx_review'}
top=markers[markers.avg_log2FC>0].sort_values(['cluster','avg_log2FC'],ascending=[True,False]).groupby('cluster').head(10)
rows=[]
for ci in range(27):
 genes=top.loc[top.cluster==ci,'gene'].tolist()
 rows.append({'cluster':ci,'working_parent_review':'Mast_competing_program_review' if ci in held else 'Mast_parent_supported',
 'working_subtype_recommendation':'hold_for_identity_QC_review' if ci in held else 'coarse_Mast_endpoint_acceptable',
 'state_program_candidate':'cycling_enriched' if ci==16 else 'not_assigned',
 'reason':held.get(ci,'Mast_identity_markers_supported; source_feature_or_state_groups_not_promoted_to_subtypes'),
 'top10_RNA_DEG':','.join(genes),'n_top10_features_measured_in_le2_Mast_sources':sum(int(av[a.var_names.get_loc(g)].sum())<=2 for g in genes if g in a.var_names)})
review=pd.DataFrame(rows).merge(qc,on='cluster',validate='one_to_one').merge(source,on='cluster',validate='one_to_one');assert review.n_cells.sum()==45036
review.to_csv(o/'mast27_cluster_review.tsv',sep='\t',index=False)
# Minimal new cell-level review fields, never replacing any frozen label.
overlay=a.obs[['major_celltype_v7','final_annotation_v7','annotation_confidence_v7']].copy();overlay['mast_review_cluster_res0.4']=cl
for c in ['working_parent_review','working_subtype_recommendation','state_program_candidate']:overlay[c]=review.set_index('cluster')[c].reindex(cl).to_numpy()
overlay['review_version']='V7_literature_guided_20261105_v1';overlay.to_pickle(o/'mast_cell_review_overlay_PRIVATE.pkl.gz')
# Source-specific genes are annotated, not silently removed/renamed/merged.
markers['n_Mast_sources_measuring_exact_feature']=markers.gene.map({g:int(av[k].sum()) for k,g in enumerate(a.var_names)})
markers.to_csv(o/'mast_initial_DEG_with_measurement_scope.tsv',sep='\t',index=False)
(o/'source_header_manifest_PRIVATE.json').write_text(json.dumps(header_manifest,indent=2))
# Scientific dot plot: positive measured expression only. Low coverage is shown as x.
show=['TPSAB1','TPSB2','CPA3','KIT','MS4A2','CD3D','CD3E','TRAC','LYZ','LST1','EPCAM','KRT18','COL1A2','DCN','MKI67','TOP2A']
p=e[e.gene.isin(show)].copy();p['x']=p.gene.map({g:i for i,g in enumerate(show)});p['v']=np.log1p(p.mean_cp10k_measured)
p['z']=p.groupby('gene')['v'].transform(lambda s:(s-s.mean())/(s.std(ddof=0) or 1)).clip(-2,2)
fig,ax=plt.subplots(figsize=(12,11));good=(p.coverage_fraction>=.8)&p.pct_expressing_measured.notna()
sc=ax.scatter(p.loc[good,'x'],p.loc[good,'cluster'],s=150*p.loc[good,'pct_expressing_measured'],c=p.loc[good,'z'],cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none')
ax.scatter(p.loc[~good,'x'],p.loc[~good,'cluster'],s=15,c='#999999',marker='x')
ax.set_xticks(range(len(show)),show,rotation=60,ha='right');ax.set_yticks(range(27),[str(k)+(' *' if k in held else '') for k in range(27)])
ax.invert_yaxis();ax.set_ylabel('Exploratory cluster; * = identity/QC review');ax.set_title('V7 Mast: full-gene lineage model and measured marker evidence',loc='left')
ax.set_xlim(-.7,len(show)-.3);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
cb=fig.colorbar(sc,ax=ax,pad=.02,shrink=.55);cb.set_label('Gene-scaled log1p(mean CP10K)')
for frac in [.1,.5,1]:ax.scatter([],[],s=150*frac,c='#777777',label=str(int(frac*100))+'% detected')
ax.legend(loc='upper left',bbox_to_anchor=(1.13,1),frameon=False)
fig.text(.08,.01,'Dot size: detection among source-measured cells. x: <80% source coverage.\nClusters are diagnostic partitions, not 27 accepted subtypes. Original V7 labels unchanged.',fontsize=9)
fig.subplots_adjust(bottom=.15,right=.78,left=.08,top=.94);fig.savefig(o/'Mast_marker_review.png',dpi=170);plt.close(fig)
manifest={'status':'completed_Mast_parent_and_QC_review_not_final_annotation','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':45036,'n_clusters_at_exploratory_res0_4':27,'n_parent_supported_candidate_cells':int(review.loc[~review.cluster.isin(held),'n_cells'].sum()),'n_competing_program_QC_review_cells':int(review.loc[review.cluster.isin(held),'n_cells'].sum()),'n_cycling_enriched_cluster_cells':int(n[16]),'n_Mast_sources_with_native_headers':len(gses),'original_V7_labels_preserved':True,'new_final_label_changes':0,'limits':['Manual review routes are not final doublet decisions or automatic exclusion.','Unintegrated source/feature effects are substantial; do not name 27 subtypes.','No outside-parent recruitment or final cleaned lineage refit yet.'],'outputs':[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir() if p.is_file()]}
(o/'complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest))
