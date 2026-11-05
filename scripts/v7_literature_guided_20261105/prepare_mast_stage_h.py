import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_h_mast_source_robust_refit_v1';o.mkdir(exist_ok=False)
g=pd.read_pickle(r/'stage_g_mast_marker_QC_v2/mast_cell_QC_marker_review_PRIVATE.pkl.gz');a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');assert g.index.equals(a.obs_names)
required=['TPSAB1','TPSB2','CPA3','KIT','MS4A2'];v=a.X[:,a.var_names.get_indexer(required)].toarray()
# Cross-functional Mast identity support: a tryptase and a second Mast component.
tryptase=(v[:,:2]>=2).any(axis=1);support=(v[:,2:]>=2).any(axis=1);parent=tryptase&support
competitors=['T_competitor','myeloid_competitor','stromal_competitor','epithelial_competitor'];anycomp=np.zeros(len(g),bool)
for c in competitors:anycomp|=g[c+'_n_genes_ge2UMI'].to_numpy()>=2
q=g.doublet_calibration_evaluable.to_numpy()&(g.doublet_within_sample_method_percentile.to_numpy()>=.95)
joint=parent&anycomp&q
weakcomp=~parent&anycomp
routes=np.select([weakcomp,~parent,joint,parent&anycomp],['parent_routing_review_competitor_with_weak_Mast','Mast_parent_low_marker_review','joint_doublet_candidate_review','Mast_cross_program_review'],default='Mast_coarse_identity_supported')
g['Mast_tryptase_plus_support_raw2']=parent;g['Mast_competing_program_raw2']=anycomp;g['review_route_v3']=routes
# Fit only one conservative sensitivity pool. All 45036 cells remain in the review overlay.
fit=parent&~anycomp&g.panMC_published_QC_pass.to_numpy()
g['Mast_source_robust_refit_candidate_v2']=fit
g.to_pickle(o/'mast_cell_review_v3_PRIVATE.pkl.gz');g.to_csv(o/'mast_allcell_metadata_PRIVATE.tsv.gz',sep='\t',index_label='cell_id')
counts=g.review_route_v3.value_counts().rename_axis('route_v3').reset_index(name='n_cells');counts.to_csv(o/'review_routes_v3.tsv',sep='\t',index=False)
summary=[]
for ci,t in g.groupby('mast_review_cluster_res0.4',observed=True):
 row={'cluster':int(ci),'n_cells':len(t),'n_Mast_cross_functional_support':int(t.Mast_tryptase_plus_support_raw2.sum()),'n_fit_pool':int(t.Mast_source_robust_refit_candidate_v2.sum())}
 for label in counts.route_v3:row[label+'_n']=int(t.review_route_v3.eq(label).sum())
 summary.append(row)
pd.DataFrame(summary).to_csv(o/'mast27_refined_parent_QC_routes.tsv',sep='\t',index=False)
assert counts.n_cells.sum()==45036 and fit.sum()>10000
manifest={'status':'H_input_prepared','UTC':datetime.now(timezone.utc).isoformat(),'n_all_retained':len(g),'n_refit_pool':int(fit.sum()),'n_parent_cross_functional_support':int(parent.sum()),'n_joint_doublet_candidates':int(joint.sum()),'n_parent_weak_competitor_review':int(weakcomp.sum()),'route_counts':dict(zip(counts.route_v3,counts.n_cells.astype(int))),'parent_rule':'At least one of TPSAB1/TPSB2 and at least one of CPA3/KIT/MS4A2, each >=2 raw UMI. Marker support, not calibrated identity probability.','refit_pool_rule':'Cross-functional parent support; no >=2-gene competing parent at >=2 UMI; published panMC numeric QC with measurable mitochondrial genes. This is a conservative sensitivity fit, not a final cell exclusion.','model_feature_rule':'Measured in >=80% of the 49 actual Mast sources; HB genes omitted from HVG candidates. Full RNA assay and marker universe retained.','new_final_annotation_changes':0}
(o/'input_complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest))
