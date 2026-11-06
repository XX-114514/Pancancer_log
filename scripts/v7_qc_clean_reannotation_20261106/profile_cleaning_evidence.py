import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
r=PROJECT_ROOT;old=r/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';run=r/'runs/20261106_v7_QC_clean_reannotation_v1';run.mkdir(exist_ok=False);o=run/'stage_a_evidence';o.mkdir()
config={'run':'20261106_v7_QC_clean_reannotation_v1','UTC_start':datetime.now(timezone.utc).isoformat(),'authorization':'User explicitly authorized cleaning and continuing lineage review on 2026-11-06.','scope':'Frozen V7; first actual cleanup/refit milestone covers Tcell_6 and Mast; global previously-unscored units audited separately.','preservation':'Original data and labels remain intact; exclusions apply to a new versioned analysis view.','interpretation':'Operational QC exclusions can be probable doublet or mixed-RNA; no ground-truth confirmation claimed.','next_decision':'Review joint candidates and independent lineage anchors before freezing primary/exploratory masks.'};(run/'plan.json').write_text(json.dumps(config,indent=2))
panels={'T':[['CD3D','CD3E','CD3G'],['TRAC','TRBC1','TRBC2']],'B':[['MS4A1','CD19','BANK1'],['CD79A','CD79B']],'Plasma':[['MZB1','DERL3'],['JCHAIN','SDC1']],'Myeloid':[['CSF1R','CD14','CD163','C1QA'],['LST1','MS4A7']],'Stromal':[['COL1A1','COL1A2'],['DCN','LUM']],'Epithelial':[['EPCAM'],['KRT8','KRT18','KRT19','KRT14','KRT17']],'Endothelial':[['VWF','CLDN5'],['PECAM1','KDR']],'Platelet':[['PF4','PPBP'],['GP9','ITGA2B']],'Mast':[['TPSAB1','TPSB2'],['CPA3','KIT','MS4A2']]}
datasets=[('T6',old/'stage_j_T6_fullgene/T6_v7_fullgene_raw.h5ad',old/'stage_n_T6_measurement_scope_v1/T6_hierarchical_measured_review_v2_PRIVATE.pkl.gz','T6_review_cluster_source_res0.4'),('Mast',old/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad',old/'stage_i_mast_extended_review_v1/Mast_annotation_review_v1_PRIVATE.pkl.gz',None)]
rows=[];cohorts=[];programs=[]
for name,rawpath,metapath,clustercol in datasets:
 a=ad.read_h5ad(rawpath);q=pd.read_pickle(metapath);assert a.obs_names.equals(q.index)
 if name=='Mast':
  md=pd.read_csv(old/'stage_c2_mast_initial_diagnostic/mast_initial_metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert md.index.equals(q.index);cl=md['RNA_snn_res.0.4'].astype(int).to_numpy();joint=q.review_route_v4.eq('joint_doublet_candidate_review').to_numpy()
 else:cl=q[clustercol].astype(int).to_numpy();joint=q.joint_competitor_doublet_candidate.to_numpy()
 genes=list(dict.fromkeys(g for groups in panels.values() for group in groups for g in group));ids=a.var_names.get_indexer(genes);assert (ids>=0).all();raw=a.X[:,ids].toarray();ix={g:i for i,g in enumerate(genes)};total=np.asarray(a.X.sum(axis=1)).ravel();strong=np.zeros(len(q),bool);compgeneUMI=np.zeros(len(q));anycomp=np.zeros(len(q),bool)
 for other,groups in panels.items():
  if (name=='T6' and other=='T') or (name=='Mast' and other=='Mast'):continue
  groupsraw=[raw[:,[ix[g] for g in group]] for group in groups];anchor=(groupsraw[0]>=1).any(axis=1)&(groupsraw[1]>=1).any(axis=1);gg=list(dict.fromkeys(g for group in groups for g in group));x=raw[:,[ix[g] for g in gg]];program=((x>=1).sum(axis=1)>=3)|((x>=2).sum(axis=1)>=2);strong|=anchor&program;anycomp|=program;compgeneUMI+=x.sum(axis=1);q['adjudication_'+other+'_anchor']=anchor&program
  programs.append({'lineage':name,'competing_program':other,'n_joint_candidates':int((joint&program).sum()),'n_joint_with_distinct_anchor':int((joint&anchor&program).sum()),'n_not_joint_with_distinct_anchor':int((~joint&anchor&program).sum())})
 q['review_cluster_before_cleanup']=cl;q['joint_review_original']=joint;q['distinct_competing_anchor_supported']=strong;q['competing_panel_RNA_fraction']=compgeneUMI/np.maximum(total,1);q['raw_full_library_total']=total;q.to_pickle(o/(name+'_adjudication_evidence_PRIVATE.pkl.gz'))
 for c in np.unique(cl):
  take=cl==c;z=q.iloc[np.flatnonzero(take)];jj=joint[take];rows.append({'lineage':name,'old_review_cluster':int(c),'n_cells':int(take.sum()),'n_joint_candidates':int(jj.sum()),'n_joint_distinct_anchor':int((joint&strong&take).sum()),'n_not_joint_distinct_anchor':int((~joint&strong&take).sum()),'doublet_percentile_median':float(pd.to_numeric(z.doublet_within_sample_method_percentile,errors='coerce').median()),'decontx_median':float(pd.to_numeric(z.decontx_contamination,errors='coerce').median()),'competing_RNA_fraction_median':float(z.competing_panel_RNA_fraction.median()),'n_sources':z.gse_id.nunique()})
 cohorts.append({'lineage':name,'n_cells':len(q),'n_joint_candidates':int(joint.sum()),'n_joint_distinct_anchor':int((joint&strong).sum()),'n_joint_without_distinct_anchor':int((joint&~strong).sum())})
pd.DataFrame(rows).to_csv(o/'cluster_adjudication_evidence.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(programs).to_csv(o/'competing_anchor_review.tsv',sep='\t',index=False);out={'status':'evidence_profile_completed_before_mask_freeze','UTC':datetime.now(timezone.utc).isoformat(),'cohorts':cohorts,'n_exclusions_applied_yet':0};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps({'summary':out,'cluster_rows':rows,'program_rows':programs}))
