import pathlib,sys,json,hashlib
import pandas as pd,numpy as np
h,r=map(pathlib.Path,sys.argv[1:3]);assert (r/'consensus_v1/DONE.json').exists()
for pool in ['Microglia','Endothelial','Myeloid']:assert (r/f'transfer_v1/{pool}.DONE.json').exists(),pool
out=r/'working_annotation_v1';out.mkdir(exist_ok=False)
L=pd.read_csv(h/'p0/cell_labels_r3.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('cell_id');assert L.index.is_unique
D=L.copy();D['working_subtype']=D.r2_label.where(D.r2_label!='',D.label_major);D['lineage']=D.major_label.str.replace('^r3_','',regex=True)
D['latest_action']=D.r2_action;D['latest_round']='r2';D['annotation_status']='outside_refined_pools';D['state']='';D['consensus_label']='';D['consensus_stability']=np.nan;D['include_main_subtype']=False;D['include_lineage_total']=D.major_label.str.startswith('r3_');D['transfer_conf']=np.nan;D['transfer_label_raw']='';D['rare_candidate']='';D['stability_scope']='not_tested'
for pool in ['Microglia','Alveolar_macrophage','Kupffer','Sinusoidal_EC']:
 ix=D.index[D.major_label=='r3_'+pool];D.loc[ix,'annotation_status']='observed_tissue_no_seed_consensus'
 low=ix[(D.loc[ix,'r1_action']=='keep_flag_low_depth')|(D.loc[ix,'r2_action']=='keep_flag_low_depth')];D.loc[low,'annotation_status']='unresolved_low_depth'
 if pool=='Microglia':
  D.loc[ix,'state']=D.loc[ix,'r2_label'].str.replace('^Microglia_?','',regex=True);D.loc[ix,'working_subtype']='Microglia'
R=pd.read_csv(r/'review_v1/lineage_r3_reviewed.tsv',sep='\t',dtype=str,keep_default_na=False)
for pool in ['Myeloid','Endothelial']:
 o=pd.read_csv(h/f'p1/obs/r3_{pool}.obs.tsv.gz',sep='\t',usecols=['obs_name','leiden_r1'],dtype=str).set_index('obs_name');rv=R[R.pool==pool].set_index('cluster')
 D.loc[o.index,'latest_round']='r3';D.loc[o.index,'latest_action']=o.leiden_r1.map(rv.action).values;D.loc[o.index,'state']=o.leiden_r1.map(rv.state).values;D.loc[o.index,'working_subtype']=o.leiden_r1.map(rv.r3_label).values
for pool in ['T_NK','B_Plasma','Stromal','Myeloid','Endothelial']:
 c=pd.read_csv(r/f'consensus_v1/{pool}.cells.tsv.gz',sep='\t').set_index('cell_id');c=c[c.retained_pool]
 ix=c.index;D.loc[ix,'consensus_label']=c.consensus_label;D.loc[ix,'consensus_stability']=c.stability;D.loc[ix,'stability_scope']='resolution1_identity_consensus';D.loc[ix,'annotation_status']='unstable_or_nonidentity_component'
 D.loc[ix,'latest_action']=c.review_action;D.loc[ix,'working_subtype']=c.review_label
 D.loc[ix,'include_main_subtype']=c.include_main_subtype
 good=c.index[c.include_main_subtype];D.loc[good,'annotation_status']='observed_stable';D.loc[good,'working_subtype']=c.loc[good,'consensus_label']
 for action,status in [('keep_state','state'),('keep_flag_source_restricted','source_restricted'),('keep_flag_sample_specific','sample_specific')]:
  ii=c.index[c.review_action==action];D.loc[ii,'annotation_status']=status
 ii=c.index[c.low_depth];D.loc[ii,'annotation_status']='unresolved_low_depth'
# R3 non-keep cells are preserved as explicit pending records, never silently deleted.
pending=D.include_lineage_total & D.annotation_status.eq('outside_refined_pools');D.loc[pending,'annotation_status']='reassigned_lineage_subtype_pending'
held=D.latest_action.eq('set_aside_review')|D.major_label.str.contains('set_aside_review');removed=D.latest_action.str.startswith('remove')|D.major_label.str.contains('remove_')|D.major_label.isin(['Doublet_like','Low_quality']);D.loc[held,'annotation_status']='set_aside_review';D.loc[removed,'annotation_status']='excluded_qc';D.loc[held|removed,['include_main_subtype','include_lineage_total']]=False
for pool in ['T_NK','B_Plasma','Stromal','Microglia','Myeloid','Endothelial']:
 base=h/'p0/transfer' if pool in ['T_NK','B_Plasma','Stromal'] else r/'transfer_v1'
 t=pd.read_csv(base/f'{pool}.cells.tsv.gz',sep='\t').set_index('cell_id');assert t.index.is_unique
 assert (D.loc[t.index,'major_label']=='r3_'+pool).all()
 assert D.loc[t.index,'include_lineage_total'].all()
 D.loc[t.index,'transfer_label_raw']=t.transfer_label;D.loc[t.index,'working_subtype']=t.transfer_label;D.loc[t.index,'transfer_conf']=t.transfer_conf;D.loc[t.index,'include_main_subtype']=False
 D.loc[t.index,'annotation_status']=np.where(t.transfer_status=='transferred','transferred','unresolved_low_depth')
# Requested high-resolution candidates are recorded only for stable observed cells.
# Their 10-seed stability was tested at resolution 1, not at the finer split:
# retain parent label for main proportions until fine-level validation.
o=pd.read_csv(h/'p1/obs/r2_T_NK.obs.tsv.gz',sep='\t',usecols=['obs_name','leiden_r3'],dtype=str).set_index('obs_name')
rare={'16':'CD8_TEMRA_like','36':'NK_cytotoxic','4':'MAIT_Th17_like','38':'ILC3_like','37':'CD8_naive_like'}
for cluster,label in rare.items():
 ix=o.index[o.leiden_r3==cluster];ix=ix[D.loc[ix,'annotation_status'].eq('observed_stable')];D.loc[ix,'rare_candidate']=label
# Cluster42 has TRDC/proliferation evidence insufficient to force NK identity: not overwritten.
assert len(D)==3347659 and D.index.equals(L.index)
assert not D.loc[D.annotation_status.isin(['transferred','unresolved_low_depth','set_aside_review','excluded_qc','state','sample_specific']),'include_main_subtype'].any()
D.to_csv(out/'cell_labels_working_v1.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
D.groupby(['lineage','annotation_status','include_main_subtype'],dropna=False).size().rename('cells').reset_index().to_csv(out/'status_summary.tsv',sep='\t',index=False)
D.groupby(['lineage','working_subtype','annotation_status'],dropna=False).size().rename('cells').reset_index().to_csv(out/'subtype_summary.tsv',sep='\t',index=False)
D[D.rare_candidate!=''].groupby('rare_candidate').size().rename('cells').to_csv(out/'rare_candidates.tsv',sep='\t')
summary={'rows':len(D),'unique_cell_id':D.index.nunique(),'main_subtype_cells':int(D.include_main_subtype.sum()),'lineage_total_cells':int(D.include_lineage_total.sum()),'status_counts':D.annotation_status.value_counts().to_dict(),'scientific_final':False,'limits':['41 R3 clusters have evidence-based working review, not independent biological validation','rare fine splits are candidates only; stability assessed at resolution1','four tissue pools lack 10-seed consensus','E1 October comparison and downstream sample thresholds not yet accepted','Microglia Mac_glioma transfers retain original pool and are sensitivity-only']}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
