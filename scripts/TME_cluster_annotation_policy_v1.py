"""Prospective cluster-level RNA proposals; display marker is not the classifier.

Fractions summarize coherent native single-cell modules. Thresholds are frozen
engineering rules for this exploratory pass, not thresholds taken from a paper.
No old cell label, external prediction or clinical endpoint enters these calls.
"""
import numpy as np

STATE_PANELS={
 'naive_memory_related':['TCF7','LEF1','CCR7','IL7R'],
 'cytotoxic_related':['GZMB','PRF1','NKG7','CTSW'],
 'dysfunction_related':['PDCD1','TOX','LAG3','HAVCR2','CXCL13'],
 'regulatory_related':['FOXP3','IL2RA','CTLA4','IKZF2'],
 'interferon_related':['IFIT1','IFIT3','ISG15','MX1'],
 'cycling_related':['MKI67','TOP2A','UBE2C','STMN1'],
 'CXCL13_helper_related':['CXCL13','ICOS','MAF','PDCD1'],
 'tissue_residency_related':['ITGAE','ZNF683','CXCR6','CD69']}

GATES={
 'T':[['CD3D','CD3E'],['TRAC','TRBC1','TRBC2','TRDC'],['CD247','LCK','LAT','ZAP70']],
 'NK':[['NKG7'],['KLRD1','NCR1','KLRF1'],['GNLY','PRF1']],
 'B':[['CD79A','CD79B'],['MS4A1','CD19','CD22']],
 'Plasma':[['MZB1','JCHAIN'],['SDC1','DERL3','TNFRSF17']],
 'Myeloid':[['LST1','LYZ'],['AIF1','TYROBP','CSF1R']],
 'Endothelial':[['PECAM1','CDH5'],['CLDN5','VWF','ESAM']],
 'Fibroblast':[['COL1A1','COL1A2','COL3A1'],['DCN','LUM']],
 'Epithelial':[['EPCAM','CDH1','KRT19'],['KRT8','KRT18','KRT5','KRT14']],
 'Mural':[['RGS5','CSPG4','PDGFRB'],['ACTA2','TAGLN','MYH11']],
 'Mast':[['TPSAB1','TPSB2','KIT'],['CPA3','HDC','MS4A2']],
 'Treg':[['FOXP3'],['IL2RA','CTLA4']],
 'CD4_T':[['CD4'],['IL7R','LTB','CCR7','FOXP3']],
 'CD8_T':[['CD8A'],['CD8B']],
 'gamma_delta_T_RNA_candidate':[['TRDC'],['TRGC1','TRGC2']]}
GATE_ORDER=list(GATES)

def coherent_flags(raw,genes):
 """Each OR module must be measurable; unmeasured modules stay unassessed."""
 by={g:j for j,g in enumerate(genes)};flags={};assessed={}
 for name,modules in GATES.items():
  columns=[[by[g] for g in m if g in by] for m in modules]
  assessed[name]=all(columns);v=np.ones(raw.shape[0],bool)
  if assessed[name]:
   for cols in columns:v&=np.asarray((raw[:,cols]>0).sum(axis=1)).ravel()>0
  else:v[:]=False
  flags[name]=v
 for name in ['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate']:
  assessed[name]&=assessed['T'];flags[name]&=flags['T']
 # Treg is nested in the CD4 lineage; it does not constitute a competing call.
 flags['CD4_T']&=~flags['Treg']
 return np.column_stack([flags[k] for k in GATE_ORDER]),assessed

def annotate_cluster(n,counts,assessed,minimum_group=30,family_fraction=.2):
 if n<=0:raise ValueError('Empty clusters cannot receive a proposal')
 fractions={k:float(counts.get(k,0)/n) if assessed.get(k,False) else None for k in GATE_ORDER}
 active=[k for k in ['T','NK','B','Plasma','Myeloid','Endothelial','Fibroblast','Epithelial','Mural','Mast'] if fractions[k] is not None and fractions[k]>=family_fraction]
 result=dict(major='unresolved',celltype='unresolved',status='insufficient_coherent_native_cluster_evidence',native_gate_fractions=fractions,coherent_supported_families=active,annotation_scope='provisional_cluster_level_RNA_not_independent_cell_truth',fine_candidates=[])
 if n<minimum_group:result['status']='small_cluster_descriptive_only';return result
 competing=[k for k in active if k not in ['T','NK']]
 if 'T' in active:
  if competing:result['status']='coherent_competing_lineages_review_required';return result
  result.update(major='T',celltype='Tcell',status='coherent_cluster_T_RNA_support')
  choices=[k for k in ['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate'] if fractions[k] is not None and fractions[k]>=.1 and fractions[k]/max(fractions['T'],1e-12)>=.5]
  result['fine_candidates']=choices
  if len(choices)==1:result.update(celltype=choices[0],status='coherent_cluster_fine_T_RNA_proposal')
  elif len(choices)>1:result['status']='fine_T_modules_mixed_keep_coarse_T'
  return result
 if 'NK' in active and not competing:
  if fractions['T'] is None:result['status']='T_module_unmeasured_NK_review_required';return result
  if fractions['T']>=.1:result.update(major='T_NK',celltype='T_NK',status='mixed_T_NK_context_keep_broad');return result
  result.update(major='NK',celltype='NKcell',status='coherent_cluster_NK_RNA_support');return result
 if len(competing)==1:result.update(major=competing[0],celltype=competing[0],status='non_T_lineage_in_common_pool_cluster_proposal')
 elif len(competing)>1:result['status']='coherent_competing_lineages_review_required'
 return result

def state_context(genes,own_log_mean,other_log_mean,own_detection,minimum_log_difference=.25):
 by={g:j for j,g in enumerate(genes)};out={}
 for name,panel in STATE_PANELS.items():
  use=[g for g in panel if g in by];cols=[by[g] for g in use];enough=len(cols)>=max(2,int(np.ceil(len(panel)*.5)))
  mean=float(np.mean(own_log_mean[cols])) if cols else None
  contrast=float(np.mean(own_log_mean[cols]-other_log_mean[cols])) if cols and other_log_mean is not None else None
  supportive=sum(own_detection[j]>=.1 and (other_log_mean is None or own_log_mean[j]-other_log_mean[j]>0) for j in cols)
  out[name]=dict(n_measured_genes=len(cols),n_panel_genes=len(panel),measured_genes=use,mean_native_logCP10K=mean,mean_log_difference_from_comparator=contrast,n_supportive_genes=int(supportive),descriptively_enriched=bool(enough and contrast is not None and contrast>=minimum_log_difference and supportive>=2),functional_or_clinical_identity_assigned=False)
 return out
