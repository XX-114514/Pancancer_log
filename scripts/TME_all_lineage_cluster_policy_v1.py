"""Frozen broad RNA proposals for all source-qualified TME root/child clusters.

Extend the accepted T branch to nonimmune tissue; never infer malignancy from
lineage. Marker ranking and descriptive programs remain separate from calls.
These thresholds are exploratory project rules, not published accuracy claims.
"""
import numpy as np
from TME_cluster_annotation_policy_v1 import GATES as T_BRANCH_GATES,state_context
from TME_assay_aware_T_anchor_policy_v1 import assay_aware_flags
from native_fine_resolution_v1 import FINE_ANCHOR_PATTERNS

GATES={k:[list(m) for m in v] for k,v in T_BRANCH_GATES.items()}
GATES['Myeloid']=[['LST1','LYZ','CLEC9A','CLEC4C','CSF3R'],['AIF1','TYROBP','CSF1R','XCR1','IL3RA','FCGR3B']]
GATES['Endothelial']=[['PECAM1','CDH5','PROX1','FLT4'],['CLDN5','VWF','ESAM','PDPN','CCL21']]
GATES.update({
 'Glial':[['OLIG1','OLIG2','SOX10'],['MBP','PLP1','MOG','MOBP','MAG']],
 'Astrocyte':[['GFAP','AQP4'],['SLC1A3','SLC1A2','SOX9']],
 'Neuron':[['SNAP25','SYT1'],['RBFOX3','STMN2']],
 'Erythroid':[['HBA1','HBA2','HBB'],['AHSP','ALAS2','GYPA','SLC4A1']],
 'Platelet':[['PF4','PPBP'],['GP9','GP1BA','ITGA2B']]})
FINE_BY_PARENT={'T':['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate'],
 'Myeloid':['Macrophage','Monocyte','cDC1','cDC2_related','pDC','Neutrophil'],
 'Endothelial':['Lymphatic_endothelial_related']}
for name in FINE_BY_PARENT['Myeloid']+FINE_BY_PARENT['Endothelial']:GATES[name]=[list(m) for m in FINE_ANCHOR_PATTERNS[name]]
GATE_ORDER=list(GATES)
BROAD=['T','NK','B','Plasma','Myeloid','Endothelial','Fibroblast','Epithelial','Mural','Mast','Glial','Astrocyte','Neuron','Erythroid','Platelet']

def coherent_flags(raw,genes):
 by={g:j for j,g in enumerate(genes)};flags={};assessed={}
 for name,modules in GATES.items():
  cols=[[by[g] for g in module if g in by] for module in modules];assessed[name]=all(cols);v=np.ones(raw.shape[0],bool)
  for c in cols:v&=np.asarray((raw[:,c]>0).sum(axis=1)).ravel()>0 if c else False
  flags[name]=v
 base,a,mode=assay_aware_flags(raw,genes);flags['T']=base[:,0];assessed['T']=a['T']
 for parent,names in FINE_BY_PARENT.items():
  for name in names:flags[name]&=flags[parent];assessed[name]&=assessed[parent]
 flags['CD4_T']&=~flags['Treg']
 return np.column_stack([flags[k] for k in GATE_ORDER]),assessed

def annotate_cluster(n,counts,assessed,minimum_group=30,family_fraction=.2):
 if n<=0:raise ValueError('Empty clusters cannot receive a proposal')
 fractions={k:float(counts.get(k,0)/n) if assessed.get(k,False) else None for k in GATE_ORDER};active=[k for k in BROAD if fractions[k] is not None and fractions[k]>=family_fraction]
 result=dict(major='unresolved',celltype='unresolved',status='insufficient_coherent_native_cluster_evidence',native_gate_fractions=fractions,coherent_supported_families=active,annotation_scope='provisional_all_lineage_cluster_RNA_not_cell_truth_or_malignancy',fine_candidates=[])
 if n<minimum_group:result['status']='small_cluster_descriptive_only';return result
 other=[k for k in active if k not in ['T','NK']]
 if 'T' in active:
  if other:result['status']='coherent_competing_lineages_review_required';return result
  parent='T';result.update(major='T',celltype='Tcell',status='coherent_cluster_T_RNA_support')
 elif 'NK' in active and not other:
  if fractions['T'] is None:result['status']='T_module_unmeasured_NK_review_required';return result
  if fractions['T']>=.1:result.update(major='T_NK',celltype='T_NK',status='mixed_T_NK_context_keep_broad');return result
  result.update(major='NK',celltype='NKcell',status='coherent_cluster_NK_RNA_support');return result
 elif len(other)==1:
  parent=other[0];result.update(major=parent,celltype=parent,status='coherent_broad_tissue_lineage_RNA_proposal')
 else:
  if len(other)>1:result['status']='coherent_competing_lineages_review_required'
  return result
 choices=[k for k in FINE_BY_PARENT.get(parent,[]) if fractions[k] is not None and fractions[k]>=.1 and fractions[k]/max(fractions[parent],1e-12)>=.5];result['fine_candidates']=choices
 if len(choices)==1:result.update(celltype=choices[0],status='coherent_cluster_fine_lineage_RNA_proposal')
 elif len(choices)>1:result['status']='fine_modules_mixed_keep_coarse_lineage'
 return result
