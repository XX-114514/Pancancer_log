"""Evidence-role proposal for an already frozen pDC disagreement cohort.

Shared APC, cytotoxic and secretory genes cannot by themselves establish a
competing lineage. Specific multi-gene competing programs remain review cases.
This is RNA-context adjudication, not purity, ontogeny or calibrated confidence.
"""
import numpy as np
CORE_PANELS={
 'B_specific':(['MS4A1','CD19','CD79B','PAX5','BANK1','CD22'],2),
 'T_receptor':(['CD3D','CD3E','CD247','TRAC','TRBC1','TRBC2'],2),
 'NK_specific':(['GNLY','KLRD1','NCR1','KLRF1'],2),
 'Endothelial':(['PECAM1','CDH5','CLDN5','VWF','EMCN'],2),
 'Fibroblast':(['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1'],3),
 'Glial':(['GFAP','AQP4','OLIG2','SOX10','PLP1'],2),
}
STATES=['context_supported_pDC_related_proposal','competing_specific_context_review',
        'unmeasured_required_context_review','additional_pDC_context_not_met_review']

def context_proposal(values,genes,measured,minimum_native_value=2):
 x=np.asarray(values);m=np.asarray(measured)
 if (x.ndim!=2 or x.shape[1]!=len(genes) or len(set(genes))!=len(genes)
     or m.dtype!=bool or m.shape!=(len(genes),) or minimum_native_value!=2
     or not np.isfinite(x[:,m]).all() or (x[:,m]<0).any()):
  raise ValueError('Unique native gene values/measurement axes and frozen threshold required')
 n=len(x);lookup={g:i for i,g in enumerate(genes)}
 def atleast(names,k):
  available=[lookup[g] for g in names if g in lookup and m[lookup[g]]]
  h=(x[:,available]>=minimum_native_value).sum(axis=1) if available else np.zeros(n,int)
  unmeasured=len(names)-len(available)
  return np.where(h>=k,1,np.where(h+unmeasured<k,0,2)).astype(np.uint8)
 def both(*states):
  a=np.stack(states);return np.where((a==0).any(axis=0),0,np.where((a==2).any(axis=0),2,1)).astype(np.uint8)
 specific=both(atleast(['LILRA4'],1),atleast(['PLD4','TSPAN13'],1))
 competing={name:atleast(names,k) for name,(names,k) in CORE_PANELS.items()}
 competing['Macrophage_specific']=both(atleast(['C1QA'],1),atleast(['C1QB'],1),atleast(['LYZ','LST1'],1))
 competing['Plasma_specific']=both(atleast(['SDC1','TNFRSF17'],1),atleast(['MZB1','JCHAIN','DERL3','XBP1'],2))
 all_core=np.stack(list(competing.values()));any_true=(all_core==1).any(axis=0);any_unknown=(all_core==2).any(axis=0)
 code=np.full(n,3,np.uint8);code[(specific==2)|any_unknown]=2;code[any_true]=1;code[(specific==1)&~any_true&~any_unknown]=0
 return {'disposition_code':code,'additional_pDC_specific_state':specific,
         'competing_core_state':np.stack(list(competing.values()),axis=1),
         'competing_core_order':list(competing),'disposition_order':list(STATES)}
