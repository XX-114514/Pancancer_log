"""Alternative native T evidence only when all TCR-chain alternatives are unmeasured.

This does not infer TCR chains or rescue measured-but-negative chains. Both
CD3D/CD3E, an accessory anchor and two signaling genes must co-occur in each
supporting cell. Thresholds/classifier stay frozen; evidence mode is explicit.
"""
import numpy as np
from TME_cluster_annotation_policy_v1 import GATES,GATE_ORDER,coherent_flags

CHAIN=['TRAC','TRBC1','TRBC2','TRDC']
ACCESSORY=['CD2','CD7','TRAT1']
SIGNAL=['CD247','LCK','LAT','ZAP70','BCL11B']

def assay_aware_flags(raw,genes):
 flags,assessed=coherent_flags(raw,genes);by={g:j for j,g in enumerate(genes)}
 mode='original_chain_and_CD3_signaling_native_module'
 if assessed['T']:return flags,assessed,mode
 if any(g in by for g in CHAIN):return flags,assessed,'incomplete_other_T_module_not_rescued'
 possible=all(g in by for g in ['CD3D','CD3E']) and any(g in by for g in ACCESSORY) and sum(g in by for g in SIGNAL)>=2
 if not possible:return flags,assessed,'chain_unmeasured_and_alternative_anchor_unassessable'
 def positive(g):return np.asarray((raw[:,[by[g]]]>0).sum(axis=1)).ravel()>0
 accessory=np.zeros(raw.shape[0],bool);signals=np.zeros(raw.shape[0],np.int16)
 for gene in ACCESSORY:
  if gene in by:accessory|=positive(gene)
 for gene in SIGNAL:
  if gene in by:signals+=positive(gene)
 alt=positive('CD3D')&positive('CD3E')&accessory&(signals>=2);assessed['T']=True;flags[:,GATE_ORDER.index('T')]=alt
 for name in ['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate']:
  modules=GATES[name];valid=all(any(g in by for g in module) for module in modules);assessed[name]=valid;value=alt.copy()
  for module in modules:
   module_positive=np.zeros(raw.shape[0],bool)
   for gene in module:
    if gene in by:module_positive|=positive(gene)
   value&=module_positive
  flags[:,GATE_ORDER.index(name)]=value
 flags[:,GATE_ORDER.index('CD4_T')]&=~flags[:,GATE_ORDER.index('Treg')]
 return flags,assessed,'T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence'
