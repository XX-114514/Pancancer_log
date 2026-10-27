"""Reuse the frozen direct fine-cue rule after qualified T coarse annotation.

No new subtype thresholds, clinical selection or group-major override. The
gamma label remains an RNA candidate. Existing fine identities stay unchanged.
"""
import numpy as np

TYPES=['T','Treg','CD4_T','CD8_T','T_gamma_delta_RNA_candidate']
FROZEN_BITS=['Treg','CD4_T','CD8_T','T_gamma_delta_RNA_candidate']
REASONS=['outside_coarse_T_target','direct_native_and_unique_reference_fine_RNA',
         'gamma_RNA_candidate_requires_orthogonal_review','no_direct_fine_intersection',
         'multiple_direct_fine_candidates','coarse_T_native_support_not_qualified','zero_native_RNA']

def direct_fine(current_major,current_fine,own,reference,native_family,coarse_review,positive):
    major=np.asarray(current_major);fine=np.asarray(current_fine)
    own=np.asarray(own);reference=np.asarray(reference);native=np.asarray(native_family);review=np.asarray(coarse_review);positive=np.asarray(positive,bool)
    if major.ndim!=1 or any(v.shape!=major.shape for v in [fine,own,reference,native,review,positive]):raise ValueError('Cell axes differ')
    for v in [own,reference,native,review]:
        if v.dtype.kind not in 'iu' or (v<0).any():raise ValueError('Nonnegative categorical bit fields required')
    target=(major=='T')&(fine=='T');supported=(native==1)|(review==6)
    both=(own&reference&15).astype(np.uint8)
    # Treg is a child of CD4_T in the already frozen hierarchy.
    both=np.where((both&1)!=0,both&13,both).astype(np.uint8)
    n=np.asarray([bin(k).count('1') for k in range(16)],np.uint8)[both]
    reason=np.zeros(len(major),np.uint8);reason[target]=3
    reason[target&(n>1)]=4;reason[target&~supported]=5;reason[target&~positive]=6
    eligible=target&supported&positive&(n==1);proposed=np.zeros(len(major),np.uint8)
    for j,name in enumerate(FROZEN_BITS):proposed[eligible&(both==(1<<j))]=TYPES.index(name)
    reason[eligible]=1;reason[eligible&(proposed==4)]=2
    return {'proposed_fine_code':proposed,'fine_review_code':reason,'direct_intersection_bitmask':both,'eligible_coarse_T':target}
