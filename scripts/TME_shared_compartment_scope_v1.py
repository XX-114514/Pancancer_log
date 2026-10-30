"""Project comparison scope, not a new cell classifier or doublet detector."""
import numpy as np

def shared_compartment(marker,low_reference,high_reference,native_bits,family_bitmasks):
    marker=np.asarray(marker);low=np.asarray(low_reference);high=np.asarray(high_reference);bits=np.asarray(native_bits);masks=np.asarray(family_bitmasks)
    if marker.ndim!=1 or marker.shape!=low.shape or marker.shape!=high.shape or marker.shape!=bits.shape:raise ValueError('Comparison-scope row axes differ')
    if masks.ndim!=1 or not len(masks) or masks.dtype.kind!='u' or bits.dtype.kind!='u':raise ValueError('Expected frozen unsigned native family bitmasks')
    for a in [marker,low,high]:
        if a.dtype.kind not in 'iu' or (a<-1).any() or (a>=len(masks)).any():raise ValueError('Reference/marker compartment code invalid')
    result=marker.copy();reason=np.where(marker>=0,1,0).astype(np.int8)
    concordant=(low>=0)&(low==high);chosen_mask=masks[np.maximum(low,0)];supported=concordant&((bits&chosen_mask)!=0)
    result[supported]=low[supported]
    reason[supported&(marker==low)]=2
    reason[supported&(marker<0)]=3
    reason[supported&(marker>=0)&(marker!=low)]=4
    return result,reason

def counterexamples():
    # T/NK, B/plasma, myeloid, mast and epithelial compartments. B/plasma
    # and myeloid may share some coherent biology; no doublet is inferred.
    masks=np.asarray([1,2,4,8,16],np.uint32)
    marker=np.asarray([3,4,4,-1,-1,0,-1,2,-1],np.int16)
    low=np.asarray([3,0,0,1,0,0,0,2,-1],np.int16)
    high=np.asarray([3,0,0,1,2,0,0,2,0],np.int16)
    native=np.asarray([8,1,16,2,1,0,0,4|2,1],np.uint32)
    got,reason=shared_compartment(marker,low,high,native,masks)
    np.testing.assert_array_equal(got,[3,0,4,1,-1,0,-1,2,-1])
    np.testing.assert_array_equal(reason,[2,4,1,3,0,1,0,2,0])
    return dict(n_hand_checked_cases=9,reference_alone_cannot_override_marker=True,unmeasured_or_negative_native_support_not_invented=True,missing_native_does_not_erase_existing_marker_scope=True,shared_RNA_not_automatically_doublet=True)
