"""Exact Boolean signal overlap, without diagnosis or inferred identity."""
import numpy as np


def contingency(first,second):
    a,b=np.asarray(first),np.asarray(second)
    if a.ndim!=1 or b.shape!=a.shape or a.dtype!=bool or b.dtype!=bool:raise ValueError('Expected matched Boolean lineage support vectors')
    both=int(np.count_nonzero(a&b));ao=int(np.count_nonzero(a&~b));bo=int(np.count_nonzero(~a&b));neither=int(np.count_nonzero(~a&~b))
    n=len(a);na=ao+both;nb=bo+both
    return {'n_cells':n,'n_first_only':ao,'n_second_only':bo,'n_both':both,'n_neither':neither,
            'n_first_supported':na,'n_second_supported':nb,'fraction_first_with_second':both/na if na else None,
            'fraction_second_with_first':both/nb if nb else None,
            'fraction_both_in_group':both/n if n else None,'signal_overlap_not_doublet_or_identity_diagnosis':True}
