"""Case-block inference for mean excess spatial marks; no cell-level P.

Whole-case sign flips preserve all within-case cellular dependencies. They
require independent case units and symmetry of case effects under the null;
calibration under fixed synthetic mechanisms is separate, not proof of either.
"""
import hashlib
import numpy as np
from CM02_marked_spatial_core_v1 import cohort_summary as descriptive_summary


def case_signflip_p(delta,case_keys,B=2000,seed=20261022):
    d=np.asarray(delta,float)
    if len(d)!=len(case_keys) or len(set(case_keys))!=len(case_keys):raise ValueError('case_block_keys_not_unique')
    if not np.isfinite(d).all() or not len(d):raise ValueError('invalid_case_effects')
    flips=[]
    for key in case_keys:
        task_seed=int.from_bytes(hashlib.sha256((str(seed)+':'+str(key)).encode()).digest()[:8],'little')
        flips.append(np.random.default_rng(task_seed).choice([-1.,1.],B))
    values=np.column_stack(flips)@d/len(d)
    return float((1+np.count_nonzero(np.abs(values)>=abs(d.mean())-1e-12))/(B+1))


def cohort_summary(observed,null,case_keys,minimum=20,inference=True,seed=20261022):
    result=descriptive_summary(observed,null,case_keys,minimum=minimum,inference=False,seed=seed)
    if not inference:return result
    if len(observed)<minimum:result['status']='COVERAGE_BELOW_PRIMARY_MINIMUM';return result
    delta=np.asarray(observed)-np.mean(np.asarray(null),axis=1)
    result.update(status='ESTIMABLE',patient_signflip_P=case_signflip_p(delta,case_keys,2000,seed),
                  inferential_unit='whole_author_patient_key_sign_flip')
    return result


def bh_fixed_family(results,family_size):
    slots=sorted([(k,v['patient_signflip_P']) for k,v in results.items() if 'patient_signflip_P' in v],key=lambda x:x[1])
    q=1.;harmonic=sum(1/k for k in range(1,family_size+1))
    for rank in range(len(slots),0,-1):
        key,p=slots[rank-1];q=min(q,p*family_size/rank)
        results[key]['BH_fixed_family']=float(min(q,1.))
        results[key]['BY_arbitrary_dependence_sensitivity']=float(min(q*harmonic,1.))
