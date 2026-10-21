"""Resolve the original unweighted network members into frozen authority IDs.

No learned NMF weights or individual records are exported. Descriptive
parent-state rates are aggregated by reviewed case before source summaries.
"""
from collections import defaultdict
from itertools import product
from pathlib import Path
import hashlib
import os

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def resolve(atlas,base,anchors):
    path=Path(os.environ['PANCANCER_NETWORK_MEMBER_FILE'])
    assert sha(path)==os.environ['PANCANCER_NETWORK_MEMBER_SHA256']
    table=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False,
                      usecols=['cm','subCluster'])
    names=sorted(table.loc[table.cm.eq('CM02'),'subCluster'])
    assert len(names)>2 and len(names)==len(set(names))
    z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    nodes=np.asarray(read_elem(z['obs/node_id'])).astype(str)
    major=np.asarray(read_elem(z['obs/major'])).astype(str)
    mapped=[]
    for field in ['subCluster','subcluster']:
        values=pd.Series(read_elem(z['obs/'+field])).astype(str)
        lookup={}
        for name in names:
            mask=values.eq(name).to_numpy();ids=np.unique(nodes[mask])
            assert len(ids)==1 and set(major[mask])=={'T'}
            lookup[name]=str(ids[0])
        mapped.append(lookup)
    assert mapped[0]==mapped[1] and len(set(mapped[0].values()))==len(names)
    assert set(anchors)<=set(mapped[0].values())
    contexts={base[0]:None,base[1]:{anchors[0]},base[2]:{anchors[1]}}
    for label,node in mapped[0].items():
        if node not in anchors:contexts['CM02_network_member_'+label]={node}
    contexts['CM02_network_member_union']=set(mapped[0].values())
    return contexts,dict(network_membership_SHA256=sha(path),module='CM02',
        scope='T_whole_TME_primary_tumor',network_member_count=len(names),
        exact_authority_label_fields_agree=True,authority_major_T_verified=True,
        original_network_major_field_not_used_for_identity=True,
        network_labels_to_authority_ID=mapped[0],weights_not_read=True,
        full_network_members_not_the_full_nonzero_NMF_basis_support=True)


def parent_profiles(frame,depth,flags,parent,contexts,metrics,phases,depth_names):
    patients=sorted(frame.patient.unique());n=len(patients)
    patient=pd.Index(patients).get_indexer(frame.patient)
    phase=pd.Index(phases).get_indexer(frame.timepoint)
    key=phase*n+patient;result=[]
    for name,di in product(contexts,range(len(depth_names))):
        chosen=np.ones(len(frame),bool)
        nodes=contexts[name]
        if nodes is not None:chosen&=frame.node_id.isin(nodes).to_numpy()
        if di:chosen&=depth==di-1
        totals=np.bincount(key[chosen],minlength=len(phases)*n).reshape(len(phases),n)
        all_counts=np.column_stack([np.bincount(key[chosen],weights=flags[chosen,j],minlength=len(phases)*n) for j in range(len(metrics))]).astype(np.int64).reshape(len(phases),n,len(metrics))
        for pi,phase_name in enumerate(phases):
            rates={}
            for metric in ['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']:
                parent_name='regulatory_only' if metric.startswith('R_') else 'CD8_only'
                denominator=all_counts[pi,:,metrics.index(parent_name)]
                numerator=all_counts[pi,:,metrics.index(metric)]
                assert np.all(numerator<=denominator)
                good=denominator>=20
                values=numerator[good]/denominator[good]
                rates[metric]=dict(parent_observations=int(denominator.sum()),
                    state_positive_observations=int(numerator.sum()),
                    reviewed_cases_parent_atleast20=int(good.sum()),
                    case_state_fraction_min_q25_median_q75_max=np.quantile(values,[0,.25,.5,.75,1]).tolist() if good.sum()>=5 else None,
                    status='DESCRIPTIVE_WITHIN_PARENT_STATE_RATE' if good.sum()>=5 else 'INSUFFICIENT_CASE_SUPPORT')
            result.append(dict(source='GSE286827',phase=phase_name,context=name,
                capture_scope=depth_names[di],regulatory_parent=parent,
                cell_observations=int(totals[pi].sum()),parent_state_rates=rates,
                no_new_CM_correlation_clinical_or_spatial_test=True))
    return result
