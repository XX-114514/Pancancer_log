"""Finite author-type tissue descriptions; no P, FDR or patient inference.

All individual arrays/edges remain source RAM. Only whole-source equal-image
summaries leave the process. Coordinates have unknown physical units.
"""
from collections import Counter
import hashlib
import resource
import numpy as np
from scipy.spatial import cKDTree


def fingerprint(payloads):
    digests=[]
    for p in payloads:
        d=hashlib.sha256()
        for a in (p['xy'],p['codes']['anno'],p['codes']['anno2'],p['codes']['bin50 structural annotation'],p['usable']):
            d.update(np.ascontiguousarray(a).tobytes())
        digests.append(d.hexdigest())
    return hashlib.sha256(''.join(digests).encode()).hexdigest()


def interval(values):
    values=np.asarray(values,float);values=values[np.isfinite(values)]
    if not len(values):return dict(estimable_images=0,median=None,IQR=None,positive_images=0,negative_images=0,leave_one_image_out_median_range=None)
    loo=[float(np.median(np.delete(values,i))) for i in range(len(values))] if len(values)>1 else []
    return dict(estimable_images=len(values),median=float(np.median(values)),IQR=np.quantile(values,[.25,.75]).tolist(),
                positive_images=int((values>0).sum()),negative_images=int((values<0).sum()),
                leave_one_image_out_median_range=[min(loo),max(loo)] if loo else None)


def odds(table):
    a,b,c,d=np.asarray(table,float)
    return float(np.log((a+.5)*(d+.5)/((b+.5)*(c+.5))))


def mh(tables):
    numerator=denominator=0.;accepted=0
    for t in tables:
        a,b,c,d=map(float,t);n=a+b+c+d
        if n<10 or a+b==0 or c+d==0:continue
        numerator+=a*d/n;denominator+=b*c/n;accepted+=1
    return (float(np.log(numerator/denominator)) if numerator>0 and denominator>0 else None),accepted


def graph(payload,k):
    ids=np.flatnonzero(payload['usable']);points=payload['xy'][ids];tree=cKDTree(points)
    assert len(ids)>k
    nearest=tree.query(points,k=2,workers=1)[0][:,1]
    gap=5*float(np.median(nearest));assert gap>0
    anchors=np.flatnonzero(payload['cd4']);query=tree.query(payload['xy'][anchors],k=k+1,workers=1)[0][:,-1]
    radii=np.minimum(query,gap)
    lists=tree.query_ball_point(payload['xy'][anchors],np.nextafter(radii,np.inf),workers=1)
    edges=[]
    for anchor, neighbors in zip(anchors,lists):
        targets=ids[np.asarray(neighbors,int)];targets=targets[targets!=anchor]
        edges.append((int(anchor),targets))
    for position in sorted(set([0,len(anchors)//2,len(anchors)-1])) if len(anchors) else []:
        anchor=anchors[position];dist=np.sqrt(np.sum((points-payload['xy'][anchor])**2,axis=1))
        dist[ids==anchor]=np.inf
        radius=min(float(np.partition(dist,k-1)[k-1]),gap)
        direct=ids[dist<=np.nextafter(radius,np.inf)]
        assert np.array_equal(np.sort(edges[position][1]),np.sort(direct))
    return edges


def one_pair(payload,edges,label,minimum):
    categories=payload['categories'];codes=payload['codes'];fine=categories['anno2'];fc=codes['anno2']
    target=payload['cd8'] & (fc==fine.index(label)) if label in fine else np.zeros(len(fc),bool)
    other=payload['cd8'] & ~target;R=payload['treg'];Rother=payload['cd4'] & ~R
    groupcounts=[int(v.sum()) for v in (R,Rother,target,other)]
    if min(groupcounts)<minimum:return dict(status='PARENT_OR_SUBTYPE_BELOW20',value=None,region_value=None,region_status='PARENT_OR_SUBTYPE_BELOW20')
    table=np.zeros(4,np.int64);byregion={};anchorcoverage=np.zeros(2,int);regioncoverage={}
    rc=codes['bin50 structural annotation'];accepted_edges=[];accepted_region_edges=[]
    for anchor, neighbors in edges:
        isR=bool(R[anchor]);row=0 if isR else 2
        e=int(target[neighbors].sum());other_e=int(other[neighbors].sum())
        table[row]+=e;table[row+1]+=other_e
        if e+other_e>0:anchorcoverage[0 if isR else 1]+=1
        accepted_edges.extend((isR,bool(target[j])) for j in neighbors if payload['cd8'][j])
        region=int(rc[anchor])
        if region<0:continue
        same=neighbors[rc[neighbors]==region]
        e=int(target[same].sum());other_e=int(other[same].sum())
        t=byregion.setdefault(region,np.zeros(4,np.int64));t[row]+=e;t[row+1]+=other_e
        ac=regioncoverage.setdefault(region,np.zeros(2,int))
        if e+other_e>0:ac[0 if isR else 1]+=1
        accepted_region_edges.extend((region,isR,bool(target[j])) for j in same if payload['cd8'][j])
    replay=Counter(accepted_edges)
    assert np.array_equal(table,[replay[(True,True)],replay[(True,False)],replay[(False,True)],replay[(False,False)]])
    region_replay=Counter(accepted_region_edges)
    for region,t in byregion.items():
        assert np.array_equal(t,[region_replay[(region,True,True)],region_replay[(region,True,False)],
                                 region_replay[(region,False,True)],region_replay[(region,False,False)]])
    if min(anchorcoverage)<10:return dict(status='ANCHORS_WITH_CD8_NEIGHBORS_BELOW10',value=None,region_value=None,region_status='ANCHORS_WITH_CD8_NEIGHBORS_BELOW10')
    region_count=np.sum(list(regioncoverage.values()),axis=0) if regioncoverage else np.zeros(2,int)
    region_value,strata=mh(list(byregion.values()))
    if min(region_count)<10:region_value=None;region_status='SAME_REGION_ANCHORS_BELOW10'
    elif region_value is None:region_status='NO_FINITE_MH_ODDS'
    else:region_status='ESTIMABLE'
    return dict(status='ESTIMABLE',value=odds(table),region_value=region_value,region_status=region_status,
                MH_eligible_regions=strata)


def regional_composition(payload,region,label,minimum):
    labels=payload['categories']['bin50 structural annotation'];rc=payload['codes']['bin50 structural annotation']
    inside=(rc==labels.index(region)) if region in labels else np.zeros(len(rc),bool)
    parent=payload['cd4'] if label=='CD4_Treg' else payload['cd8']
    fine=payload['categories']['anno2'];fc=payload['codes']['anno2']
    target=parent & (fc==fine.index(label)) if label in fine else np.zeros(len(rc),bool)
    if min(int((parent&inside).sum()),int((parent&~inside).sum()))<minimum:return None
    return odds([int((target&inside).sum()),int((parent&inside&~target).sum()),
                 int((target&~inside).sum()),int((parent&~inside&~target).sum())])


def run_description(payloads,contract):
    synthetic()
    assert contract['no_P_CI_FDR_or_patient_claim'] and contract['minimum_fourgroup']==20
    assert contract['k_values']==[5,10,20] and contract['primary_k']==10
    assert all(label.startswith('CD8_') for label in contract['all_CD8_labels'])
    assert fingerprint(payloads)==contract['combined_source_RAM_array_SHA256']
    all_pairs=[];regional=[]
    for k in contract['k_values']:
        records={label:[] for label in contract['all_CD8_labels']}
        for payload in payloads:
            edges=graph(payload,k)
            for label in records:records[label].append(one_pair(payload,edges,label,20))
            assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < contract['RSS_cap_KiB']
        for label,items in records.items():
            all_pairs.append(dict(CD8_author_label=label,k=k,primary=(k==10),
                global_native_graph=interval([v['value'] for v in items if v['value'] is not None]),
                same_author_region_MH=interval([v['region_value'] for v in items if v['region_value'] is not None]),
                missing_reasons=dict(Counter(v['status'] for v in items if v['value'] is None)),
                region_missing_reasons=dict(Counter(v['region_status'] for v in items if v['region_value'] is None)),
                ambiguous_CD8_Treg=(label=='CD8_Treg')))
    for region in contract['all_qualified_region_labels']:
        for label in ['CD4_Treg']+contract['all_CD8_labels']:
            vals=[regional_composition(p,region,label,20) for p in payloads]
            regional.append(dict(author_region_label=region,author_type=label,
                relative_to_same_parent_elsewhere=interval([v for v in vals if v is not None]),
                unavailable_images=sum(v is None for v in vals)))
    assert len(all_pairs)==3*len(contract['all_CD8_labels'])
    return dict(status='COMPLETE_AUTHOR_TYPE_SPATIAL_DESCRIPTION_NO_INFERENCE',source='GSE328481',
        image_count=len(payloads),all_spatial_slots=all_pairs,all_regional_composition_slots=regional,
        independent_flat_edge_replay_all_estimable_designs=True,
        direct_dense_neighbor_checks_first_middle_last_CD4_per_image_all_k=True,
        pseudocount_half_for_crude_and_region_composition_only=True,MH_no_pseudocount=True,
        image_median_IQR_and_LOO_are_not_patient_intervals=True,
        annotation_transfer_not_RNA_state_or_original_CM_exact_replication=True,
        selected_immature_TLS_source_not_population_estimate=True,
        no_clinical_patient_phase_physical_units_or_direct_contact_claim=True)


def synthetic():
    assert abs(odds([2,1,1,2])-np.log(25/9))<1e-14
    assert odds([0,0,0,0])==0
    value,count=mh([[2,1,1,2],[4,2,2,4]])
    # First stratum is below the fixed10 edge threshold and remains excluded.
    assert count==1 and abs(value-np.log(4))<1e-14
    assert interval([])['median'] is None
    v=interval([-2,0,2]);assert v['median']==0 and v['positive_images']==1 and v['negative_images']==1
    p=dict(usable=np.ones(25,bool),xy=np.column_stack([np.arange(25),np.zeros(25)]),cd4=np.array([True]+[False]*24))
    e=graph(p,5);assert e[0][0]==0 and list(e[0][1])==list(range(1,6))
    p=dict(usable=np.ones(5,bool),xy=np.array([[0,0],[1,0],[-1,0],[0,1],[0,-1.]]),cd4=np.array([True,False,False,False,False]))
    assert len(graph(p,1)[0][1])==4  # All ties retained, never arbitrary first.


def full_family_synthetic():
    n=400;cd4=np.zeros(n,bool);cd4[:200:2]=True
    cd8=np.zeros(n,bool);cd8[1:200:2]=True
    R=cd4 & (np.arange(n)<100);fc=np.full(n,4,dtype=np.int16)
    fc[R]=0;fc[cd4 & ~R]=1;fc[cd8 & (np.arange(n)<100)]=2;fc[cd8 & (np.arange(n)>=100)]=3
    p=dict(xy=np.column_stack([np.arange(n),np.zeros(n)]),usable=np.ones(n,bool),
           cd4=cd4,cd8=cd8,treg=R,
           categories={'anno2':['CD4_Treg','CD4_IL7R','CD8_Tex','CD8_Tem','cancer_cell'],
                       'bin50 structural annotation':['Tumor','stroma']},
           codes={'anno':np.zeros(n,np.int16),'anno2':fc,'bin50 structural annotation':(np.arange(n)>=100).astype(np.int16)})
    opposite={**p,'codes':{**p['codes'],'anno2':fc.copy(),'bin50 structural annotation':np.zeros(n,np.int16)}}
    opposite['codes']['anno2'][cd8]=5-fc[cd8]
    missing={**p,'treg':np.zeros(n,bool)}
    contract=dict(no_P_CI_FDR_or_patient_claim=True,minimum_fourgroup=20,k_values=[5,10,20],primary_k=10,
                  all_CD8_labels=['CD8_Tex','CD8_Tem'],all_qualified_region_labels=['Tumor','stroma'],RSS_cap_KiB=1024**2,
                  combined_source_RAM_array_SHA256=fingerprint([p,opposite,missing]))
    out=run_description([p,opposite,missing],contract)
    assert len(out['all_spatial_slots'])==6 and len(out['all_regional_composition_slots'])==6
    for slot in out['all_spatial_slots']:
        v=slot['global_native_graph']
        assert v['estimable_images']==2 and v['positive_images']==1 and v['negative_images']==1
        assert abs(v['median'])<1e-14 and slot['missing_reasons']=={'PARENT_OR_SUBTYPE_BELOW20':1}
        assert slot['same_author_region_MH']['estimable_images']==1
        assert slot['region_missing_reasons']=={'NO_FINITE_MH_ODDS':1,'PARENT_OR_SUBTYPE_BELOW20':1}


if __name__=='__main__':
    synthetic();full_family_synthetic();print('PASS_TIE_SAFE_GEOMETRY_ODDS_MH_AND_FULL_FAMILY')
