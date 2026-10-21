"""Replay all nine primary spatial questions and three source-history scopes.

Uses independent tuple strata and 8x8 joint-state edge histograms for every
one of the 2000 randomizations; no cell or case arrays are persisted.
"""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from itertools import product
from pathlib import Path
import hashlib
import json
import os

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
import xenium_observed_state_context_v2 as context

B=2000
BITS=np.array([[(i>>j)&1 for j in range(3)] for i in range(8)],np.int64)


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def groups(data,masks,parent):
    xy=data['coordinates'][['x_centroid','y_centroid']].to_numpy(float)
    selected=masks[parent];p=xy[selected]
    density=np.array([len(v)-1 for v in cKDTree(xy[masks['T_support']]).query_ball_point(p,50)])
    capture=data['human_gene_molecules'][selected]
    db=np.digitize(density,np.unique(np.quantile(density,[1/3,2/3])),right=True)
    cb=np.digitize(capture,np.unique(np.quantile(capture,[1/3,2/3])),right=True)
    epcam=context.REQUIRED.index('EPCAM')
    epi=masks['qc']&~masks['T_support']&(data['counts'][:,epcam]>=1)
    close=cKDTree(xy[epi]).query(p,k=1)[0]<=50 if epi.any() else np.zeros(len(p),bool)
    by_key=defaultdict(list)
    for k,(point,d,c,e) in enumerate(zip(p,db,cb,close)):
        by_key[(int(np.floor(point[0]/500)),int(np.floor(point[1]/500)),int(d),int(c),int(e))].append(k)
    return [np.array(by_key[k],np.int64) for k in sorted(by_key)]


def probabilities(labels,group):
    result=np.empty(labels.shape,float)
    for g in group:result[g]=labels[g].sum(axis=0)/len(g)
    return result


def cross(left_code,right_code,rows,cols):
    joint=np.bincount(left_code[rows]*8+right_code[cols],minlength=64).reshape(8,8)
    return BITS.T@joint@BITS


def one(ctx,ordinal):
    meta=ctx['pairs'][ordinal]['source_metadata']
    if meta['diagnosis']!='HNSCC':return dict(ordinal=ordinal,excluded=True)
    d=context.load_pair(ctx,ordinal);m=context.observed_masks(d,20,'dual_regulatory_evidence')
    left=m['R_states'][m['regulatory']].astype(np.int64)
    right=m['E_states'][m['CD8']].astype(np.int64)
    lg=groups(d,m,'regulatory');rg=groups(d,m,'CD8')
    xy=d['coordinates'][['x_centroid','y_centroid']].to_numpy(float)
    neighbour_lists=cKDTree(xy[m['regulatory']]).query_ball_tree(cKDTree(xy[m['CD8']]),50)
    lengths=np.array([len(v) for v in neighbour_lists],np.int64)
    rows=np.repeat(np.arange(len(left)),lengths)
    cols=np.array([k for v in neighbour_lists for k in v],np.int64)
    graph=csr_matrix((np.ones(len(rows),np.int64),(rows,cols)),shape=(len(left),len(right)))
    lp=probabilities(left,lg);rp=probabilities(right,rg)
    independent_expected=lp[rows].T@rp[cols]
    # The independent edge expectation verifies the shared numerical sparse
    # convention; use that convention thereafter to avoid floating tie drift.
    expected=lp.T@(graph@rp)
    np.testing.assert_allclose(expected,independent_expected,rtol=1e-12,atol=1e-8)
    lc=left@np.array([1,2,4]);rc=right@np.array([1,2,4])
    observed=cross(lc,rc,rows,cols)
    np.testing.assert_array_equal(observed,left.T@(graph@right))
    draws=np.empty((B,3,3),np.int64);rng=np.random.default_rng([20261021,ordinal,0,1])
    for b in range(B):
        li=np.arange(len(left));ri=np.arange(len(right))
        for g in lg:li[g]=rng.permutation(g)
        for g in rg:ri[g]=rng.permutation(g)
        draws[b]=cross(lc[li],rc[ri],rows,cols)
    ready=(expected>0)&(np.ptp(draws,axis=0)>0)
    for i,j in product(range(3),range(3)):
        ready[i,j]&=len(left)>=20 and len(right)>=20 and 5<=left[:,i].sum()<len(left) and 5<=right[:,j].sum()<len(right)
    return dict(ordinal=ordinal,excluded=False,metadata=meta,ready=ready,
        effects=(observed-expected)/(expected+.5),
        null_effects=(draws-expected[None,...])/(expected[None,...]+.5),
        log2=np.log2((observed+.5)/(expected+.5)))


def compare(sections,original):
    proofs=[];computed=[]
    for si,scope in enumerate(['all_HNSCC','primary','recurrent']):
        current=[s for s in sections if scope=='all_HNSCC' or s['metadata']['recurrence']==scope]
        for i,j in product(range(3),range(3)):
            by_case=defaultdict(list)
            for s in current:
                if s['ready'][i,j]:by_case[s['metadata']['case']].append(s)
            target=next(r for r in original if r['scope']==scope and r['radius_um']==50
                and r['minimum_human_gene_molecules']==20 and r['regulatory_parent']=='dual_regulatory_evidence'
                and r['null']=='local_matched' and r['R_observed_state']==context.R_NAMES[i] and r['E_observed_state']==context.E_NAMES[j])
            assert target['supporting_cases']==len(by_case)
            if len(by_case)<5:
                assert target['status']=='INSUFFICIENT_SOURCE_CASES_FOR_FROZEN_SPATIAL_SUMMARY'
                proofs.append(dict(scope=scope,i=i,j=j,status='INSUFFICIENT_SCOPE_REPLAYED'));continue
            ordered=[by_case[k] for k in sorted(by_case)]
            effect=np.array([np.mean([s['effects'][i,j] for s in v]) for v in ordered])
            null=np.column_stack([np.mean([s['null_effects'][:,i,j] for s in v],axis=0) for v in ordered])
            logs=np.array([np.mean([s['log2'][i,j] for s in v]) for v in ordered])
            point=float(effect.mean());p=float((1+np.sum(np.abs(null.mean(axis=1))>=abs(point)))/(B+1))
            boot=np.random.default_rng([20261021,0,1,1,i,j,si]).integers(0,len(effect),size=(B,len(effect)))
            ci=np.quantile(effect[boot].mean(axis=1),[.025,.975])
            log_ci=np.quantile(logs[boot].mean(axis=1),[.025,.975])
            loo=np.array([(effect.sum()-v)/(len(effect)-1) for v in effect])
            assert abs(point-target['conditional_excess'])<1e-12
            assert abs(p-target['permutation_P'])<1e-15
            np.testing.assert_allclose(ci,target['CI95_case_bootstrap'],rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(log_ci,target['log2_CI95_case_bootstrap'],rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose([loo.min(),loo.max()],target['leave_one_case_range'],rtol=1e-12,atol=1e-12)
            assert int((effect>0).sum())==target['positive_case_effects']
            assert int((effect<0).sum())==target['negative_case_effects']
            assert bool(np.all(np.sign(loo)==np.sign(point)))==target['leave_one_all_same_sign_as_full']
            proofs.append(dict(scope=scope,i=i,j=j,status='ALL2000_NULL_BOOTSTRAP_LOO_AND_SIGNS_REPLAYED',
                point_absolute_error=abs(point-target['conditional_excess']),P_absolute_error=abs(p-target['permutation_P'])))
            if scope=='all_HNSCC':computed.append((p,target))
    ordered=sorted(computed,key=lambda p:p[0]);q=1.
    for rank in range(len(ordered),0,-1):
        p,t=ordered[rank-1];q=min(q,9*p/rank)
        assert abs(q-t['BH_nine_primary_opportunities'])<1e-15
    return proofs


def main():
    parent=Path(os.environ['PANCANCER_SPATIAL_PRODUCTION']);accepted=context.verified(parent)
    assert accepted['source_SHA256']==os.environ['PANCANCER_SPATIAL_SOURCE_SHA256']
    assert accepted['context_SHA256']==sha(context.__file__)
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
    save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('plan_before_primary_spatial_replay.json',dict(started_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(__file__),context_SHA256=sha(context.__file__),parent_terminal_SHA256=sha(parent/'private_terminal.json'),
        scope='All9 primary fixed radius/QC/parent/null, all/primary/recurrent source history. Independent tuple strata, state-pattern edge histograms, all2000 randomizations/bootstrap/LOO/BH.',
        no_new_effect_family_or_threshold_selection=True,no_cell_case_arrays_saved=True))
    original=json.loads((parent/'private_all648_source_spatial_opportunities.json').read_text())
    ctx=context.load_context();rows=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        work={pool.submit(one,ctx,i):i for i in range(18)}
        for f in as_completed(work):
            rows.append(f.result());save('private_checkpoint.json',dict(completed_source_slots=len(rows)))
    sections=sorted([s for s in rows if not s['excluded']],key=lambda s:s['ordinal'])
    assert len(rows)==18 and len(sections)==17
    proofs=compare(sections,original);assert len(proofs)==27
    save('private_primary_source_replay_proofs.json',proofs)
    save('private_terminal.json',dict(status='COMPLETE_ALL9_PRIMARY_SPATIAL_QUESTIONS_FULL_RANDOM_STREAM_AND27_SOURCE_SCOPES_REPLAYED',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),context_SHA256=sha(context.__file__),
        independent_tuple_strata_and_8x8_joint_state_edge_counts=True,
        all2000_nulls_bootstraps_leave_one_and_fixed_nine_BH_replayed=True,
        observed_and_analytic_expectation_independently_edge_checked=True,
        source_decoder_and_final_float_expectation_convention_shared_and_explicit=True,
        no_cell_case_arrays_or_new_effect_family_saved=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_PRIMARY_SPATIAL_INDEPENDENT_REPLAY',terminal_SHA256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
