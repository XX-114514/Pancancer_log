"""Full independent OLS replay of fixed native interface effects and streams."""
import hashlib
import json
import multiprocessing
import os
import time
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist,rankdata,spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_external_exact_rank_guards_v1 import rational_rank,any_exactly_fitted
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_native_interface_association_v1 import inputs

STATE=None


def seed(record,stream):
    return int.from_bytes(hashlib.sha256((record+'|'+stream+'|native_continuous_prior_association_20261016_v1').encode()).digest()[:8],'little')


def fit(a,b,quality,blocks,minimum=50):
    block=np.asarray(blocks,str);counts=pd.Series(block).value_counts()
    keep=np.array([counts[s]>=2 for s in block],bool)
    if keep.sum()<minimum:raise ValueError('INSUFFICIENT_INFORMATIVE_UNITS')
    a,b,q,block=np.asarray(a)[keep],np.asarray(b)[keep],np.asarray(quality)[keep],block[keep]
    if not(np.isfinite(a).all() and np.isfinite(b).all()):raise ValueError('NONFINITE_PROGRAMME')
    if not np.isfinite(q).all():raise ValueError('NONFINITE_QUALITY')
    qr=rankdata(q,axis=0);sd=qr.std(axis=0)
    dummies=[(block==s).astype(int) for s in sorted(set(block))]
    z=np.column_stack([np.ones(len(a)),*dummies,(qr-qr.mean(axis=0))/np.where(sd>0,sd,1)])
    y=rankdata(np.column_stack([a,b]),axis=0);coef,_,rank,singular=np.linalg.lstsq(z,y,rcond=1e-10)
    if len(a)-rank<10:raise ValueError('QUALITY_AND_BLOCK_DESIGN_INSUFFICIENT_RESIDUAL_DF')
    residual=y-z@coef;norm=np.linalg.norm(residual,axis=0)
    if norm.min()<1e-6:
        zi=np.column_stack([np.ones(len(a)),*dummies,2*qr]).astype(np.int64)
        if rational_rank(zi)==rank and any_exactly_fitted(zi,(2*y).astype(np.int64)):raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    if (norm<1e-10).any():raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    machine=(singular>singular[0]*np.finfo(float).eps*max(z.shape)).sum()
    if rank!=machine or singular[rank-1]/singular[0]<1e-8:raise ValueError('UNSTABLE_FULL_NUISANCE_SPAN')
    return dict(a=a,b=b,quality=q,blocks=block,keep=keep,z=z,y=y,residual=residual,rank=int(rank),rho=float(np.clip(np.corrcoef(residual.T)[0,1],-1,1)))


def secondary(a,b,q,block):
    try:return fit(a,b,q,block)['rho'],'ESTIMABLE'
    except ValueError as error:return np.nan,str(error)


def make_context(data,z,clinical):
    cases=data['cases'].sampleID.to_numpy(str);case_index={s:i for i,s in enumerate(cases)};block=data['cases'].holdout_block.to_numpy(str)
    parts={};composition={};total=np.zeros(len(cases),float)
    for u in data['units'].itertuples():
        if u.sampleID in case_index:total[case_index[u.sampleID]]+=u.n_cells
    for sv,sender in enumerate(VIEWS):
        labels={s:i for i,s in enumerate(data['labels'][sender])};counts=np.zeros((len(cases),len(labels)),float)
        for u in data['units'].itertuples():
            if u.sampleID in case_index and u.sub_cell_type in labels:counts[case_index[u.sampleID],labels[u.sub_cell_type]]+=u.n_cells
        assert np.array_equal(counts.sum(axis=1),z['cells'][sv]);composition[sender]=counts/z['cells'][sv,:,None]
        for rv,receiver in enumerate(VIEWS):
            if sv==rv:continue
            quality=np.column_stack([clinical,data['prepared'][sv]['quality'],data['prepared'][rv]['quality']])
            complete=np.isfinite(quality).all(axis=1)
            for name,mask in [('internal_holdout',~data['train']),('discovery',data['train']),('all_frozen_cases',np.ones(len(cases),bool))]:
                ix=np.flatnonzero(complete&mask);take=[i for i in ix if np.sum(block[ix]==block[i])>=2];ix=np.array(take,int)
                parts[(sender,receiver,name)]=(ix,quality[ix],block[ix])
    assert (total>0).all()
    return dict(data=data,z=z,parts=parts,composition=composition,fractions=z['cells']/total[None,:],views={v:i for i,v in enumerate(VIEWS)})


def check_one(item):
    state=STATE;z=state['z'];sv=state['views'][item['sender_view']];rv=state['views'][item['receiver_view']]
    slot=item['slot_index'];j=item['ligand_axis_index'];direction=item['direction'];k=int(item['factor'].rsplit('F',1)[1])-1
    ix,q,block=state['parts'][(item['sender_view'],item['receiver_view'],'internal_holdout')]
    x=z['ligand_mean'][sv,ix,j];y=direction*z['response'][slot,ix];p=fit(x,y,q,block)
    i=item['test_index'];row=state['results'].iloc[i]
    assert row.record==item['record'] and row.active_cases==len(ix) and row.nuisance_rank==p['rank'] and row.residual_df==len(ix)-p['rank']
    assert np.isclose(row.partial_rank_rho,p['rho'],atol=1e-10,rtol=1e-10)
    inverse=np.linalg.pinv(p['z'],rcond=1e-10);mu=p['y'][:,1]-p['residual'][:,1];groups=[np.flatnonzero(block==s) for s in sorted(set(block))]
    rng=np.random.default_rng(seed(row.record,'Freedman_Lane'));count=0
    for start in range(0,199999,256):
        m=min(256,199999-start);take=np.broadcast_to(np.arange(len(ix)),(m,len(ix))).copy()
        for g in groups:take[:,g]=g[np.argsort(rng.random((m,len(g))),axis=1)]
        yy=mu[None,:]+p['residual'][take,1];error=yy-(yy@inverse.T)@p['z'].T
        values=(error@p['residual'][:,0])/(np.linalg.norm(error,axis=1)*np.linalg.norm(p['residual'][:,0]))
        count+=int((abs(values)>=abs(p['rho'])-1e-12).sum())
    assert count==row.permutation_exceedances and np.isclose(row.P,(1+count)/200000,atol=1e-15)
    upper=1. if count==199999 else float(beta_dist.ppf(1-.01/state['family_size'],count+1,199999-count))
    assert np.isclose(row.MC_upper_P,min(1.,(1+199999*upper)/200000),atol=1e-10,rtol=1e-10)
    rng=np.random.default_rng(seed(row.record,'bootstrap'));boot=np.full(2000,np.nan)
    for b in range(2000):
        take=np.concatenate([g[rng.integers(len(g),size=len(g))] for g in groups])
        try:boot[b]=fit(x[take],y[take],q[take],block[take])['rho']
        except ValueError:pass
    reference=state['bootstrap'][i];assert np.allclose(boot,reference,atol=1e-9,rtol=1e-9,equal_nan=True)
    finite=np.isfinite(boot);assert row.valid_bootstraps==finite.sum()
    ci=np.quantile(boot[finite],[.025,.975]) if finite.sum()>=1900 else [np.nan,np.nan]
    assert np.allclose([row.CI_low,row.CI_high],ci,atol=1e-9,rtol=1e-9,equal_nan=True)
    loo=np.full(len(ix),np.nan)
    for n in range(len(ix)):
        keep=np.arange(len(ix))!=n
        try:loo[n]=fit(x[keep],y[keep],q[keep],block[keep],minimum=3)['rho']
        except ValueError:pass
    offsets=state['LOO_offsets'][i:i+2]
    assert np.allclose(loo,state['LOO'][offsets[0]:offsets[1]],atol=1e-9,rtol=1e-9,equal_nan=True)
    assert row.LOO_valid==np.isfinite(loo).sum() and row.LOO_all_valid==np.isfinite(loo).all()
    assert row.LOO_same_sign==bool(np.isfinite(loo).all() and (np.sign(loo)==np.sign(p['rho'])).all())
    limits=[np.nanmin(loo),np.nanmax(loo)] if np.isfinite(loo).any() else [np.nan,np.nan]
    assert np.allclose([row.LOO_min,row.LOO_max],limits,atol=1e-9,rtol=1e-9,equal_nan=True)
    assert np.isclose(row.raw_Spearman,spearmanr(x,y).statistic,atol=1e-12)
    for part in ['discovery','all_frozen_cases']:
        ii,qq,bb=state['parts'][(item['sender_view'],item['receiver_view'],part)]
        value,status=secondary(z['ligand_mean'][sv,ii,j],direction*z['response'][slot,ii],qq,bb)
        assert getattr(row,part+'_status')==status and getattr(row,part+'_cases')==len(ii)
        assert np.isclose(getattr(row,part+'_rho'),value,atol=1e-9,rtol=1e-9,equal_nan=True)
    extra=np.column_stack([state['composition'][item['sender_view']][ix],state['composition'][item['receiver_view']][ix]])
    for name,yy,qq in [('original_receiver_programme',direction*z['original_scores'][rv,ix,k],q),('sender_programme_adjusted',y,np.column_stack([q,z['original_scores'][sv,ix,k]])),('native_immune_abundance_adjusted',y,np.column_stack([q,state['fractions'][[sv,rv]][:,ix].T])),('complete_author_subtype_adjusted',y,np.column_stack([q,extra]))]:
        value,status=secondary(x,yy,qq,block)
        assert getattr(row,name+'_status')==status and np.isclose(getattr(row,name+'_rho'),value,atol=1e-9,rtol=1e-9,equal_nan=True)
    valid=finite&np.isfinite(reference);error=float(np.max(abs(boot[valid]-reference[valid]))) if valid.any() else 0.
    return dict(maximum_effect_error=abs(row.partial_rank_rho-p['rho']),maximum_bootstrap_error=error,LOO_replayed=len(ix))


def fdr(p,by=False):
    order=np.argsort(p);n=len(p);scale=np.sum(1/np.arange(1,n+1)) if by else 1.
    values=np.minimum.accumulate((p[order]*n*scale/np.arange(1,n+1))[::-1])[::-1];out=np.empty(n);out[order]=np.minimum(values,1.)
    return out


def main():
    global STATE
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);src=root/'NSCLC_native_interface_association_results_20261018_v1'
    job=root/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_native_interface_association_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+43200
    while not (src/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists(),'Bound production missing without terminal'
        values=p.read_text().rsplit(')',1)[1].split();assert values[0]!='Z' and values[19]==str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired; no production termination inferred';time.sleep(5)
    resource(root,1048576);terminal=json.loads((src/'private_terminal.json').read_text())
    assert terminal['source_sha256']==receipt['source_sha256']
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    data=load();pre,proof,contract,z,clinical,scope=inputs(root);plan=json.loads((src/'plan_before_all_native_interface_effects.json').read_text())
    assert plan['contract']==contract and plan['input_audit_sha256']==sha(proof)
    active=scope[scope.test_index.ge(0)];result=pd.read_csv(src/'private_all_primary_and_sensitivity_effects.tsv.gz',sep='\t')
    assert result.record.tolist()==active.record.tolist() and len(result)==plan['family_size_before_effects']
    for col in active.columns.drop('status'):pd.testing.assert_series_equal(result[col],active[col].reset_index(drop=True),check_dtype=False,check_names=False)
    with np.load(src/'private_all_bootstrap_and_LOO.npz',allow_pickle=False) as saved:streams={key:saved[key] for key in saved.files}
    assert streams['records'].tolist()==result.record.tolist() and streams['bootstrap'].shape==(len(result),2000)
    STATE=make_context(data,z,clinical);STATE.update(results=result,family_size=len(result),**streams)
    completed=[];checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('fork')) as executor:
        futures=[executor.submit(check_one,item) for item in active.to_dict('records')]
        for future in as_completed(futures):
            completed.append(future.result())
            if len(completed)%10==0 or len(completed)==len(result):checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_primary_tests=len(completed),total_primary_tests=len(result),stage='independent_full_random_replay',workers=4))+'\n')
    for target,source,by in [('BH_q','P',False),('BY_q','P',True),('MC_upper_BH_q','MC_upper_P',False)]:assert np.allclose(result[target],fdr(result[source].to_numpy(),by),atol=1e-12,rtol=1e-12)
    candidate=result.partial_rank_rho.gt(0)&result.BH_q.lt(.05)&result.MC_upper_BH_q.lt(.05)&result.CI_low.gt(0)&result.LOO_same_sign&result.discovery_rho.gt(0)
    assert np.array_equal(result.internal_candidate_under_prespecified_rules,candidate)
    report=dict(status='PASS_COMPLETE_FROZEN_NATIVE_INTERFACE_FAMILY_ALL_REAL_EFFECTS_STREAMS_SENSITIVITIES_AND_CORRECTIONS',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),all_planned_positions=4800,primary_tests_checked=len(result),permutations_replayed=len(result)*199999,bootstraps_replayed=len(result)*2000,LOO_replayed=sum(r['LOO_replayed'] for r in completed),maximum_effect_error=max(r['maximum_effect_error'] for r in completed),maximum_bootstrap_error=max(r['maximum_bootstrap_error'] for r in completed),clinical_outcomes_or_direct_communication_validated=False)
    proof=src/'private_independent_native_interface_result_audit.json';assert not proof.exists();proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
