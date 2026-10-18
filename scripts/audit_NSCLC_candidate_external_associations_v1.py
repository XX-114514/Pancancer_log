"""Complete OLS replay of every external effect, random stream and correction.

The OLS oracle is the previously accepted source-coordination audit, with the
same predeclared full-span stability guard. No source/candidate selection.
"""
import hashlib
import json
import os
import time
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist,spearmanr

import audit_NSCLC_external_programme_coordination_v3 as oracle
from test_NSCLC_candidate_external_associations_v1 import inputs,OUT
from prepare_NSCLC_candidate_external_associations_v1 import QCOLS
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def seed(model,stream):
    return int.from_bytes(hashlib.sha256((model+'|'+stream+'|external_fixed_native_ligand_programme_20261018_v1').encode()).digest()[:8],'little')


def fit(a,b,q,blocks,minimum=20):
    p=oracle.fit(a,b,q,blocks,minimum=minimum)
    s=np.linalg.svd(p['z'],compute_uv=False)
    rank=int((s>s[0]*np.finfo(float).eps*max(p['z'].shape)).sum())
    if rank!=p['rank'] or s[p['rank']-1]/s[0]<1e-8:raise ValueError('UNSTABLE_FULL_NUISANCE_SPAN')
    return p


def secondary(a,b,q,blocks):
    try:return fit(a,b,q,blocks)['rho'],'ESTIMABLE'
    except ValueError as e:return np.nan,str(e)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1024**2)
    job=root/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_candidate_external_associations_v1'
    launch=json.loads((job/'launch_receipt.json').read_text());out=root/OUT;deadline=time.monotonic()+21600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'Bound inferential producer missing without terminal'
        f=p.read_text().rsplit(')',1)[1].split();assert f[0]!='Z' and f[19]==str(launch['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired, no restart inferred';time.sleep(3)
    terminal=json.loads((out/'private_terminal.json').read_text());assert terminal['source_sha256']==launch['source_sha256']
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    src,proof,contract,eligibility,cases,frozen=inputs(root)
    active=eligibility[eligibility.test_index.ge(0)]
    result=pd.read_csv(out/'private_all_estimable_external_effects_and_sensitivities.tsv.gz',sep='\t',keep_default_na=False,na_values=[''])
    planned=pd.read_csv(out/'private_complete_planned_eligibility.tsv.gz',sep='\t',keep_default_na=False)
    pd.testing.assert_frame_equal(planned,eligibility)
    assert result.model.tolist()==active.model.tolist()
    for col in active.columns:
        if col in ['candidate_units','active_units','active_reviewed_units','pooled_sender_positive_fraction','sender_positive_units','nuisance_rank','residual_df','test_index']:
            assert np.allclose(pd.to_numeric(result[col]),pd.to_numeric(active[col]),atol=1e-12,rtol=1e-12,equal_nan=True)
        else:assert np.array_equal(result[col].astype(str),active[col].astype(str))
    with np.load(out/'private_all_bootstrap_and_LOO.npz',allow_pickle=False) as f:streams={k:f[k] for k in f.files}
    assert streams['bootstrap'].shape==(len(result),2000)
    maximum_effect_error=maximum_boot_error=0.;loo_checked=0
    checkpoint=root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for row in result.itertuples():
        g=cases.iloc[json.loads(row.active_case_rows_json)]
        x,y,q,blocks=g.ligand.to_numpy(float),g.response.to_numpy(float),g[QCOLS].to_numpy(float),g.block_key.to_numpy(str)
        p=fit(x,y,q,blocks);assert p['keep'].all() and p['rank']==row.nuisance_rank
        assert np.isclose(p['rho'],row.partial_rank_rho,atol=1e-10,rtol=1e-10)
        maximum_effect_error=max(maximum_effect_error,abs(p['rho']-row.partial_rank_rho))
        groups=[np.flatnonzero(blocks==b) for b in sorted(set(blocks))]
        inverse=np.linalg.pinv(p['z'],rcond=1e-10);mu=p['y'][:,1]-p['residual'][:,1]
        rng=np.random.default_rng(seed(row.model,'permutation'));count=0;draws=contract['permutations']
        for start in range(0,draws,256):
            m=min(256,draws-start);ix=np.broadcast_to(np.arange(len(g)),(m,len(g))).copy()
            for group in groups:ix[:,group]=group[np.argsort(rng.random((m,len(group))),axis=1)]
            response=mu[None,:]+p['residual'][ix,1];residual=response-(response@inverse.T)@p['z'].T
            values=(residual@p['residual'][:,0])/(np.linalg.norm(residual,axis=1)*np.linalg.norm(p['residual'][:,0]))
            assert np.isfinite(values).all();count+=int((abs(values)>=abs(p['rho'])-1e-12).sum())
        assert count==row.permutation_exceedances and np.isclose(row.P,(1+count)/(1+draws),atol=1e-15)
        size=frozen['eligible_family_sizes'][row.family]
        upper=float(beta_dist.ppf(1-.01/size,count+1,draws-count)) if count<draws else 1.
        assert np.isclose(row.MC_upper_P,upper,atol=1e-12,rtol=1e-12)
        rng=np.random.default_rng(seed(row.model,'bootstrap'));boot=np.full(2000,np.nan)
        for j in range(2000):
            ix=np.concatenate([group[rng.integers(len(group),size=len(group))] for group in groups])
            try:boot[j]=fit(x[ix],y[ix],q[ix],blocks[ix])['rho']
            except ValueError:pass
        ref=streams['bootstrap'][row.test_index]
        assert np.allclose(boot,ref,atol=1e-9,rtol=1e-9,equal_nan=True)
        finite=np.isfinite(boot);assert int(finite.sum())==row.valid_bootstraps
        ci=np.quantile(boot[finite],[.025,.975]) if finite.sum()>=1900 else [np.nan,np.nan]
        assert np.allclose(ci,[row.CI_low,row.CI_high],atol=1e-9,rtol=1e-9,equal_nan=True)
        if finite.any():maximum_boot_error=max(maximum_boot_error,float(np.max(abs(boot[finite]-ref[finite]))))
        loo=np.full(len(g),np.nan)
        for j in range(len(g)):
            take=np.arange(len(g))!=j
            try:loo[j]=fit(x[take],y[take],q[take],blocks[take],minimum=3)['rho']
            except ValueError:pass
        assert np.allclose(loo,streams['loo_'+str(row.test_index)],atol=1e-9,rtol=1e-9,equal_nan=True)
        finite=np.isfinite(loo);assert int(finite.sum())==row.LOO_valid and bool(finite.all())==row.LOO_all_valid
        assert bool(finite.all() and (np.sign(loo)==np.sign(p['rho'])).all())==row.LOO_same_sign
        bounds=[float(np.min(loo[finite])),float(np.max(loo[finite]))] if finite.any() else [np.nan,np.nan]
        assert np.allclose(bounds,[row.LOO_min,row.LOO_max],atol=1e-9,rtol=1e-9,equal_nan=True)
        loo_checked+=len(loo)
        assert np.isclose(row.raw_Spearman,spearmanr(x,y).statistic,atol=1e-12)
        for label,qq in [('block_only',np.empty((len(g),0))),('ten_QC_without_parent_fractions',q[:,:10]),
                         ('receiver_own_ligand_adjusted',np.column_stack([q,g.receiver_own_ligand.to_numpy(float)])),
                         ('whole_annotation_fractions_adjusted',np.column_stack([q,g[['fraction_whole_sender','fraction_whole_receiver']].to_numpy(float)]))]:
            value,status=secondary(x,y,qq,blocks)
            assert getattr(row,label+'_status')==status and np.isclose(getattr(row,label+'_rho'),value,atol=1e-9,rtol=1e-9,equal_nan=True)
        checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed=row.test_index+1,total=len(result),stage='complete_external_OLS_and_stream_replay'))+'\n')
    for family,size in frozen['eligible_family_sizes'].items():
        take=result.family.eq(family);assert take.sum()==size
        if size:
            for col,source,by in [('BH_q','P',False),('BY_q','P',True),('MC_upper_BH_q','MC_upper_P',False)]:
                assert np.allclose(result.loc[take,col],oracle.fdr(result.loc[take,source].to_numpy(float),by),atol=1e-12,rtol=1e-12)
    report=dict(status='PASS_ALL_FROZEN_EXTERNAL_CANDIDATE_EFFECTS_FULL_STREAMS_SENSITIVITIES_AND_FAMILY_CORRECTIONS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),OLS_oracle_sha256=sha(oracle.__file__),
        producer_terminal_sha256=sha(out/'private_terminal.json'),all_planned_positions=len(eligibility),all_estimable_models=len(result),
        permutations_replayed=len(result)*contract['permutations'],bootstraps_replayed=len(result)*2000,LOO_replayed=loo_checked,
        maximum_effect_error=maximum_effect_error,maximum_bootstrap_error=maximum_boot_error,
        limitations='Same accepted paired-case I/O and exact integer-rank primitives; independent OLS projection, full random streams and multiplicity. Not untouched independent cohorts, direct communication or clinical validation.')
    p=out/'private_independent_external_candidate_result_audit.json';assert not p.exists();p.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],completed_UTC=report['completed_UTC'],proof_sha256=sha(p))))


if __name__=='__main__':main()
