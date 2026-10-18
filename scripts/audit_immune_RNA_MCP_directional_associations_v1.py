"""All directional positions, complete random streams and independent batched SVD."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from audit_immune_RNA_MCP_stability_and_coordination_v1 import independent_fit, family_q, rho, equal


def boot_svd(x,y,C,takes):
    """Independent pandas tied ranks and batched full SVD nuisance projection."""
    count,n=takes.shape;k=C.shape[1]+1
    xx=x[takes];yy=y[takes];cc=C[takes]
    status=np.full(count,'ESTIMABLE',dtype='<U90')
    constant=(np.ptp(xx,axis=1)==0)|(np.ptp(yy,axis=1)==0);status[constant]='CONSTANT_SCORE'
    rankc=np.stack([pd.DataFrame(cc[:,:,j]).rank(axis=1).to_numpy() for j in range(C.shape[1])],axis=2)
    sd=rankc.std(axis=1,ddof=1);constcov=np.any(sd==0,axis=1)
    status[(status=='ESTIMABLE')&constcov]='CONSTANT_REQUIRED_COVARIATE_NO_AUTOMATIC_DELETION'
    scaled=(rankc-rankc.mean(axis=1)[:,None,:])/np.where(sd==0,1,sd)[:,None,:]
    D=np.concatenate([np.ones((count,n,1)),scaled],axis=2)
    U,S,_=np.linalg.svd(D,full_matrices=False)
    rank=(S>S[:,0,None]*max(n,k)*np.finfo(float).eps).sum(axis=1)
    status[(status=='ESTIMABLE')&(rank!=k)]='COLLINEAR_REQUIRED_DESIGN_NO_AUTOMATIC_DELETION'
    with np.errstate(divide='ignore',invalid='ignore'):condition=S[:,0]/S[:,-1]
    status[(status=='ESTIMABLE')&(condition>=1e8)]='ILL_CONDITIONED_REQUIRED_DESIGN'
    if n-k-1<10:status[status=='ESTIMABLE']='FEWER_THAN10_FULL_MODEL_RESIDUAL_DF'
    rx=pd.DataFrame(xx).rank(axis=1).to_numpy();ry=pd.DataFrame(yy).rank(axis=1).to_numpy()
    rx-=np.einsum('bnk,bk->bn',U,np.einsum('bnk,bn->bk',U,rx))
    ry-=np.einsum('bnk,bk->bn',U,np.einsum('bnk,bn->bk',U,ry))
    norm=np.sqrt((rx*rx).sum(axis=1)*(ry*ry).sum(axis=1))
    status[(status=='ESTIMABLE')&(norm<1e-12)]='SCORE_EXPLAINED_BY_REQUIRED_COVARIATES'
    values=np.full(count,np.nan);valid=status=='ESTIMABLE'
    values[valid]=(rx[valid]*ry[valid]).sum(axis=1)/norm[valid]
    return values,status


def independent_leave(x,y,C):
    results=[]
    for i in range(len(x)):
        take=np.arange(len(x))!=i;f,_=independent_fit(x[take],y[take],C[take]);results.append(f['effect'] if f else np.nan)
    v=np.array(results);valid=np.isfinite(v)
    return int(valid.sum()),float(np.min(v[valid])) if valid.any() else np.nan,float(np.max(v[valid])) if valid.any() else np.nan


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,16*1024**2)
    source=root/'immune_RNA_MCP_directional_associations_20261018_v1'
    terminal=json.loads((source/'private_terminal.json').read_text())
    assert terminal['source_sha256']=='b68c7e42f5a3a5cb23dbf27d1fc6cf084a846437d1f057091503e08888caf4f3'
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    preflight=root/'immune_RNA_MCP_disjoint_receiver_preflight_20261018_v1'
    assert sha(preflight/'private_terminal.json')=='469df1863a6770168a7c67fd78078455b586f654ef7be3e350ede9a7f3592f2b'
    for name,item in json.loads((preflight/'private_terminal.json').read_text())['files'].items():assert sha(preflight/name)==item['sha256']
    out=root/'immune_RNA_MCP_directional_associations_audit_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_all_random_stream_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        producer_terminal_sha256=sha(source/'private_terminal.json'),
        scope='All394 main positions,1970 point sensitivities,788 discovery descriptions and every main residual permutation/bootstrap value. Independent pandas ranks/OLS, batched-SVD bootstrap projection and explicit whole-family q formula. No partial effect interpretation.'),indent=2)+'\n')
    registry=pd.read_csv(preflight/'private_all_unique_directional_endpoint_slots.tsv',sep='\t')
    inputs=pd.read_csv(preflight/'private_all_slot_case_model_inputs.tsv.gz',sep='\t')
    main=pd.read_csv(source/'private_all394_main_directional_tests.tsv',sep='\t')
    sensitivity=pd.read_csv(source/'private_all1970_point_sensitivities.tsv',sep='\t')
    discovery=pd.read_csv(source/'private_all788_discovery_source_descriptions.tsv',sep='\t')
    assert len(main)==394 and len(sensitivity)==1970 and len(discovery)==788
    assert main.slot.tolist()==registry.slot.tolist()
    controls=[c+'_'+s for s in ['S','R'] for c in ['fraction','logN','cellQ','geneQ']]+['sender_composition','receiver_composition','receptor_RNA']
    permutations=0;bootstraps=0;max_perm_error=0.;max_boot_error=0.;available=0
    for item in registry.itertuples(index=False):
        all_cases=inputs[inputs.slot.eq(item.slot)];hold=all_cases[all_cases.dataset.eq('GSE305078')];d=hold[hold.primary_eligible]
        r=main[main.slot.eq(item.slot)].iloc[0]
        x=d.ligand_RNA.to_numpy();y=d.receiver_score.to_numpy();C=d[controls].to_numpy(float)
        fitted,status=independent_fit(x,y,C)
        if not item.all3_source_fidelity_pass:fitted=None;status='DISJOINT_RECEIVER_NOT_FAITHFUL'
        assert status==r.status and r.source_units==len(hold) and r.complete_units==len(d) and r.reviewed_patient_units==d.reviewed_patient.sum()
        equal(r.unadjusted_rho,rho(x,y));equal(r.adjusted_rho,fitted['effect'] if fitted else np.nan)
        if fitted is not None:
            z=np.load(source/('private_slot_%04d_random_stream.npz'%item.slot),allow_pickle=False)
            assert z['samples'].tolist()==d['sample'].tolist() and z['control_names'].tolist()==controls
            assert np.allclose(z['x'],x) and np.allclose(z['y'],y) and np.allclose(z['controls'],C)
            equal(r.design_condition,fitted['condition']);assert r.residual_df==fitted['residual_df']
            rng=np.random.default_rng(np.random.SeedSequence([2026101811,int(item.slot),0]))
            indices=np.stack([rng.permutation(len(d)) for _ in range(19999)])
            assert np.array_equal(indices,z['permutation_indices'])
            values=np.empty(19999)
            for start in range(0,19999,1024):
                Y=fitted['fit_y'][:,None]+fitted['ry'][indices[start:start+1024]].T
                residual=Y-fitted['D']@np.linalg.lstsq(fitted['D'],Y,rcond=None)[0]
                values[start:start+residual.shape[1]]=fitted['rx']@residual/(np.linalg.norm(fitted['rx'])*np.linalg.norm(residual,axis=0))
            error=float(np.max(np.abs(values-z['permutation_effects'])));assert error<1e-9;max_perm_error=max(max_perm_error,error)
            equal(r.permutation_p,(1+np.count_nonzero(np.abs(values)>=abs(fitted['effect'])-1e-12))/20000)
            rng=np.random.default_rng(np.random.SeedSequence([2026101811,int(item.slot),1]))
            takes=rng.integers(0,len(d),size=(2000,len(d)),dtype=np.int16);assert np.array_equal(takes,z['bootstrap_indices'])
            boot,reasons=boot_svd(x,y,C,takes)
            assert reasons.tolist()==z['bootstrap_reasons'].tolist()
            assert np.allclose(boot,z['bootstrap_effects'],atol=1e-9,rtol=0,equal_nan=True)
            valid=np.isfinite(boot);assert valid.sum()==r.valid_bootstrap
            if valid.any():max_boot_error=max(max_boot_error,float(np.max(np.abs(boot[valid]-z['bootstrap_effects'][valid]))))
            ci=np.quantile(boot[valid],[.025,.975]) if valid.sum()>=1000 else [np.nan,np.nan]
            equal(r.lower,ci[0]);equal(r.upper,ci[1])
            count,low,high=independent_leave(x,y,C);assert count==r.leave_one_estimable_units;equal(low,r.leave_one_min);equal(high,r.leave_one_max)
            permutations+=len(values);bootstraps+=len(boot);available+=1
        else:
            assert pd.isna(r.permutation_p) and pd.isna(r.lower) and pd.isna(r.upper)
            assert not (source/('private_slot_%04d_random_stream.npz'%item.slot)).exists()
        for s in sensitivity[sensitivity.slot.eq(item.slot)].itertuples(index=False):
            subset=d.copy();names=controls.copy();xx=x
            if s.mode=='original8_same_cases':names=controls[:8]
            elif s.mode=='original10_same_cases':names=controls[:10]
            elif s.mode=='LR_product_proxy':xx=np.log1p(subset.ligand_RNA.to_numpy()*subset.receptor_RNA.to_numpy())
            elif s.mode=='receptor10pct':subset=subset[subset.receptor_support10pct];xx=subset.ligand_RNA.to_numpy()
            elif s.mode=='reviewed_patient_only':subset=subset[subset.reviewed_patient];xx=subset.ligand_RNA.to_numpy()
            else:raise AssertionError('UNDECLARED_SENSITIVITY_MODE')
            f,reason=independent_fit(xx,subset.receiver_score.to_numpy(),subset[names].to_numpy(float))
            assert reason==s.status and len(subset)==s.complete_units;equal(s.point_effect,f['effect'] if f else np.nan)
        for s in discovery[discovery.slot.eq(item.slot)].itertuples(index=False):
            g=all_cases[all_cases.dataset.eq(s.dataset)];subset=g[g.primary_eligible]
            xx=subset.ligand_RNA.to_numpy();yy=subset.receiver_score.to_numpy();cc=subset[controls].to_numpy(float)
            f,status=independent_fit(xx,yy,cc);assert status==s.status and len(subset)==s.complete_units and len(g)==s.source_units
            equal(s.unadjusted_rho,rho(xx,yy));equal(s.adjusted_rho,f['effect'] if f else np.nan)
            count,low,high=independent_leave(xx,yy,cc);assert count==s.leave_one_estimable_units;equal(low,s.leave_one_min);equal(high,s.leave_one_max)
        if int(item.slot)%10==0 or int(item.slot)==393:
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),
                completed_main_positions=int(item.slot)+1,planned_main_positions=394,permutation_values_replayed=permutations,bootstrap_values_replayed=bootstraps,
                partial_scientific_effect_review=False))+'\n')
    p=main.permutation_p.fillna(1).to_numpy(float);assert np.allclose(p,main.family_adjustment_p)
    for name,by in [('BH394',False),('BY394',True)]:
        expected=family_q(p,by=by)
        for i,r in enumerate(main.itertuples(index=False)):equal(getattr(r,name),expected[i] if np.isfinite(r.permutation_p) else np.nan)
    assert available==terminal['estimable_main_slots']
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ACCEPTED_ALL394_DIRECTIONAL_TESTS_SENSITIVITIES_AND_COMPLETE_RANDOM_STREAMS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        main_positions=394,estimable_main_positions=available,sensitivity_positions=1970,discovery_positions=788,
        permutation_values_replayed=permutations,bootstrap_values_replayed=bootstraps,max_permutation_error=max_perm_error,max_bootstrap_error=max_boot_error,
        clinical_or_causal_validation=False),indent=2)+'\n')
    print('ALL_DIRECTIONAL_POSITIONS_AND_COMPLETE_RANDOM_STREAMS_INDEPENDENTLY_ACCEPTED')


if __name__=='__main__':main()
