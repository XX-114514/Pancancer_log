"""Whole frozen native ligand-programme family; no effect-based endpoint changes."""
import json
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist,spearmanr
from statsmodels.stats.multitest import multipletests

import NSCLC_native_interface_association_math_v1 as kernel
from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from preflight_NSCLC_continuous_prior_native_association_v1 import make_case_partitions
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

CONTEXT=None


def inputs(root):
    src=root/'NSCLC_continuous_prior_native_association_inputs_20261016_v1'
    proof=src/'private_independent_association_input_audit.json';audit=json.loads(proof.read_text())
    assert audit['status'].startswith('PASS_ALL4800_OUTCOME_FREE_') and audit['producer_terminal_sha256']==sha(src/'private_terminal.json')
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    plan=json.loads((src/'plan_before_any_native_interface_association.json').read_text())
    frozen=json.loads((src/'private_frozen_family_before_effects.json').read_text())
    measurements=root/'NSCLC_continuous_prior_native_interfaces_20261016_v2'
    assert sha(measurements/'private_independent_native_interface_audit.json')==frozen['input_audit_sha256']
    with np.load(measurements/'private_all_native_interface_measurements.npz',allow_pickle=False) as saved:z={k:saved[k] for k in saved.files}
    with np.load(src/'private_all208_covariates_only.npz',allow_pickle=False) as saved:clinical=saved['clinical']
    scope=pd.read_csv(src/'private_all4800_pre_effect_eligibility.tsv.gz',sep='\t',keep_default_na=False)
    active=scope[scope.test_index.ge(0)]
    assert len(scope)==4800 and len(active)==frozen['eligible_primary_records']==audit['eligible_primary_records']
    assert active.test_index.tolist()==list(range(len(active)))
    return src,proof,plan,z,clinical,scope


def context(data,z,clinical):
    parts=make_case_partitions(data,clinical);cases=data['cases'].sampleID.to_numpy(str)
    assert np.array_equal(cases,z['cases']) and np.array_equal(data['train'],z['train'])
    units=data['units'];total=units.groupby('sampleID').n_cells.sum().reindex(cases).to_numpy(float)
    assert (total>0).all();fractions=z['cells']/total[None,:];composition={}
    for vi,view in enumerate(VIEWS):
        table=units[units.sub_cell_type.isin(data['labels'][view])].pivot(index='sampleID',columns='sub_cell_type',values='n_cells')
        counts=table.reindex(index=cases,columns=data['labels'][view]).fillna(0).to_numpy(float)
        assert np.array_equal(counts.sum(axis=1),z['cells'][vi]);composition[view]=counts/z['cells'][vi,:,None]
    return dict(data=data,z=z,parts=parts,fractions=fractions,composition=composition,view_index={v:i for i,v in enumerate(VIEWS)})


def one(item):
    state=CONTEXT;z=state['z'];sv=state['view_index'][item['sender_view']];rv=state['view_index'][item['receiver_view']]
    slot=item['slot_index'];ligand=item['ligand_axis_index'];direction=item['direction'];factor=int(item['factor'].rsplit('F',1)[1])-1
    ix,quality,blocks=state['parts'][(item['sender_view'],item['receiver_view'],'internal_holdout')]
    x=z['ligand_mean'][sv,ix,ligand];y=direction*z['response'][slot,ix]
    p=kernel.prepare(x,y,quality,blocks)
    assert len(p['a'])==item['active_cases'] and p['design_rank']==item['nuisance_rank'] and p['residual_df']==item['residual_df']
    assert p['keep'].all()
    rho=kernel.effect(p);count=kernel.permutation_count(p,199999,kernel.seed(item['record'],'Freedman_Lane'))
    pv=(1+count)/200000;m=state['family_size']
    upper=1. if count==199999 else float(beta_dist.isf(.01/m,count+1,199999-count));upper=min(1.,(1+199999*upper)/200000)
    boot=kernel.bootstrap(p,2000,kernel.seed(item['record'],'bootstrap'));valid=np.isfinite(boot)
    ci=np.quantile(boot[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
    loo=kernel.leave_one(p)
    row=dict(item,partial_rank_rho=rho,P=pv,permutation_exceedances=count,MC_upper_P=upper,CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstraps=int(valid.sum()),
             raw_Spearman=float(spearmanr(x,y).statistic),LOO_valid=int(np.isfinite(loo).sum()),LOO_all_valid=bool(np.isfinite(loo).all()),
             LOO_same_sign=bool(np.isfinite(loo).all() and (np.sign(loo)==np.sign(rho)).all()),LOO_min=float(np.nanmin(loo)) if np.isfinite(loo).any() else np.nan,LOO_max=float(np.nanmax(loo)) if np.isfinite(loo).any() else np.nan)
    for part in ['discovery','all_frozen_cases']:
        ii,q,b=state['parts'][(item['sender_view'],item['receiver_view'],part)]
        value,status=kernel.secondary(z['ligand_mean'][sv,ii,ligand],direction*z['response'][slot,ii],q,b)
        row[part+'_rho']=value;row[part+'_status']=status;row[part+'_cases']=len(ii)
    sender_score=z['original_scores'][sv,ix,factor]
    abundance=state['fractions'][[sv,rv]][:,ix].T
    composition=np.column_stack([state['composition'][item['sender_view']][ix],state['composition'][item['receiver_view']][ix]])
    for label,yy,q in [
        ('original_receiver_programme',direction*z['original_scores'][rv,ix,factor],quality),
        ('sender_programme_adjusted',y,np.column_stack([quality,sender_score])),
        ('native_immune_abundance_adjusted',y,np.column_stack([quality,abundance])),
        ('complete_author_subtype_adjusted',y,np.column_stack([quality,composition]))]:
        value,status=kernel.secondary(x,yy,q,blocks);row[label+'_rho']=value;row[label+'_status']=status
    return row,boot,loo


def main():
    global CONTEXT
    data=load();root=data['root'];budget=resource(root,128*1024**2)
    src,proof,contract,z,clinical,scope=inputs(root)
    checked=root/'NSCLC_native_interface_association_kernel_check_20261018_v1/private_independent_native_interface_kernel_check.json'
    verification=json.loads(checked.read_text());assert verification['status'].startswith('PASS_INDEPENDENT_OLS_') and verification['math_sha256']==sha(kernel.__file__)
    active=scope[scope.test_index.ge(0)];out=root/'NSCLC_native_interface_association_results_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(kernel.__file__),input_audit_sha256=sha(proof),input_plan_sha256=sha(src/'plan_before_any_native_interface_association.json'),kernel_check_sha256=sha(checked),contract=contract,
              all_planned_positions=4800,family_size_before_effects=len(active),budget=budget,workers=4,BLAS_threads_per_worker=1,
              schedule='Four bounded fork workers; each fixed record has independent SHA seed, completion order never changes test_index or family.',
              primary_outcome_is_not_clinical=True,clinical_outcomes_loaded=False,direct_communication_or_causality_validated=False)
    (out/'plan_before_all_native_interface_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    CONTEXT=context(data,z,clinical);CONTEXT['family_size']=len(active)
    rows=[];boots=np.full((len(active),2000),np.nan);leaves=[None]*len(active)
    with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('fork')) as executor:
        futures=[executor.submit(one,item) for item in active.to_dict('records')]
        for future in as_completed(futures):
            row,boot,loo=future.result();i=row['test_index'];rows.append(row);boots[i]=boot;leaves[i]=loo
            if len(rows)%10==0 or len(rows)==len(active):
                (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_primary_tests=len(rows),total_primary_tests=len(active),stage='full_fixed_family_resampling',workers=4))+'\n')
    result=pd.DataFrame(rows).sort_values('test_index');assert result.test_index.tolist()==list(range(len(active)))
    for target,source,method in [('BH_q','P','fdr_bh'),('BY_q','P','fdr_by'),('MC_upper_BH_q','MC_upper_P','fdr_bh')]:result[target]=multipletests(result[source].to_numpy(),method=method)[1]
    result['internal_candidate_under_prespecified_rules']=result.BH_q.lt(.05)&result.MC_upper_BH_q.lt(.05)&result.CI_low.gt(0)&result.LOO_same_sign&result.discovery_rho.gt(0)
    result['status']='ESTIMATED_OBSERVATIONAL_INTERNAL_CASE_LIGAND_PROGRAMME_ASSOCIATION_NOT_DIRECT_COMMUNICATION'
    result.to_csv(out/'private_all_primary_and_sensitivity_effects.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
    offsets=np.r_[0,np.cumsum([len(v) for v in leaves])]
    np.savez_compressed(out/'private_all_bootstrap_and_LOO.npz',records=result.record.to_numpy(str),bootstrap=boots,LOO=np.concatenate(leaves) if leaves else np.empty(0),LOO_offsets=offsets)
    report=dict(status='COMPLETE_FIXED_NATIVE_INTERFACE_ASSOCIATION_FAMILY_INDEPENDENT_RESULT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
                files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},new_direct_communication_clinical_prediction_or_causal_validation=False)
    assert sum(item['bytes'] for item in report['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
