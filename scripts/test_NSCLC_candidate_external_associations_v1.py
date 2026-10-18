"""All frozen external ligand-programme tests; no endpoint or precision rescue."""
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

import NSCLC_candidate_external_association_math_v1 as math
from prepare_NSCLC_candidate_external_associations_v1 import OUT as INPUT,QCOLS,SCOPE_SHA
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

OUT='NSCLC_candidate_external_association_results_20261018_v1'
STATE=None


def inputs(root):
    src=root/INPUT;proof=src/'private_independent_external_association_input_audit.json'
    assert sha(proof)=='1719ebcf16d4634bb803620c048e028605e1ab076356a255d694bc56b22dd46a'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(src/'private_terminal.json')
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    scope=root/'NSCLC_candidate_external_association_scope_20261018_v1.json';assert sha(scope)==SCOPE_SHA
    contract=json.loads(scope.read_text());assert (contract['permutations'],contract['bootstrap'],contract['min_informative_units'],contract['min_residual_df'])==(19999,2000,20,8)
    eligibility=pd.read_csv(src/'private_all_planned_model_eligibility.tsv.gz',sep='\t',keep_default_na=False)
    cases=pd.read_csv(src/'private_all_case_measurement_units.tsv.gz',sep='\t',keep_default_na=False)
    for c in QCOLS+['ligand','response','receiver_own_ligand','fraction_whole_sender','fraction_whole_receiver']:cases[c]=pd.to_numeric(cases[c],errors='raise')
    assert cases.case_row.tolist()==list(range(len(cases)))
    frozen=json.loads((src/'private_frozen_families_before_effects.json').read_text())
    active=eligibility[eligibility.test_index.ge(0)]
    assert active.test_index.tolist()==list(range(len(active))) and len(active)==sum(frozen['eligible_family_sizes'].values())
    return src,proof,contract,eligibility,cases,frozen


def one(item):
    item=dict(item);g=STATE['cases'].iloc[json.loads(item['active_case_rows_json'])]
    assert g.record.eq(item['record']).all() and g.context_id.eq(item['context_id']).all()
    if item['family']=='reviewed_patients':assert g.reviewed_patient.all()
    x,y,q,blocks=g.ligand.to_numpy(float),g.response.to_numpy(float),g[QCOLS].to_numpy(float),g.block_key.to_numpy(str)
    p=math.prepare(x,y,q,blocks);assert p['keep'].all() and p['design_rank']==item['nuisance_rank'] and p['residual_df']==item['residual_df']
    rho=math.base.effect(p);draws=STATE['contract']['permutations']
    count=math.base.permutation_count(p,draws,math.seed(item['model'],'permutation'))
    pv=(1+count)/(1+draws);size=STATE['sizes'][item['family']]
    upper=float(beta_dist.ppf(1-.01/size,count+1,draws-count)) if count<draws else 1.
    boot=math.bootstrap(p,2000,math.seed(item['model'],'bootstrap'));valid=np.isfinite(boot)
    ci=np.quantile(boot[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
    loo=math.leave_one(p);finite=np.isfinite(loo)
    row=dict(item,partial_rank_rho=rho,P=pv,permutation_exceedances=count,MC_upper_P=upper,
        CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstraps=int(valid.sum()),
        raw_Spearman=float(spearmanr(x,y).statistic),LOO_valid=int(finite.sum()),LOO_all_valid=bool(finite.all()),
        LOO_same_sign=bool(finite.all() and (np.sign(loo)==np.sign(rho)).all()),
        LOO_min=float(np.min(loo[finite])) if finite.any() else np.nan,LOO_max=float(np.max(loo[finite])) if finite.any() else np.nan)
    for label,qq in [('block_only',np.empty((len(g),0))),('ten_QC_without_parent_fractions',q[:,:10]),
                     ('receiver_own_ligand_adjusted',np.column_stack([q,g.receiver_own_ligand.to_numpy(float)])),
                     ('whole_annotation_fractions_adjusted',np.column_stack([q,g[['fraction_whole_sender','fraction_whole_receiver']].to_numpy(float)]))]:
        value,status=math.secondary(x,y,qq,blocks);row[label+'_rho']=value;row[label+'_status']=status
    return row,boot,loo


def main():
    global STATE
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,64*1024**2)
    src,proof,contract,eligibility,cases,frozen=inputs(root)
    original=json.loads((src/'plan_before_external_association_qualification.json').read_text())
    assert original['kernel_sha256']==sha(math.base.__file__)
    active=eligibility[eligibility.test_index.ge(0)]
    out=root/OUT;assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(math.__file__),base_math_sha256=sha(math.base.__file__),
        scope_sha256=SCOPE_SHA,input_audit_sha256=sha(proof),budget=budget,
        frozen_families=frozen,contract=contract,workers=2,BLAS_threads_per_worker=1,
        selection='All and only independently accepted pre-effect eligible designs. Fixed SHA seeds per model/stream; completion order does not affect results.',
        receiver_direction='Already original-prior-oriented in input response; no second sign multiplication.',
        independent_result_and_stream_replay_still_required=True,new_clinical_outcomes=False)
    (out/'plan_before_any_external_candidate_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    STATE=dict(cases=cases,contract=contract,sizes=frozen['eligible_family_sizes'])
    rows=[];boots=np.full((len(active),2000),np.nan);leaves={}
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('fork')) as pool:
        jobs=[pool.submit(one,row) for row in active.to_dict(orient='records')]
        for job in as_completed(jobs):
            row,boot,loo=job.result();rows.append(row);boots[row['test_index']]=boot;leaves['loo_'+str(row['test_index'])]=loo
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed=len(rows),total=len(active),partial_effects_not_for_release=True))+'\n')
    result=pd.DataFrame(rows).sort_values('test_index') if rows else pd.DataFrame(columns=['test_index','family','P','MC_upper_P'])
    for name in ['BH_q','BY_q','MC_upper_BH_q']:result[name]=np.nan
    for family,size in frozen['eligible_family_sizes'].items():
        take=result.family.eq(family);assert take.sum()==size
        if size:
            result.loc[take,'BH_q']=multipletests(result.loc[take,'P'],method='fdr_bh')[1]
            result.loc[take,'BY_q']=multipletests(result.loc[take,'P'],method='fdr_by')[1]
            result.loc[take,'MC_upper_BH_q']=multipletests(result.loc[take,'MC_upper_P'],method='fdr_bh')[1]
    result.to_csv(out/'private_all_estimable_external_effects_and_sensitivities.tsv.gz',sep='\t',index=False,compression='gzip')
    eligibility.to_csv(out/'private_complete_planned_eligibility.tsv.gz',sep='\t',index=False,compression='gzip')
    np.savez_compressed(out/'private_all_bootstrap_and_LOO.npz',bootstrap=boots,**leaves)
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<64*1024**2
    terminal=dict(status='COMPLETE_ALL_FROZEN_EXTERNAL_CANDIDATE_ASSOCIATIONS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),math_sha256=sha(math.__file__),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        full_result_audit_or_biological_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
