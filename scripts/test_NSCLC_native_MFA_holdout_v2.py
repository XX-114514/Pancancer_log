"""All60 frozen internal-holdout view-pair slots, no endpoint selection."""
import itertools,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta as beta_dist,spearmanr
from statsmodels.stats.multitest import multipletests
from NSCLC_native_MFA_math_v1 import sha
from NSCLC_native_MFA_holdout_math_v1 import point,permutation_count,bootstrap,seed
from prepare_NSCLC_native_MFA_inputs_v1 import quota_free


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=root/'NSCLC_native_MFA_discovery_20261016_v3';p=source/'private_independent_MFA_discovery_audit.json';receipt=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_native_MFA_discovery_v4/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not p.exists():
        stat=Path('/proc')/str(receipt['pid'])/'stat';assert stat.exists();a=stat.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(receipt['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(3)
    proof=json.loads(p.read_text());assert proof['status'].startswith('PASS_ALL_THREE_') and not proof['primary_maximum_iteration_reached'];assert proof['producer_terminal_sha256']==sha(source/'private_terminal.json')
    inputs=root/'NSCLC_native_MFA_inputs_20261016_v1';frozen=json.loads((inputs/'plan_before_native_expression_and_MFA_fit.json').read_text());contract=frozen['holdout_before_any_fit'];assert contract['family']==60 and contract['permutations']==19999 and contract['bootstrap']==2000
    scores=np.load(source/'seed_20261017/private_frozen_parameters_and_own_view_scores.npz',allow_pickle=False);axis=pd.read_csv(inputs/'private_case_axis.tsv',sep='\t',keep_default_na=False,dtype={'sampleID':str});assert scores['cases'].astype(str).tolist()==axis.sampleID.tolist();train=scores['train'];assert np.array_equal(train,axis.partition.eq('discovery')) and (~train).sum()==63;blocks=axis.loc[~train,'holdout_block'].to_numpy(str);views=scores['views'].astype(str).tolist();all_scores=[scores['own_view_scores_'+str(m)][~train] for m in range(4)]
    out=root/'NSCLC_native_MFA_holdout_20261016_v2';assert not out.exists();free=quota_free(root);out.mkdir(mode=0o700);plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),statistical_helper_sha256=sha(Path(__file__).with_name('NSCLC_native_MFA_holdout_math_v1.py')),discovery_audit_sha256=sha(p),primary_score_archive_sha256=sha(source/'seed_20261017/private_frozen_parameters_and_own_view_scores.npz'),frozen_pre_fit_contract_sha256=sha(inputs/'plan_before_native_expression_and_MFA_fit.json'),contract=contract,author_cases_not_reviewed_independent_patients=True,quota_free_bytes=free,output_cap_bytes=268435456,clinical_outcomes_used=False);(out/'plan_before_holdout_associations.json').write_text(json.dumps(plan,indent=2)+'\n');rows=[];boot_rows=[];loo_rows=[]
    for k in range(10):
        for m,n in itertools.combinations(range(4),2):
            factor='NSCLC_MFA_20261016_F'+str(k+1).zfill(2);identifier=factor+'|'+views[m]+'|'+views[n];a=all_scores[m][:,k];b=all_scores[n][:,k];B=19999;boots=np.full(2000,np.nan);leave=np.full(63,np.nan);row=dict(slot=identifier,factor=factor,view_A=views[m],view_B=views[n],author_cases=63,active_block_cases=0,status='UNESTIMABLE',raw_Spearman=np.nan,block_rank_rho=np.nan,P=np.nan,MC_upper_P=np.nan,permutation_exceedances=-1,valid_bootstraps=0,CI_low=np.nan,CI_high=np.nan,LOO_valid=0,LOO_same_sign=False,LOO_min=np.nan,LOO_max=np.nan)
            try:rho,u,v,denom,active=point(a,b,blocks)
            except ValueError as error:row['status']='UNESTIMABLE_'+str(error)
            else:
                count=permutation_count(u,v,denom,rho,blocks,B,seed(identifier,'permutation'));pv=(count+1)/(B+1);upper=1. if count==B else beta_dist.isf(.01/60,count+1,B-count);upper=min(1.,(1+B*upper)/(B+1));boots=bootstrap(a,b,blocks,2000,seed(identifier,'bootstrap'));good=np.isfinite(boots);lo,hi=np.quantile(boots[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan]
                for j in range(63):
                    keep=np.arange(63)!=j
                    try:leave[j]=point(a[keep],b[keep],blocks[keep])[0]
                    except ValueError:pass
                valid=np.isfinite(leave);row.update(status='ESTIMATED_EXPLORATORY_INTERNAL_CASE_LEVEL_COORDINATION',active_block_cases=int(active),raw_Spearman=float(spearmanr(a,b).statistic),block_rank_rho=rho,P=pv,MC_upper_P=float(upper),permutation_exceedances=count,valid_bootstraps=int(good.sum()),CI_low=float(lo),CI_high=float(hi),LOO_valid=int(valid.sum()),LOO_same_sign=bool(valid.all() and (np.sign(leave)==np.sign(rho)).all()),LOO_min=float(np.nanmin(leave)) if valid.any() else np.nan,LOO_max=float(np.nanmax(leave)) if valid.any() else np.nan)
            rows.append(row);boot_rows.append(boots);loo_rows.append(leave)
    table=pd.DataFrame(rows);assert len(table)==60 and table.slot.is_unique;ready=table.status.str.startswith('ESTIMATED_').to_numpy()
    for field,method,pcol in [('BH_q','fdr_bh','P'),('BY_q','fdr_by','P'),('MC_upper_BH_q','fdr_bh','MC_upper_P')]:table[field]=multipletests(table[pcol].fillna(1),method=method)[1];table.loc[~ready,field]=np.nan
    table.to_csv(out/'private_all60_internal_holdout_coordination.tsv',sep='\t',index=False,na_rep='NA');np.savez_compressed(out/'private_all60_bootstrap_and_LOO.npz',bootstrap_rho=np.array(boot_rows),LOO_rho=np.array(loo_rows),slots=table.slot.to_numpy(str));assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<plan['output_cap_bytes'];terminal=dict(status='COMPLETE_ALL60_NATIVE_MFA_INTERNAL_HOLDOUT_COORDINATION_BOOTSTRAP_LOO_AND_FULL_FAMILY_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},independent_patient_or_clinical_or_direct_communication_validation=False);(out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
