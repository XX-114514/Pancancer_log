"""Independent pandas-rank replay of all holdout permutations, intervals and FDR."""
import hashlib,itertools,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import beta as beta_dist


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rngseed(identifier,stream):return int.from_bytes(hashlib.sha256((identifier+'|'+stream+'|NSCLC_MFA_v1').encode()).digest()[:8],'little')


def centred(x,blocks):
    vector=np.asarray(x).ndim==1;array=np.asarray(x).reshape(1,-1) if vector else np.asarray(x);ranked=pd.DataFrame(array).rank(axis=1,method='average').to_numpy()
    for level in sorted(set(blocks)):
        ids=np.where(blocks==level)[0];ranked[:,ids]-=ranked[:,ids].mean(axis=1)[:,None]
    return ranked[0] if vector else ranked


def effect(a,b,blocks):
    active=sum(int(np.sum(blocks==level)) for level in set(blocks) if np.sum(blocks==level)>=2)
    if active<30:return np.nan,None,None,None,active
    u=centred(a,blocks);v=centred(b,blocks);denom=np.linalg.norm(u)*np.linalg.norm(v)
    if not np.isfinite(denom) or denom<1e-12:return np.nan,None,None,None,active
    return float(u@v/denom),u,v,float(denom),active


def bh(p,by=False):
    p=np.asarray(p,float);n=len(p);order=np.argsort(p,kind='stable');sorted_q=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
    if by:sorted_q*=sum(1/k for k in range(1,n+1))
    result=np.empty(n);result[order]=np.minimum(1.,sorted_q);return result


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=root/'NSCLC_native_MFA_holdout_20261016_v2';receipt=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_native_MFA_holdout_v2/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (source/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(receipt['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(3)
    terminal=json.loads((source/'private_terminal.json').read_text());assert terminal['code_sha256']==receipt['source_sha256']
    for n,h in terminal['files'].items():assert sha(source/n)==h
    plan=json.loads((source/'plan_before_holdout_associations.json').read_text());assert plan['contract']['family']==60 and plan['contract']['permutations']==19999 and plan['contract']['bootstrap']==2000
    model=root/'NSCLC_native_MFA_discovery_20261016_v3';s=model/'seed_20261017/private_frozen_parameters_and_own_view_scores.npz';assert sha(s)==plan['primary_score_archive_sha256'];archive=np.load(s,allow_pickle=False);axis=pd.read_csv(root/'NSCLC_native_MFA_inputs_20261016_v1/private_case_axis.tsv',sep='\t',dtype={'sampleID':str});assert archive['cases'].astype(str).tolist()==axis.sampleID.tolist();keep=~archive['train'];blocks=axis.loc[keep,'holdout_block'].to_numpy(str);assert len(blocks)==63;views=archive['views'].astype(str).tolist();values=[archive['own_view_scores_'+str(m)][keep] for m in range(4)]
    table=pd.read_csv(source/'private_all60_internal_holdout_coordination.tsv',sep='\t');saved=np.load(source/'private_all60_bootstrap_and_LOO.npz',allow_pickle=False);boots=saved['bootstrap_rho'];leave_saved=saved['LOO_rho'];assert len(table)==60 and table.slot.is_unique and saved['slots'].astype(str).tolist()==table.slot.tolist() and boots.shape==(60,2000) and leave_saved.shape==(60,63);counter=0;total_perm=0;total_boot=0;maximum=0.
    for k in range(10):
        for m,n in itertools.combinations(range(4),2):
            row=table.iloc[counter];factor='NSCLC_MFA_20261016_F'+str(k+1).zfill(2);identifier=factor+'|'+views[m]+'|'+views[n];assert row.slot==identifier and row.factor==factor and row.view_A==views[m] and row.view_B==views[n] and row.author_cases==63;a=values[m][:,k];b=values[n][:,k];rho,u,v,denom,active=effect(a,b,blocks)
            if not np.isfinite(rho):
                assert not row.status.startswith('ESTIMATED_') and np.isnan(row.P) and np.isnan(boots[counter]).all() and np.isnan(leave_saved[counter]).all();counter+=1;continue
            assert row.status.startswith('ESTIMATED_') and row.active_block_cases==active and np.isclose(row.block_rank_rho,rho,atol=1e-12);raw=pd.Series(a).rank().corr(pd.Series(b).rank());assert np.isclose(row.raw_Spearman,raw,atol=1e-12);rng=np.random.default_rng(rngseed(identifier,'permutation'));count=0
            for start in range(0,19999,256):
                size=min(256,19999-start);take=np.tile(np.arange(63),(size,1))
                for level in sorted(set(blocks)):
                    ids=np.where(blocks==level)[0];keys=rng.random((size,len(ids)));take[:,ids]=ids[np.argsort(keys,axis=1)]
                statistic=np.einsum('ij,j->i',v[take],u)/denom;count+=int((abs(statistic)>=abs(rho)-1e-12).sum())
            assert row.permutation_exceedances==count and np.isclose(row.P,(count+1)/20000,atol=1e-14);upper=1. if count==19999 else beta_dist.ppf(1-.01/60,count+1,19999-count);upper=min(1.,(1+19999*upper)/20000);assert np.isclose(row.MC_upper_P,upper,rtol=1e-10,atol=1e-12);total_perm+=19999
            rng=np.random.default_rng(rngseed(identifier,'bootstrap'));take=np.empty((2000,63),int);labels=[];offset=0
            for level in sorted(set(blocks)):
                ids=np.where(blocks==level)[0];take[:,offset:offset+len(ids)]=ids[rng.integers(0,len(ids),(2000,len(ids)))];offset+=len(ids);labels.extend([level]*len(ids))
            u2=centred(a[take],np.array(labels));v2=centred(b[take],np.array(labels));den=np.linalg.norm(u2,axis=1)*np.linalg.norm(v2,axis=1);expected=np.full(2000,np.nan);good=np.isfinite(den)&(den>1e-12);expected[good]=np.einsum('ij,ij->i',u2[good],v2[good])/den[good];assert np.allclose(expected,boots[counter],rtol=1e-10,atol=1e-12,equal_nan=True);maximum=max(maximum,float(np.nanmax(abs(expected-boots[counter]))));assert row.valid_bootstraps==good.sum();lo,hi=np.quantile(expected[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan];assert np.isclose(row.CI_low,lo,atol=1e-12,equal_nan=True) and np.isclose(row.CI_high,hi,atol=1e-12,equal_nan=True);total_boot+=2000
            loo=np.array([effect(a[np.arange(63)!=j],b[np.arange(63)!=j],blocks[np.arange(63)!=j])[0] for j in range(63)]);valid=np.isfinite(loo);assert np.allclose(loo,leave_saved[counter],rtol=1e-10,atol=1e-12,equal_nan=True) and row.LOO_valid==valid.sum() and row.LOO_same_sign==bool(valid.all() and (np.sign(loo)==np.sign(rho)).all());assert np.isclose(row.LOO_min,np.nanmin(loo),atol=1e-12,equal_nan=True) and np.isclose(row.LOO_max,np.nanmax(loo),atol=1e-12,equal_nan=True);counter+=1
    valid=table.status.str.startswith('ESTIMATED_').to_numpy();p=table.P.fillna(1).to_numpy();up=table.MC_upper_P.fillna(1).to_numpy()
    for col,wanted in [('BH_q',bh(p)),('BY_q',bh(p,True)),('MC_upper_BH_q',bh(up))]:wanted[~valid]=np.nan;assert np.allclose(table[col],wanted,rtol=1e-10,atol=1e-12,equal_nan=True)
    proof=dict(status='PASS_ALL60_NATIVE_MFA_HOLDOUT_SLOTS_FULL_PERMUTATION_BOOTSTRAP_LOO_AND_COMPLETE_FAMILY_CORRECTIONS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),all_slots_checked=counter,permutations_replayed=total_perm,bootstrap_draws_replayed=total_boot,maximum_bootstrap_rho_difference=maximum,reviewed_patient_or_direct_LR_validation=False);p=source/'private_independent_MFA_holdout_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
