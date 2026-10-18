"""All20 post-result uncertainty, redundancy and source-background diagnostics.

Original main P/BH/BY, original invalid percentile intervals and all models stay
unchanged. Supplementary intervals and backgrounds are explicitly post-result.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,estimate,adjust
from NSCLC_tensor_followup_math_v1 import seed,full_model_wild_t,studentized_interval


def fixture():
    from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate
    errors=[]
    for scenario in range(3):
        rng=np.random.default_rng(2026101927+scenario);x=np.r_[np.zeros(44),np.ones(43)];rng.shuffle(x);C=rng.normal(size=(87,3));C[:,2]=0;C[:2 if scenario else 10,2]=1;y=.7*x+C@np.array([.3,.2,.5])+rng.normal(size=87)*(np.exp(C[:,0]/2) if scenario==2 else 1)
        m,status=design(x,y,C,['a','b','rare'],'clinical');assert m is not None;g,_=estimate(m);values=full_model_wild_t(g,128,2026101928+scenario)
        audit,status=prepare(x,y,C,['a','b','rare'],'clinical');assert audit is not None;D=audit['D'];q=audit['q'];fitted=q@(q.T@y);h=audit['h'];residual=(y-fitted)/(1-h);r=np.random.default_rng(2026101928+scenario);signs=np.where(r.random((128,len(y)))<.5,-1.,1.);base=evaluate(audit)
        independent=np.array([(evaluate(audit,y=fitted+residual*s)['beta']-base['beta'])/evaluate(audit,y=fitted+residual*s)['se'] for s in signs]);error=float(np.max(np.abs(values-independent)));assert error<1e-8;errors.append(error)
    return dict(scenarios=3,direct_QR_full_refits_per_scenario=128,max_studentized_statistic_error=max(errors),rare_control_cases_in_two_scenarios=2,scope='Numerical identity only, not a promise of finite-sample95percent coverage')


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);out=root/'NSCLC_tensor_complete_family_diagnostics_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'NSCLC_communication_tensor_family_20261019_v1';assert sha(root/'NSCLC_communication_tensor_family_audit_20261019_v2/private_terminal.json')=='e83389e75443b9168cd9e3d56a63980c14d42021ddfe1faedb3d3371c6d37101'
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'math_sha256':sha(Path(__file__).with_name('NSCLC_tensor_followup_math_v1.py')),'resource_guard':guard,'scope':'POST-RESULT diagnostic of all20previously audited slots, not a new primary family. The completed original results and insufficient valid case-bootstrap interval draws were observed; complete original outputs are immutable.','intervals':'Original1900valid gate retained, no retries or lowering. Supplementary pointwise fixed-design full-model HC3-adjusted Rademacher wild-t studentized interval,19999draws, entirely separate from original restricted-null P. Also analytic HC3 Student-t CI and secondary analytic BH/BY20 for transparent calibration comparison. All conditional on fixed discovery factors/observed author-record design; no independent-patient or exact finite-sample coverage claim.','redundancy':'All190factor pairs, cosine of all4axes; complete-link clustering with distance=max(1-mean4cos,.5*(1-min3noncontextcos)),cut.10. Thus every grouped pair mean4>=.9 and minimum LR/sender/receiver cosine>=.8. Pure descriptive dedup, no cluster-level FDR or clinical selection. Coarse roles use original author labels, no new formal states.','backgrounds':'All20models additionally condition on observed mean totalLRscore plus numbercapturedtypes; or sender/receiver factorweighted meanlogRNA/loggenes; or both. Report point and analytic HC3 CI only; no new primary P/selection. Endogenous overall communication background changes the relative-pattern question, not proof of technical confounding.','failure_diagnostics':'All40000original case-bootstrap statuses; rarebinary/dummy levels and resampling missing/singleton contrasts diagnosed in aggregate, no controls removed. Original leave-one-author-record effects/invalid positions summarized.','unmodified':'No new cell selection, source download, factor training, rank/seed choice, main P or annotation changes.'}
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n');save('plan_before_all20_post_result_diagnostics.json',plan);save('private_full_model_wild_numerical_fixture.json',fixture())
    read=lambda p:pd.read_csv(p,sep='\t');results=read(src/'private_all20_heldout_clinical_results.tsv').set_index('model_id');assert len(results)==20;slots=results.index.tolist();axes=np.load(src/'private_tensor_axes.npz');states=axes['states'].tolist();cases=axes['heldout_ids'].tolist();available=axes['availability'][129:];rms=float(axes['global_RMS']);native=root/'NSCLC_full_native_gene_units_20261016_v1';u=read(native/'private_all_native_units_and_quality.tsv');units={(r.sampleID,r.sub_cell_type):int(r.unit_index) for r in u.itertuples()};mmap=np.load(native/'mean_logCP10K.npy',mmap_mode='r');vocab=json.loads((src/'private_selected_LR_vocabulary.json').read_text());N=len(cases);S=len(states);L=len(vocab);burden=np.zeros(N);quality=np.zeros((N,S,2))
    for i,case in enumerate(cases):
        columns=axes['native_gene_columns'];reverse={int(c):j for j,c in enumerate(columns)};M=np.zeros((S,len(columns)))
        for j,state in enumerate(states):
            unit=units.get((case,state))
            if unit is not None:M[j]=mmap[unit,columns];quality[i,j]=u.loc[unit,['mean_logRNA','mean_loggenes']].to_numpy(float)
        lig=np.array([M[:,[reverse[c] for c in d['native_ligand_columns']]].min(axis=1) for d in vocab],np.float32);rec=np.array([M[:,[reverse[c] for c in d['native_receptor_columns']]].min(axis=1) for d in vocab],np.float32)
        y=np.log1p(lig.astype(float)[:,:,None]*rec.astype(float)[:,None,:]).astype(np.float32);a=available[i];burden[i]=np.log1p(float(y[:,a][:,:,a].astype(float).mean())/rms)
    del mmap;capture=available.sum(1).astype(float);np.savez_compressed(out/'private_all87_source_background_inputs.npz',author_record_ids=np.asarray(cases),log1p_mean_observed_LRscore=burden,captured_author_types=capture,unit_meanlogRNA_loggenes=quality)
    factors={};state_meta=read(root/'NSCLC_communication_tensor_inputs_20261019_v1/private_all51_author_type_coverage.tsv').set_index('original_author_type');role=[]
    def coarse(label):
        if label.startswith('CD4T_Treg_'):return 'Treg'
        if label.startswith('CD8T_'):return 'CD8'
        major=state_meta.loc[label,'major'].lower()
        return 'Myeloid' if 'myeloid' in major else 'B' if major.startswith('b') else 'Other_T_NK'
    for slot in slots:
        k=int(slot[2:4]);j=int(slot[-2:])-1;b=dict(np.load(src/('private_selected_K%d_factors_and_projection.npz'%k)));factors[slot]=[b[key][:,j] for key in ['context','LR','sender','receiver']]
        categories={}
        for kind,a in [('sender',factors[slot][2]),('receiver',factors[slot][3])]:
            mass={group:float(sum(v for label,v in zip(states,a) if coarse(label)==group)) for group in sorted(set(map(coarse,states)))};largest=max(mass,key=mass.get);categories[kind+'_coarse_category']=largest if mass[largest]/sum(mass.values())>=.5 else 'Mixed';categories[kind+'_highest_loading_original_author_label']=states[int(np.argmax(a))]
        role.append(dict(model_id=slot,**categories))
    distances=np.zeros((20,20));pairs=[]
    for i in range(20):
        for j in range(i+1,20):
            cos=[float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))) for a,b in zip(factors[slots[i]],factors[slots[j]])];distance=max(0.,1-np.mean(cos),.5*(1-min(cos[1:])));distances[i,j]=distances[j,i]=distance;pairs.append(dict(left=slots[i],right=slots[j],context_cosine=cos[0],LR_cosine=cos[1],sender_cosine=cos[2],receiver_cosine=cos[3],mean4_cosine=float(np.mean(cos)),distance=distance))
    rawgroups=fcluster(linkage(squareform(distances),method='complete'),t=.1,criterion='distance');members=sorted([sorted([slots[i] for i in range(20) if rawgroups[i]==g]) for g in set(rawgroups)]);mapping={slot:'CPG%02d'%(i+1) for i,g in enumerate(members) for slot in g};pd.DataFrame(pairs).to_csv(out/'private_all190_cross_rank_factor_pairs.tsv',sep='\t',index=False);roles=pd.DataFrame(role);roles['descriptive_redundancy_group']=roles.model_id.map(mapping);roles.to_csv(out/'private_all20_original_author_role_and_redundancy.tsv',sep='\t',index=False)
    intervals=[];sensitivity=[];diagnostics=[];global_row=None;binary_names={'male','smoker','LUSC','chemo_binary','targeted_binary'}
    for slot in slots:
        a=dict(np.load(src/('private_'+slot+'_clinical_inputs.npz')));x,y,C,names=a['x'],a['y'],a['C'],a['names'].tolist();assert a['author_record_ids'].tolist()==cases;m,status=design(x,y,C,names,'clinical');assert m is not None;g,effect=estimate(m);np.testing.assert_allclose([effect,np.sqrt(g['variance'])],[results.loc[slot].effect,results.loc[slot].HC3_se],rtol=1e-10,atol=1e-10)
        values=full_model_wild_t(g,19999,seed(slot,'full_model_wild_t'));ci=studentized_interval(g,values);critical=student_t.ppf(.975,m['df']);se=np.sqrt(g['variance']);analytic_P=float(2*student_t.sf(abs(g['t']),m['df']));intervals.append(dict(model_id=slot,source='GSE243013',effect=effect,original_wild_P=results.loc[slot].wild_P,original_BH20=results.loc[slot].BH20,original_BY20=results.loc[slot].BY20,original_percentile_CI_status='UNAVAILABLE_ORIGINAL_VALID_DRAW_GATE_UNMET',supplementary_fixed_design_wild_t_CI_low=ci[0],supplementary_fixed_design_wild_t_CI_high=ci[1],supplementary_analytic_HC3_t_CI_low=effect-critical*se,supplementary_analytic_HC3_t_CI_high=effect+critical*se,secondary_analytic_HC3_t_P=analytic_P,residual_df=m['df']))
        np.savez_compressed(out/('private_'+slot+'_fixed_design_wild_t.npz'),values=values)
        U,V=factors[slot][2:];sw=available*U[None,:];rw=available*V[None,:];sq=np.einsum('ns,nsq->nq',sw,quality)/sw.sum(1)[:,None];rq=np.einsum('ns,nsq->nq',rw,quality)/rw.sum(1)[:,None];q=np.column_stack([sq,rq]);total=np.column_stack([burden,capture]);extra_names=['sender_factor_meanlogRNA','sender_factor_meanloggenes','receiver_factor_meanlogRNA','receiver_factor_meanloggenes']
        for name,extra,en in [('total_LR_and_capture',total,['total_LRscore','captured_type_count']),('endpoint_RNA_quality',q,extra_names),('both',np.column_stack([total,q]),['total_LRscore','captured_type_count']+extra_names)]:
            cc=np.column_stack([C,extra]);mm,why=design(x,y,cc,names+en,'clinical');row=dict(model_id=slot,sensitivity=name,status=why,effect=np.nan,HC3_t_CI_low=np.nan,HC3_t_CI_high=np.nan,nuisance_rank=np.nan)
            if mm is not None:
                gg,ee=estimate(mm);critical=student_t.ppf(.975,mm['df']);se=np.sqrt(gg['variance']);row.update(effect=ee,HC3_t_CI_low=ee-critical*se,HC3_t_CI_high=ee+critical*se,nuisance_rank=mm['nuisance_rank'])
            sensitivity.append(row)
        streams=dict(np.load(src/('private_'+slot+'_all_inference_streams.npz')));keys,counts=np.unique(streams['bootstrap_status'],return_counts=True);loo=streams['leave_one_author_effect'];Q=m['Q'];h0=np.sum(Q*Q,axis=1);minor=[]
        for column,name in enumerate(names):
            if name not in binary_names and not name.startswith('agent:'):continue
            unique,count=np.unique(C[:,column],return_counts=True)
            if len(unique)!=2:continue
            value=unique[np.argmin(count)];flag=C[:,column]==value;absent=0;singleton=0
            for indices,reason in zip(streams['bootstrap_indices'],streams['bootstrap_status']):
                if reason=='ESTIMABLE':continue
                take=indices[flag[indices]]
                absent+=int(len(take)==0);singleton+=int(len(take)==1)
            minor.append(dict(control=name,original_minority_records=int(flag.sum()),failed_draws_with_minority_absent=absent,failed_draws_with_single_minority_record=singleton))
        diagnostics.append(dict(model_id=slot,original_bootstrap_reasons=dict(zip(keys,counts.astype(int).tolist())),original_nuisance_max_leverage=float(h0.max()),nuisance_records_over_half_leverage=int((h0>.5).sum()),leave_one_valid=int(np.isfinite(loo).sum()),leave_one_min=float(np.nanmin(loo)),leave_one_max=float(np.nanmax(loo)),leave_one_same_sign=int((np.sign(loo[np.isfinite(loo)])==np.sign(effect)).sum()),minority_control_diagnostics=minor))
        if global_row is None:
            mm,why=design(x,burden,C[:,:-2],names[:-2],'clinical');assert mm is not None;gg,ee=estimate(mm);critical=student_t.ppf(.975,mm['df']);se=np.sqrt(gg['variance']);global_row=dict(status=why,author_records=N,effect=ee,HC3_t_CI_low=ee-critical*se,HC3_t_CI_high=ee+critical*se,scope='POST_RESULT_DESCRIPTIVE_TOTAL_NORMALIZED_LR_SCORE_READOUT_NOT_NEW_MAIN_ENDPOINT')
        save('private_progress.json',dict(phase='ALL20_POST_RESULT_DIAGNOSTICS',model_id=slot,observed_UTC=datetime.now(timezone.utc).isoformat()))
    iv=pd.DataFrame(intervals);iv['secondary_analytic_HC3_BH20']=adjust(iv.secondary_analytic_HC3_t_P);iv['secondary_analytic_HC3_BY20']=adjust(iv.secondary_analytic_HC3_t_P,by=True);iv.to_csv(out/'private_all20_supplementary_fixed_design_intervals.tsv',sep='\t',index=False);pd.DataFrame(sensitivity).to_csv(out/'private_all60_total_signal_and_endpoint_RNA_sensitivities.tsv',sep='\t',index=False);save('private_all20_original_case_bootstrap_diagnostics.json',diagnostics);save('private_total_observed_LR_score_description.json',global_row)
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    receipt=dict(status='COMPLETE_ALL20_POST_RESULT_INTERVALS_190_REDUNDANCY_PAIRS_AND60_BACKGROUND_DIAGNOSTICS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('NSCLC_tensor_followup_math_v1.py')),original_family_audit_sha256=sha(root/'NSCLC_communication_tensor_family_audit_20261019_v2/private_terminal.json'),files=files,original_main_results_changed=False,original_bootstrap_CI_gate_changed=False,scientific_outputs_public=False)
    save('private_terminal.json',receipt);print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
