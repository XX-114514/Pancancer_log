"""Outcome-blind communication CP discovery followed by all20 heldout questions.

No tensor or cell matrix is persisted. All records, factors and effects private.
"""
import hashlib
import json
import os
import sys
import time
from datetime import datetime,timezone
from pathlib import Path

CONTRACT={
    'name':'2026-10-19 communication tensor programmes v1 GSE243013',
    'scores':'log1p(min ligand component mean_logCP10K * min receptor component mean_logCP10K), all native cells per author unit; expression cooccurrence proxy, not full LIANA consensus scores or binding/flux',
    'scope':'Frozen129discovery/87heldout and38eligible author types; full metadata-primary4520LR vocabulary, every position retained in support manifest. Global LRsupport: all ligand components detected>=10percent in at least one available sender AND all receptor components in at least one available receiver within same discovery context, >=ceil129/3 contexts. No heldout or outcome used. All measurable per-context values retained after global filter, including true zero/subthreshold detection.',
    'missing':'Unit absent or fewer20cells: entire corresponding sender/receiver slice masked missing. Every retained LR native component required. Diagonal sender=receiver entries retained as within-population interfaces, not demonstrated single-cell autocrine.',
    'scaling':'Only one global RMS over observed discovery tensor entries; no per-LR/clinical feature scaling. Exact same scalar for heldout. Noncontext factor columns each unit L2; absorb product of their norms and CP weights into context coefficients. Factor numbering decreasing discovery context column L2, deterministic tie order.',
    'fit':'Official Tensorly0.9.0 non_negative_parafac, same core called by cell2cell0.8.0; random init seeds2026101920,2026101921,2026101922 for each rank2..6. Exactly400iterations/tol0, float32 GPU, TF32off, no outcome convergence/seed/rank selection. Retain every15model and choose minimum observed relative reconstruction error per rank. No single selected rank; all20factor clinical slots enter main BH/BY20.',
    'stability':'Hungarian match each same-rank restart using mean cosine similarity of four normalized factor axes; report every matched similarity, median>=.8 descriptive stability criterion. No clinical test filtered by stability.',
    'projection':'Fixed LR/sender/receiver factors. Per heldout context compute observed NNLS through Gram=(B.T@B)*(U_available.T@U_available)*(V_available.T@V_available), rhs exact scored tensor contraction. Cholesky normal equations converted to scipy NNLS, independent fixture equality required. >=3available author types, positive observed score norm, full Gram rank and sqrt(condition)<1e8. Unidentifiable remains missing. No relearning factors or rank on heldout.',
    'clinical':'After all15fits and all20heldout projections saved: original pCR/MPR1 vs nonMPR0; postneoadjuvant surgery association only. Activity=log1p(context coefficient / median positive discovery coefficient), centered/scaled using discovery activity mean/SD only. Continuous activity is model readout; pathology response binary target regressor. Main covariates original8clinical + canonical agent fixed effects + log1p all immune cells + all-cell meanlogRNA/meanloggenes + sender/receiver factor-L1-weighted state fractions of all source immune cells. Complete nuisance column space retained and rank/span reported; no outcome/clinical P variable selection.',
    'clinical_gates':'At least50complete author records, >=10each response, >=10residual df, target identifiable, nuisance and complete design stable/condition<1e8/leverage<.999999. Heldout87nonoutcome profile blocks eachsingleton; no reviewed independent-patient alias assertion. If future non-singleton block discovered, stop before clinical inference rather than sample-level randomization.',
    'statistics':'Reuse previously independently accepted NSCLC cohort HC3 restricted-null Rademacher wild t,19999draws; two-sided add-one P. Fixed2000within cancer x canonical agent resamples, report valid; percentile CI only if>=1900valid, no resampling until desired result. All20main slots BH and BY with unestimable P=1 included. Fixed500seedoffsets defined in code/hash; no seed search. Full random vectors persisted privately for independent replay.',
    'sensitivity':'All20unadjusted group contrasts, original8+agent clinical-only, clinical+RNA/count without factor abundance: point effect and asymptotic HC3 CI only, no extra significance selection. All20unadjusted two-sided Mann-Whitney U (Wilcoxon rank-sum) with tie/asymptotic correction, separate BH/BY20 as secondary only. Leave-one-author-record main effects; no claim independent-patient leaveout.',
    'limits':'Source previously explored; heldout is separation of model learning, not untouched external validation. Immune-only labels are original author subtypes, not new unified formal states. Distinct from 2026 October9CM and 2026-10-18 RNA MCP. Original closed clinical/function families not reopened. No direct interaction, resistance mechanism or baseline prediction claim.',
    'resources':'No downloads. Tensor CPU/GPU RAM only, no full tensor persisted. Max CPUavailable>8*tensorbytes+4GiB, GPUfree>6*tensorbytes+2GiB; otherwise explicit stop/no feature downselection or silent fallback. Output cap32MiB; group reserve500000000000+1073741824 maintained.'}


def seed(slot,stream):
    return int.from_bytes(hashlib.sha256(('2026101915|'+slot+'|'+stream).encode()).digest()[:8],'little')


def worker():
    import gc
    import numpy as np
    import pandas as pd
    import torch
    from scipy.optimize import nnls,linear_sum_assignment
    from scipy.stats import mannwhitneyu
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
    from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,estimate,wild_values,takes,bootstrap,leave_one,adjust
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_TENSOR_OUT'])
    sys.path[:0]=[str(root/('Tensor_cell2cell_GPU_runtime_20261019_'+v)/'site_packages') for v in ['v3','v2']]
    import tensorly as tl
    from tensorly.decomposition import non_negative_parafac
    tl.set_backend('pytorch');assert torch.cuda.is_available();torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    save=lambda n,d:(out/n).write_text(json.dumps(d,indent=2)+'\n')
    progress=lambda **k:save('private_progress.json',dict(observed_UTC=datetime.now(timezone.utc).isoformat(),**k))
    # Qualify the actual float32 solver and the low-memory NNLS shortcut before
    # touching real expression. The previous full gate used float64/direct NNLS.
    rng=np.random.default_rng(2026101916);fixture=rng.uniform(.1,1.,(7,9,4,3)).astype(np.float32);fixture[:,-1]=0
    mask=(rng.random(fixture.shape)>.2).astype(np.float32);mask[:,-1]=1
    def tiny(device,fill):
        x=tl.tensor(np.where(mask,fixture,fill),device=device,dtype=torch.float32);m=tl.tensor(mask,device=device,dtype=torch.float32)
        return tl.to_numpy(tl.cp_to_tensor(non_negative_parafac(x,2,n_iter_max=64,init='random',random_state=2026101917,tol=0,mask=m)))
    gpu=tiny('cuda:0',0);other=tiny('cuda:0',37);cpu=tiny('cpu',0)
    np.testing.assert_array_equal(gpu,other);np.testing.assert_allclose(gpu,cpu,atol=2e-5,rtol=2e-5)
    max_projection_difference=0.
    for k in range(2,7):
        B=rng.uniform(.1,1.,(11,k));U=rng.uniform(.1,1.,(4,k));V=rng.uniform(.1,1.,(4,k))
        for ntypes in [3,4]:
            a=np.arange(4)<ntypes;D=np.einsum('lk,sk,rk->lsrk',B,U[a],V[a]).reshape(-1,k);coef=rng.uniform(.1,1.,k);coef[-1]=0
            yy=np.maximum(0,D@coef+rng.normal(0,.01,len(D)));direct,_=nnls(D,yy);G=(B.T@B)*(U[a].T@U[a])*(V[a].T@V[a]);R=np.linalg.cholesky(G).T;shortcut,_=nnls(R,np.linalg.solve(R.T,D.T@yy))
            np.testing.assert_allclose(shortcut,direct,atol=1e-8,rtol=1e-8);max_projection_difference=max(max_projection_difference,float(np.max(np.abs(shortcut-direct))))
    save('private_actual_float32_and_Gram_projection_fixture.json',dict(status='ACCEPTED_BEFORE_REAL_EXPRESSION',completed_UTC=datetime.now(timezone.utc).isoformat(),float32_placeholder_invariance_exact=True,CPU_GPU_max_error=float(np.max(np.abs(gpu-cpu))),all_ranks2to6_boundary_NNLS_max_error=max_projection_difference,source_sha256=sha(__file__)))
    del gpu,other,cpu;torch.cuda.empty_cache()
    read=lambda p:pd.read_csv(p,sep='\t',keep_default_na=False)
    src=root/'NSCLC_communication_tensor_inputs_20261019_v1';native=root/'NSCLC_full_native_gene_units_20261016_v1'
    blind=read(src/'private_all243_blind_record_partitions.tsv').set_index('sampleID');train=blind.index[blind.partition.eq('DISCOVERY')].tolist();test=blind.index[blind.partition.eq('HELDOUT')].tolist();ids=train+test
    assert len(train)==129 and len(test)==87 and blind.loc[test].profile_block.nunique()==len(test)
    st=read(src/'private_all51_author_type_coverage.tsv');states=st.loc[st.eligible,'original_author_type'].tolist();assert len(states)==38
    axis=read(src/'private_complete_primary_LR_gene_native_axis.tsv');cols=axis.native_column.to_numpy(int);colmap={v:i for i,v in enumerate(cols)}
    rawlr=json.loads((src/'private_all_unique_LR_complex_positions.json').read_text());lr=[d for d in rawlr if d['metadata_primary_eligible']];assert len(lr)==4520
    units=read(native/'private_all_native_units_and_quality.tsv');index={(r.sampleID,r.sub_cell_type):int(r.unit_index) for r in units.itertuples()};N=len(ids);S=len(states)
    lookup=np.array([[index.get((case,state),-1) for state in states] for case in ids]);present=lookup>=0;ncell=np.zeros((N,S),int);ncell[present]=units.n_cells.to_numpy()[lookup[present]];available=ncell>=20
    means=np.zeros((N,S,len(cols)),float);detect=np.zeros_like(means)
    for name,target in [('mean_logCP10K',means),('positive_cell_counts',detect)]:
        mm=np.load(native/(name+'.npy'),mmap_mode='r');target[present]=mm[lookup[present]][:,cols]
        if name=='positive_cell_counts':target[present]/=ncell[present,None]
        del mm
    lig=np.zeros((N,len(lr),S),np.float32);rec=np.zeros_like(lig);support=[]
    for j,d in enumerate(lr):
        L=[colmap[v] for v in d['native_ligand_columns']];R=[colmap[v] for v in d['native_receptor_columns']]
        lig[:,j,:]=means[:,:,L].min(axis=2);rec[:,j,:]=means[:,:,R].min(axis=2)
        left=(detect[:len(train),:,L].min(axis=2)>=.1)&available[:len(train)];right=(detect[:len(train),:,R].min(axis=2)>=.1)&available[:len(train)]
        count=int((left.any(axis=1)&right.any(axis=1)).sum());support.append(dict(canonical_resource_LR_id=d['canonical_resource_LR_id'],discovery_contexts_supported=count,eligible=count>=43))
    pd.DataFrame(support).to_csv(out/'private_all4520_discovery_LR_support.tsv',sep='\t',index=False);keep=np.array([r['eligible'] for r in support]);assert keep.sum()>=6
    chosen=[d for d,k in zip(lr,keep) if k];lig=lig[:,keep];rec=rec[:,keep];L=int(keep.sum());del means,detect;gc.collect()
    contextmask=available[:,:,None]&available[:,None,:];bytes_tensor=len(train)*L*S*S*4
    memory={s.split(':')[0]:s.split(':')[1].strip() for s in Path('/proc/meminfo').read_text().splitlines()};mem=int(memory['MemAvailable'].split()[0])*1024;gpu_free,gpu_total=torch.cuda.mem_get_info()
    assert mem>8*bytes_tensor+4*1024**3 and gpu_free>6*bytes_tensor+2*1024**3,'FULL_SCOPE_RAM_OR_GPU_GUARD'
    # Compute in float64, then store float32. Missing slices are zero placeholders only.
    def score(ci):return np.log1p(lig[ci].astype(float)[:,:,None]*rec[ci].astype(float)[:,None,:]).astype(np.float32)*contextmask[ci][None,:,:]
    Y=np.empty((len(train),L,S,S),np.float32)
    for i in range(len(train)):Y[i]=score(i)
    observed=int(contextmask[:len(train)].sum())*L;scale=float(np.sqrt(np.sum(Y.astype(float)**2)/observed));assert np.isfinite(scale) and scale>0
    Y/=scale
    save('private_tensor_input_dimensions.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),discovery_contexts=len(train),heldout_contexts=len(test),author_types=S,LRs=L,tensor_float32_bytes=bytes_tensor,observed_training_entries=observed,global_discovery_RMS=scale,CPU_available_bytes=mem,GPU_free_bytes=gpu_free,GPU_total_bytes=gpu_total,all4520support_manifest_sha256=sha(out/'private_all4520_discovery_LR_support.tsv'),source_native_replay_sha256=sha(root/'NSCLC_complete_native_LR_statistics_audit_20261019_v1/private_terminal.json')))
    save('private_selected_LR_vocabulary.json',chosen);pd.DataFrame({'sampleID':ids,'partition':['DISCOVERY']*len(train)+['HELDOUT']*len(test),'available_author_types':available.sum(axis=1)}).to_csv(out/'private_all216_context_coverage.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_tensor_axes.npz',states=np.asarray(states),discovery_ids=np.asarray(train),heldout_ids=np.asarray(test),LR_ids=np.asarray([d['canonical_resource_LR_id'] for d in chosen]),native_gene_columns=cols,availability=available,global_RMS=np.array(scale))
    t=tl.tensor(Y,device='cuda:0',dtype=torch.float32);m=tl.tensor(np.broadcast_to(contextmask[:len(train),None,:,:],Y.shape).copy(),device='cuda:0',dtype=torch.float32)
    start=time.monotonic();allmodels=[];selected={};stability=[]
    def normalize(cp):
        weights=tl.to_numpy(cp.weights).astype(float);f=[tl.to_numpy(a).astype(float) for a in cp.factors];C=f[0]*weights
        for j in range(1,4):
            norm=np.linalg.norm(f[j],axis=0);assert (norm>0).all();f[j]/=norm;C*=norm
        order=np.argsort(-np.linalg.norm(C,axis=0),kind='stable');return [C[:,order]]+[a[:,order] for a in f[1:]]
    def error(f):
        sse=0.;norm=0.
        for j in range(0,L,32):
            y=Y[:,j:j+32].astype(float);recon=np.einsum('ak,bk,ck,dk->abcd',f[0],f[1][j:j+32],f[2],f[3],optimize=True);delta=(recon-y)*contextmask[:len(train),None,:,:];sse+=float(np.sum(delta**2));norm+=float(np.sum(y**2))
        return float(np.sqrt(sse/norm))
    for rank in range(2,7):
        runs=[]
        for restart,random_seed in enumerate([2026101920,2026101921,2026101922]):
            progress(phase='FIT_ALL_FIXED_CP_MODELS',rank=rank,restart=restart,elapsed_seconds=time.monotonic()-start)
            cp=non_negative_parafac(t,rank=rank,n_iter_max=400,init='random',random_state=random_seed,tol=0,mask=m,normalize_factors=False);f=normalize(cp);err=error(f)
            name='private_CP_K%d_start%d.npz'%(rank,restart);np.savez_compressed(out/name,context=f[0],LR=f[1],sender=f[2],receiver=f[3],rank=rank,random_seed=random_seed,observed_relative_error=err)
            allmodels.append(dict(rank=rank,restart=restart,random_seed=random_seed,observed_relative_error=err,filename=name,sha256=sha(out/name)));runs.append((err,restart,f));del cp;gc.collect();torch.cuda.empty_cache()
        best=min(runs,key=lambda x:(x[0],x[1]));selected[rank]=best
        for err,restart,f in runs:
            similarity=np.zeros((rank,rank))
            for a,b in zip(best[2],f):similarity+=(a.T@b)/(np.linalg.norm(a,axis=0)[:,None]*np.linalg.norm(b,axis=0)[None,:])/4
            a,b=linear_sum_assignment(-similarity)
            for ai,bi in zip(a,b):stability.append(dict(rank=rank,best_restart=best[1],restart=restart,factor=int(ai+1),matched_factor=int(bi+1),mean_four_axis_cosine=float(similarity[ai,bi])))
        del runs;gc.collect()
    pd.DataFrame(allmodels).to_csv(out/'private_all15_CP_models.tsv',sep='\t',index=False);pd.DataFrame(stability).to_csv(out/'private_all_restart_factor_stability.tsv',sep='\t',index=False)
    del t,m;torch.cuda.empty_cache();del Y;gc.collect();projection=[];factor_registry=[]
    for rank,(besterr,restart,f) in selected.items():
        C,B,U,V=f;rawC=np.full((len(test),rank),np.nan);status=[];err=[];condition=[]
        for i in range(len(test)):
            ci=len(train)+i;avail=available[ci];G=(B.T@B)*(U[avail].T@U[avail])*(V[avail].T@V[avail]);eig=np.linalg.eigvalsh(G);cond=float(np.sqrt(eig[-1]/eig[0])) if eig[0]>0 else float('inf');why='ESTIMABLE';relative=float('nan')
            if avail.sum()<3:why='FEWER_THAN3_CAPTURED_AUTHOR_TYPES'
            elif np.linalg.matrix_rank(G)<rank or cond>=1e8:why='UNIDENTIFIABLE_FIXED_FACTOR_PROJECTION'
            else:
                y=score(ci).astype(float)/scale;norm=float(np.sum(y*y))
                if norm<=0:why='NO_OBSERVED_EXPRESSION_NORM'
                else:
                    rhs=np.einsum('lsr,lk,sk,rk->k',y,B,U,V,optimize=True);R=np.linalg.cholesky(G).T;z=np.linalg.solve(R.T,rhs);solution,_=nnls(R,z,maxiter=100*rank);rawC[i]=solution;sse=max(0.,norm-2*solution@rhs+solution@G@solution);relative=float(np.sqrt(sse/norm))
            status.append(why);err.append(relative);condition.append(cond)
        transformed=np.zeros_like(C);held=np.full_like(rawC,np.nan)
        for j in range(rank):
            slot='TC%02dF%02d'%(rank,j+1);positive=C[:,j][C[:,j]>0];median=float(np.median(positive)) if len(positive) else float('nan');activity=np.log1p(C[:,j]/median);mean=float(activity.mean());sd=float(activity.std())
            if not np.isfinite(sd) or sd<=1e-12:held[:,j]=np.nan;why='CONSTANT_DISCOVERY_ACTIVITY'
            else:held[:,j]=(np.log1p(rawC[:,j]/median)-mean)/sd;why='FROZEN_TRAINING_ACTIVITY_SCALE'
            factor_registry.append(dict(model_id=slot,rank=rank,factor=j+1,best_restart=restart,training_relative_error=besterr,positive_training_coefficient_median=median,training_log_activity_mean=mean,training_log_activity_sd=sd,status=why))
            for i,case in enumerate(test):projection.append(dict(model_id=slot,sampleID=case,rank=rank,factor=j+1,projection_status=status[i],projection_condition=condition[i],relative_observed_projection_error=err[i],coefficient=rawC[i,j],activity=held[i,j]))
        np.savez_compressed(out/('private_selected_K%d_factors_and_projection.npz'%rank),context=C,LR=B,sender=U,receiver=V,heldout_coefficient=rawC,heldout_activity=held,heldout_status=np.asarray(status),heldout_error=np.asarray(err),heldout_condition=np.asarray(condition))
    pd.DataFrame(projection).to_csv(out/'private_all1740_factor_context_projections.tsv',sep='\t',index=False);pd.DataFrame(factor_registry).to_csv(out/'private_all20_factor_registry.tsv',sep='\t',index=False)
    progress(phase='ALL_FITS_AND_HELDOUT_PROJECTIONS_SAVED_BEFORE_CLINICAL',elapsed_seconds=time.monotonic()-start)
    # Only now load the clinical outcome for the complete predeclared family.
    clinical=read(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv').set_index('sampleID').loc[test]
    response=pd.to_numeric(clinical.response,errors='coerce').to_numpy(float);assert blind.loc[test].profile_block.nunique()==len(test)
    totals=units.groupby('sampleID').n_cells.sum();quality=np.column_stack([np.log1p(totals.reindex(test)),*((units.assign(z=units.n_cells*units[field]).groupby('sampleID').z.sum()/totals).reindex(test).to_numpy() for field in ['mean_logRNA','mean_loggenes'])])
    agent=pd.get_dummies(blind.loc[test,'agent_identity'],prefix='agent:',drop_first=True,dtype=float);base=np.column_stack([clinical[COMMON].to_numpy(float),agent.to_numpy()]);base_names=COMMON+agent.columns.tolist();fractions=ncell[len(train):]/totals.reindex(test).to_numpy()[:,None]
    clinical_rows=[];sensitivity=[];wilcoxon=[]
    for rank,(besterr,restart,f) in selected.items():
        bundle=dict(np.load(out/('private_selected_K%d_factors_and_projection.npz'%rank)))
        for j in range(rank):
            slot='TC%02dF%02d'%(rank,j+1);y=bundle['heldout_activity'][:,j];U=f[2][:,j];V=f[3][:,j];abundance=np.column_stack([fractions@(U/U.sum()),fractions@(V/V.sum())]);C=np.column_stack([base,quality,abundance]);names=base_names+['logN_all_immune','all_immune_mean_logRNA','all_immune_mean_loggenes','sender_factor_weighted_fraction','receiver_factor_weighted_fraction'];valid=np.isfinite(np.c_[response,y,C]).all(axis=1);x=response[valid];yy=y[valid];cc=C[valid];strata=(blind.loc[test,'cancer_type']+'|'+blind.loc[test,'agent_identity']).to_numpy(str)[valid]
            progress(phase='ALL20_FIXED_HELDOUT_CLINICAL_QUESTIONS',model_id=slot,elapsed_seconds=time.monotonic()-start)
            np.savez_compressed(out/('private_'+slot+'_clinical_inputs.npz'),x=x,y=yy,C=cc,names=np.asarray(names),strata=strata,author_record_ids=np.asarray(test)[valid])
            model,why=design(x,yy,cc,names,'clinical');group0=int((x==0).sum());group1=int((x==1).sum());row=dict(model_id=slot,author_records=len(x),nonMPR=group0,pCR_or_MPR=group1,reviewed_independent_patients=0,status=why,effect=float('nan'),CI_low=float('nan'),CI_high=float('nan'),wild_P=1.,wild_extreme=0,bootstrap_valid=0)
            if min(group0,group1)<10:model=None;row['status']='FEWER_THAN10_IN_ORIGINAL_RESPONSE_GROUP'
            if model is not None:
                g,effect=estimate(model);null=wild_values(g,19999,seed(slot,'wild'));extreme=int((np.abs(null)>=abs(g['t'])-1e-12).sum());ix=takes(strata,2000,seed(slot,'bootstrap'));boots,reasons=bootstrap(x,yy,cc,names,'clinical',ix,model['nuisance_rank']);loo=leave_one(x,yy,cc,names,'clinical',model['nuisance_rank']);good=np.isfinite(boots);ci=np.quantile(boots[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan]
                row.update(effect=effect,CI_low=float(ci[0]),CI_high=float(ci[1]),wild_P=(1+extreme)/20000,wild_extreme=extreme,bootstrap_valid=int(good.sum()),HC3_t=g['t'],HC3_se=float(np.sqrt(g['variance'])),nuisance_rank=model['nuisance_rank'],residual_df=model['df'],design_condition=model['condition'],max_leverage=model['max_leverage'])
                np.savez_compressed(out/('private_'+slot+'_all_inference_streams.npz'),wild_t=null,bootstrap_indices=ix,bootstrap_effect=boots,bootstrap_status=reasons,leave_one_author_effect=loo)
            clinical_rows.append(row)
            for label,cov,covnames in [('unadjusted',np.empty((len(x),0)),[]),('clinical_only',base[valid],base_names),('clinical_RNA_no_factor_abundance',np.column_stack([base,quality])[valid],base_names+['logN_all_immune','all_immune_mean_logRNA','all_immune_mean_loggenes'])]:
                mod,why=design(x,yy,cov,covnames,'clinical');rr=dict(model_id=slot,sensitivity=label,status=why,effect=np.nan,HC3_CI_low=np.nan,HC3_CI_high=np.nan)
                if mod is not None:
                    fit,effect=estimate(mod);se=np.sqrt(fit['variance']);rr.update(effect=effect,HC3_CI_low=effect-1.959963984540054*se,HC3_CI_high=effect+1.959963984540054*se)
                sensitivity.append(rr)
            if min(group0,group1)>=10:
                Ustat,p=mannwhitneyu(yy[x==1],yy[x==0],alternative='two-sided',method='asymptotic',use_continuity=True);w=dict(model_id=slot,status='DESCRIPTIVE_SECONDARY_UNADJUSTED_RANK_SUM',U=float(Ustat),rank_biserial=float(2*Ustat/(group0*group1)-1),P=float(p))
            else:w=dict(model_id=slot,status='INSUFFICIENT_RESPONSE_GROUP_COVERAGE',U=np.nan,rank_biserial=np.nan,P=1.)
            wilcoxon.append(w)
    results=pd.DataFrame(clinical_rows);assert len(results)==20;results['BH20']=adjust(results.wild_P);results['BY20']=adjust(results.wild_P,by=True);results.to_csv(out/'private_all20_heldout_clinical_results.tsv',sep='\t',index=False)
    pd.DataFrame(sensitivity).to_csv(out/'private_all60_background_sensitivity_effects.tsv',sep='\t',index=False);w=pd.DataFrame(wilcoxon);w['secondary_BH20']=adjust(w.P);w['secondary_BY20']=adjust(w.P,by=True);w.to_csv(out/'private_all20_secondary_Wilcoxon_results.tsv',sep='\t',index=False)
    receipt=dict(status='COMPLETE_ALL15_CP_MODELS_ALL20_HELDOUT_CLINICAL_SLOTS_FULL_INDEPENDENT_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),elapsed_seconds=time.monotonic()-start,model_count=15,clinical_slots=20,source_author_types=S,retained_LR_count=L,discovery_records=len(train),heldout_records=len(test),files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()},full_tensor_persisted=False)
    save('private_worker_receipt.json',receipt)


def main():
    import subprocess
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
    from immune_RNA_MCP_new_source_io_v1 import wait
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,32*1024**2);out=root/'NSCLC_communication_tensor_family_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    assert sha(root/'NSCLC_communication_tensor_inputs_audit_20261019_v2/private_terminal.json')=='dfbc55bd2f2fc408afb8780c96f583d24754f5826665be32d5edef9134c1ab51'
    assert sha(root/'Tensor_cell2cell_synthetic_gate_20261019_v1/private_terminal.json')=='e1743140d72f35dcb45791586b40976b40925fcc78f2dabffea2087cacc63a7e'
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,contract=CONTRACT,reused_inference_math_sha256=sha(Path(__file__).parent/'NSCLC_immune_RNA_MCP_cohort_math_v1.py'),no_real_tensor_or_effect_at_freeze=True)
    (out/'plan_before_real_tensor_or_clinical_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    accepted=wait(root,out,'audit_NSCLC_complete_native_LR_statistics_v1',root/'NSCLC_complete_native_LR_statistics_audit_20261019_v1/private_terminal.json',{'ACCEPTED_ALL1840_NATIVE_LR_GENES_ALL11701_UNITS_FULL_RAW_COUNT_DETECTION_LOGCP10K_REPLAY'})
    assert accepted['native_LR_genes']==1840 and accepted['author_units']==11701
    python=Path(sys.executable).resolve().parents[2]/'scvi/bin/python';env=os.environ.copy();env.update(PANCANCER_TENSOR_OUT=str(out),PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(out/'numba_cache'),MPLCONFIGDIR=str(out/'matplotlib_cache'),MPLBACKEND='Agg')
    with (out/'private_worker_stdout.log').open('xb') as so,(out/'private_worker_stderr.log').open('xb') as se:r=subprocess.run([str(python),'-u',str(Path(__file__).resolve()),'--worker'],env=env,stdout=so,stderr=se,timeout=86400)
    assert r.returncode==0,'REAL_TENSOR_FAMILY_WORKER_FAILED';receipt=json.loads((out/'private_worker_receipt.json').read_text());assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<32*1024**2
    receipt.update(native_LR_audit_terminal_sha256=sha(root/'NSCLC_complete_native_LR_statistics_audit_20261019_v1/private_terminal.json'),worker_receipt_sha256=sha(out/'private_worker_receipt.json'))
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:worker()
    else:main()
