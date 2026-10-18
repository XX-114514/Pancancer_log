"""Independent complete tensor-score, all-start objective, projection and stats replay."""
import json
import os
import sys
from datetime import datetime,timezone
from pathlib import Path


def worker():
    import gc
    import numpy as np
    import pandas as pd
    import torch
    from scipy.optimize import nnls,linear_sum_assignment
    from scipy.stats import mannwhitneyu
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
    from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate,wild,boot
    from NSCLC_immune_RNA_MCP_cohort_math_v1 import seed as old_seed
    from run_NSCLC_communication_tensor_family_v1 import seed
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=Path(os.environ['PANCANCER_TENSOR_AUDIT_OUT']);src=root/'NSCLC_communication_tensor_family_20261019_v1';scope=root/'NSCLC_communication_tensor_inputs_20261019_v1';native=root/'NSCLC_full_native_gene_units_20261016_v1'
    read=lambda p:pd.read_csv(p,sep='\t',keep_default_na=False)
    write=lambda n,x:(out/n).write_text(json.dumps(x,indent=2)+'\n')
    def progress(phase,**kw):write('private_progress.json',dict(phase=phase,observed_UTC=datetime.now(timezone.utc).isoformat(),**kw))
    t=json.loads((src/'private_terminal.json').read_text())
    for name,r in t['files'].items():assert sha(src/name)==r['sha256']
    plan=json.loads((src/'plan_before_real_tensor_or_clinical_effects.json').read_text());assert plan['no_real_tensor_or_effect_at_freeze'] and plan['source_sha256']=='43b659116c33ddf6644fc8ac1459b8bb2ce5f40312ca22245306865e4af6ffdd'
    actualfixture=json.loads((src/'private_actual_float32_and_Gram_projection_fixture.json').read_text());assert actualfixture['status']=='ACCEPTED_BEFORE_REAL_EXPRESSION'
    table=read(scope/'private_all243_blind_record_partitions.tsv').set_index('sampleID');train=table.index[table.partition.eq('DISCOVERY')].tolist();test=table.index[table.partition.eq('HELDOUT')].tolist();ids=train+test;assert table.loc[test].profile_block.nunique()==len(test)==87
    states=read(scope/'private_all51_author_type_coverage.tsv');states=states[states.eligible].original_author_type.tolist();S=len(states);u=read(native/'private_all_native_units_and_quality.tsv');unitmap={(r.sampleID,r.sub_cell_type):int(r.unit_index) for r in u.itertuples()}
    sourcegenes=pd.read_csv(root/'GSE243013_genes.csv.gz',keep_default_na=False).geneSymbol.str.upper().str.strip().tolist();lookup={g:i for i,g in enumerate(sourcegenes)}
    definitions=[d for d in json.loads((scope/'private_all_unique_LR_complex_positions.json').read_text()) if d['metadata_primary_eligible']];axis=read(scope/'private_complete_primary_LR_gene_native_axis.tsv');columns=axis.native_column.to_numpy(int);wanted={g:i for i,g in enumerate(axis.gene)}
    means=np.zeros((len(ids),S,len(columns)),float);detection=np.zeros_like(means);N=np.zeros((len(ids),S),int)
    meanm=np.load(native/'mean_logCP10K.npy',mmap_mode='r');detm=np.load(native/'positive_cell_counts.npy',mmap_mode='r')
    for i,case in enumerate(ids):
        for j,state in enumerate(states):
            unit=unitmap.get((case,state))
            if unit is None:continue
            N[i,j]=int(u.iloc[unit].n_cells);means[i,j]=meanm[unit,columns];detection[i,j]=detm[unit,columns]/N[i,j]
    available=N>=20;del meanm,detm
    Ls=[];Rs=[];eligible=[];support=read(src/'private_all4520_discovery_LR_support.tsv');assert len(support)==len(definitions)==4520
    for j,d in enumerate(definitions):
        li=[wanted[g] for g in d['ligand_components']];ri=[wanted[g] for g in d['receptor_components']];left=means[:,:,li].min(2).astype(np.float32);right=means[:,:,ri].min(2).astype(np.float32)
        ld=detection[:129,:,li].min(2);rd=detection[:129,:,ri].min(2);count=sum(bool(np.any((ld[i]>=.1)&available[i])) and bool(np.any((rd[i]>=.1)&available[i])) for i in range(129));ok=count>=43
        r=support.iloc[j];assert [r.canonical_resource_LR_id,r.discovery_contexts_supported,r.eligible]==[d['canonical_resource_LR_id'],count,ok]
        if ok:Ls.append(left);Rs.append(right);eligible.append(d)
    np.testing.assert_equal(json.loads((src/'private_selected_LR_vocabulary.json').read_text()),eligible)
    L=len(eligible);left=np.stack(Ls,axis=1);right=np.stack(Rs,axis=1);del Ls,Rs,means,detection;gc.collect()
    def score(i):
        a=left[i].astype(float);b=right[i].astype(float);y=np.log1p(a[:,:,None]*b[:,None,:]).astype(np.float32);y[:,~available[i],:]=0;y[:,:,~available[i]]=0;return y
    norm=0.;observed=0
    for i in range(129):y=score(i);norm+=float(np.sum(y.astype(float)**2));observed+=L*int(available[i].sum())**2
    rms=float(np.sqrt(norm/observed));dims=json.loads((src/'private_tensor_input_dimensions.json').read_text());np.testing.assert_allclose(rms,dims['global_discovery_RMS'],rtol=1e-12,atol=1e-12);assert dims['observed_training_entries']==observed
    scale=dims['global_discovery_RMS'];axes=np.load(src/'private_tensor_axes.npz');assert axes['states'].tolist()==states and axes['discovery_ids'].tolist()==train and axes['heldout_ids'].tolist()==test;np.testing.assert_array_equal(axes['availability'],available)
    models=read(src/'private_all15_CP_models.tsv');assert len(models)==15 and set(zip(models['rank'],models.restart))=={(k,j) for k in range(2,7) for j in range(3)}
    parsed=[]
    for r in models.itertuples():
        p=src/r.filename;assert sha(p)==r.sha256;f=dict(np.load(p));assert int(f['rank'])==r.rank and int(f['random_seed'])==2026101920+r.restart
        for key,n in [('context',129),('LR',L),('sender',S),('receiver',S)]:assert f[key].shape==(n,r.rank) and np.isfinite(f[key]).all() and (f[key]>=0).all()
        for key in ['LR','sender','receiver']:np.testing.assert_allclose(np.linalg.norm(f[key],axis=0),1,atol=1e-12,rtol=1e-12)
        assert np.all(np.diff(np.linalg.norm(f['context'],axis=0))<=0);parsed.append(f)
    assert torch.cuda.is_available();torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    tensors=[{k:torch.tensor(f[k],device='cuda:0',dtype=torch.float64) for k in ['LR','sender','receiver']} for f in parsed];prednorm=np.zeros(15);cross=np.zeros(15);ynorm=0.
    for i in range(129):
        y=score(i);y/=scale;yy=torch.tensor(y.astype(float),device='cuda:0');ynorm+=float(np.sum(y.astype(float)**2));a=available[i]
        for j,f in enumerate(parsed):
            B,U,V=f['LR'],f['sender'],f['receiver'];G=(B.T@B)*(U[a].T@U[a])*(V[a].T@V[a]);c=f['context'][i];prednorm[j]+=float(c@G@c)
            z=tensors[j];rhs=torch.einsum('lsr,lk,sk,rk->k',yy,z['LR'],z['sender'],z['receiver']).cpu().numpy();cross[j]+=float(c@rhs)
        if i%10==0:progress('CHECK_ALL15_MODEL_OBSERVED_OBJECTIVES',contexts_done=i+1)
    errs=np.sqrt(np.maximum(0,ynorm-2*cross+prednorm)/ynorm);np.testing.assert_allclose(errs,models.observed_relative_error,atol=1e-8,rtol=1e-8)
    del tensors;torch.cuda.empty_cache();selected={};stability=read(src/'private_all_restart_factor_stability.tsv');project=read(src/'private_all1740_factor_context_projections.tsv');registry=read(src/'private_all20_factor_registry.tsv').set_index('model_id');assert len(project)==1740 and len(registry)==20
    for rank in range(2,7):
        rows=np.flatnonzero(models['rank'].to_numpy()==rank);best=min(rows,key=lambda j:(errs[j],models.iloc[j].restart));f=parsed[best];selected[rank]=f;bundle=dict(np.load(src/('private_selected_K%d_factors_and_projection.npz'%rank)))
        for key in ['context','LR','sender','receiver']:np.testing.assert_array_equal(f[key],bundle[key])
        for j in rows:
            g=parsed[j];similarity=np.zeros((rank,rank))
            for key in ['context','LR','sender','receiver']:
                a,b=f[key],g[key];similarity+=(a.T@b)/np.outer(np.linalg.norm(a,axis=0),np.linalg.norm(b,axis=0))/4
            aa,bb=linear_sum_assignment(-similarity)
            r=stability[(stability['rank']==rank)&(stability.restart==models.iloc[j].restart)].set_index('factor')
            for a,b in zip(aa,bb):assert r.loc[a+1].matched_factor==b+1;np.testing.assert_allclose(r.loc[a+1].mean_four_axis_cosine,similarity[a,b],atol=1e-12,rtol=1e-12)
        B,U,V=[torch.tensor(f[key],device='cuda:0',dtype=torch.float64) for key in ['LR','sender','receiver']];template=torch.einsum('lk,sk,rk->lsrk',B,U,V)
        for i in range(87):
            a=available[129+i];ix=torch.as_tensor(np.flatnonzero(a),device='cuda:0');D=template.index_select(1,ix).index_select(2,ix).reshape(-1,rank);y=score(129+i).astype(float)/scale;v=torch.tensor(y[:,a][:,:,a].ravel(),device='cuda:0',dtype=torch.float64)
            G=(D.T@D).cpu().numpy();rhs=(D.T@v).cpu().numpy();eig=np.linalg.eigvalsh(G);condition=float(np.sqrt(eig[-1]/eig[0])) if eig[0]>0 else float('inf');why='ESTIMABLE';coef=np.full(rank,np.nan);error=np.nan
            if a.sum()<3:why='FEWER_THAN3_CAPTURED_AUTHOR_TYPES'
            elif np.linalg.matrix_rank(G)<rank or condition>=1e8:why='UNIDENTIFIABLE_FIXED_FACTOR_PROJECTION'
            elif float(v@v)<=0:why='NO_OBSERVED_EXPRESSION_NORM'
            else:
                chol=np.linalg.cholesky(G);coef,_=nnls(chol.T,np.linalg.solve(chol,rhs),maxiter=rank*100);residual=v-D@torch.tensor(coef,device='cuda:0');error=float(torch.linalg.vector_norm(residual)/torch.linalg.vector_norm(v))
            assert bundle['heldout_status'][i]==why;np.testing.assert_allclose(coef,bundle['heldout_coefficient'][i],atol=1e-6,rtol=1e-6,equal_nan=True);np.testing.assert_allclose(error,bundle['heldout_error'][i],atol=1e-8,rtol=1e-8,equal_nan=True)
            for j in range(rank):
                slot='TC%02dF%02d'%(rank,j+1);r=registry.loc[slot];positive=f['context'][:,j][f['context'][:,j]>0];median=np.median(positive);act=np.log1p(f['context'][:,j]/median);mu=act.mean();sd=act.std();activity=(np.log1p(coef[j]/median)-mu)/sd
                np.testing.assert_allclose([median,mu,sd],[r.positive_training_coefficient_median,r.training_log_activity_mean,r.training_log_activity_sd],atol=1e-12,rtol=1e-12);np.testing.assert_allclose(activity,bundle['heldout_activity'][i,j],atol=1e-6,rtol=1e-6,equal_nan=True)
                row=project[(project.model_id==slot)&project.sampleID.eq(test[i])];assert len(row)==1 and row.projection_status.iloc[0]==why
                np.testing.assert_allclose(row.activity.iloc[0],activity,atol=1e-6,rtol=1e-6,equal_nan=True)
            if i%10==0:progress('INDEPENDENT_EXPLICIT_TEMPLATE_HELDOUT_NNLS',rank=rank,contexts_done=i+1)
        del template;torch.cuda.empty_cache()
    # Rebuild every clinical input directly, then independent QR/random replay.
    clinical=read(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv').set_index('sampleID').loc[test];x=pd.to_numeric(clinical.response,errors='coerce').to_numpy(float);totals=u.groupby('sampleID').n_cells.sum()
    quality=np.column_stack([np.log1p(totals.reindex(test)),*((u.assign(q=u.n_cells*u[key]).groupby('sampleID').q.sum()/totals).reindex(test).to_numpy() for key in ['mean_logRNA','mean_loggenes'])]);agents=pd.get_dummies(table.loc[test,'agent_identity'],prefix='agent:',drop_first=True,dtype=float);base=np.column_stack([clinical[COMMON].to_numpy(float),agents.to_numpy()]);base_names=COMMON+agents.columns.tolist();fraction=N[129:]/totals.reindex(test).to_numpy()[:,None]
    results=read(src/'private_all20_heldout_clinical_results.tsv').set_index('model_id');sens=read(src/'private_all60_background_sensitivity_effects.tsv').set_index(['model_id','sensitivity']);secondary=read(src/'private_all20_secondary_Wilcoxon_results.tsv').set_index('model_id');pvals=[];wvals=[]
    for rank,f in selected.items():
        bundle=dict(np.load(src/('private_selected_K%d_factors_and_projection.npz'%rank)))
        for j in range(rank):
            slot='TC%02dF%02d'%(rank,j+1);y=bundle['heldout_activity'][:,j];U=f['sender'][:,j];V=f['receiver'][:,j];fra=np.column_stack([fraction@(U/U.sum()),fraction@(V/V.sum())]);C=np.column_stack([base,quality,fra]);names=base_names+['logN_all_immune','all_immune_mean_logRNA','all_immune_mean_loggenes','sender_factor_weighted_fraction','receiver_factor_weighted_fraction'];valid=np.isfinite(np.c_[x,y,C]).all(axis=1);xx=x[valid];yy=y[valid];cc=C[valid];strata=(table.loc[test].cancer_type+'|'+table.loc[test].agent_identity).to_numpy(str)[valid];inp=np.load(src/('private_'+slot+'_clinical_inputs.npz'))
            for key,current in [('x',xx),('y',yy),('C',cc),('strata',strata),('names',np.asarray(names)),('author_record_ids',np.asarray(test)[valid])]:np.testing.assert_array_equal(inp[key],current)
            model,why=prepare(xx,yy,cc,names,'clinical');r=results.loc[slot];n0=int((xx==0).sum());n1=int((xx==1).sum())
            if min(n0,n1)<10:model=None;why='FEWER_THAN10_IN_ORIGINAL_RESPONSE_GROUP'
            assert r.status==why and r.author_records==len(xx) and r.nonMPR==n0 and r.pCR_or_MPR==n1 and r.reviewed_independent_patients==0
            P=1.
            if model is not None:
                observed,null=wild(model,19999,seed(slot,'wild'));streams=np.load(src/('private_'+slot+'_all_inference_streams.npz'));np.testing.assert_allclose(null,streams['wild_t'],atol=1e-7,rtol=1e-7)
                extreme=int((np.abs(null)>=abs(observed['t'])-1e-12).sum());P=(1+extreme)/20000;assert r.wild_extreme==extreme;np.testing.assert_allclose([r.wild_P,r.effect,r.HC3_t,r.HC3_se],[P,observed['effect'],observed['t'],observed['se']],atol=1e-8,rtol=1e-8)
                rng=np.random.default_rng(seed(slot,'bootstrap'));indices=np.empty((2000,len(xx)),np.int16);offset=0
                for level in sorted(set(strata)):
                    ix=np.flatnonzero(strata==level);indices[:,offset:offset+len(ix)]=ix[rng.integers(0,len(ix),size=(2000,len(ix)))];offset+=len(ix)
                np.testing.assert_array_equal(indices,streams['bootstrap_indices']);boots,reasons=boot(xx,yy,cc,names,'clinical',indices,model['nuisance_rank']);np.testing.assert_allclose(boots,streams['bootstrap_effect'],atol=1e-7,rtol=1e-7,equal_nan=True);np.testing.assert_array_equal(reasons,streams['bootstrap_status']);good=np.isfinite(boots);assert r.bootstrap_valid==int(good.sum());ci=np.quantile(boots[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan];np.testing.assert_allclose([r.CI_low,r.CI_high],ci,atol=1e-7,rtol=1e-7,equal_nan=True)
                loo=np.full(len(xx),np.nan)
                for a in range(len(xx)):
                    keep=np.arange(len(xx))!=a;mod,status=prepare(xx[keep],yy[keep],cc[keep],names,'clinical',expected_rank=model['nuisance_rank'])
                    if mod is not None:loo[a]=evaluate(mod)['effect']
                np.testing.assert_allclose(loo,streams['leave_one_author_effect'],atol=1e-7,rtol=1e-7,equal_nan=True)
            pvals.append(P)
            for label,cov,cnames in [('unadjusted',np.empty((len(xx),0)),[]),('clinical_only',base[valid],base_names),('clinical_RNA_no_factor_abundance',np.column_stack([base,quality])[valid],base_names+['logN_all_immune','all_immune_mean_logRNA','all_immune_mean_loggenes'])]:
                m,status=prepare(xx,yy,cov,cnames,'clinical');rr=sens.loc[(slot,label)];assert rr.status==status
                if m is not None:
                    e=evaluate(m);np.testing.assert_allclose([rr.effect,rr.HC3_CI_low,rr.HC3_CI_high],[e['effect'],e['effect']-1.959963984540054*e['se'],e['effect']+1.959963984540054*e['se']],atol=1e-7,rtol=1e-7)
            if min(n0,n1)>=10:
                U,Pw=mannwhitneyu(yy[xx==1],yy[xx==0],alternative='two-sided',method='asymptotic',use_continuity=True);np.testing.assert_allclose([secondary.loc[slot].U,secondary.loc[slot].rank_biserial,secondary.loc[slot].P],[U,2*U/(n0*n1)-1,Pw],atol=1e-12,rtol=1e-12)
            else:Pw=1.
            wvals.append(Pw);progress('ALL20_INDEPENDENT_COMPLETE_INFERENCE_STREAMS',model_id=slot)
    def correction(values,by=False):
        p=np.asarray(values);order=np.argsort(p,kind='stable');q=p[order]*20/np.arange(1,21);q*=sum(1/i for i in range(1,21)) if by else 1.;q=np.minimum.accumulate(q[::-1])[::-1];out=np.empty(20);out[order]=np.minimum(q,1);return out
    for by in [False,True]:np.testing.assert_allclose(results['BY20' if by else 'BH20'],correction(pvals,by),atol=1e-12,rtol=1e-12);np.testing.assert_allclose(secondary['secondary_BY20' if by else 'secondary_BH20'],correction(wvals,by),atol=1e-12,rtol=1e-12)
    receipt=dict(status='ACCEPTED_COMPLETE_LR_SUPPORT_ALL15_OBJECTIVES_ALL20_HELDOUT_PROJECTIONS_INPUTS_AND_INFERENCE_REPLAY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),author_types=S,retained_LRs=L,model_objectives=15,heldout_factor_records=1740,clinical_slots=20,all_fixed_random_streams_replayed=True,scientific_effects_public=False)
    write('private_worker_receipt.json',receipt)


def main():
    import subprocess
    from NSCLC_native_MFA_interpretation_inputs_v1 import sha
    from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
    from immune_RNA_MCP_new_source_io_v1 import wait
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);out=root/'NSCLC_communication_tensor_family_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_independent_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All4520support positions, all15model observed objectives and selected starts, every20x87projection using explicitly constructed masked templates/Gram, all20complete clinical inputs and every19999wild/2000bootstrap/LOO stream. No partial family interpretation.'),indent=2)+'\n')
    wait(root,out,'run_NSCLC_communication_tensor_family_v1',root/'NSCLC_communication_tensor_family_20261019_v1/private_terminal.json',{'COMPLETE_ALL15_CP_MODELS_ALL20_HELDOUT_CLINICAL_SLOTS_FULL_INDEPENDENT_AUDIT_PENDING'})
    python=Path(sys.executable).resolve().parents[2]/'scvi/bin/python';env=os.environ.copy();env.update(PANCANCER_TENSOR_AUDIT_OUT=str(out),PYTHONDONTWRITEBYTECODE='1',NUMBA_CACHE_DIR=str(out/'numba_cache'),MPLCONFIGDIR=str(out/'matplotlib_cache'),MPLBACKEND='Agg')
    with (out/'private_worker_stdout.log').open('xb') as so,(out/'private_worker_stderr.log').open('xb') as se:r=subprocess.run([str(python),'-u',str(Path(__file__).resolve()),'--worker'],env=env,stdout=so,stderr=se,timeout=86400)
    assert r.returncode==0,'COMPLETE_TENSOR_AUDIT_WORKER_FAILED';d=json.loads((out/'private_worker_receipt.json').read_text());d.update(resource_guard=guard,worker_receipt_sha256=sha(out/'private_worker_receipt.json'))
    assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<4*1024**2
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':
    if sys.argv[1:]==['--worker']:worker()
    else:main()
