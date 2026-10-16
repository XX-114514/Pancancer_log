"""Independent direct pooling, OLS/QR, ridge and rank replay of complete scope."""
import itertools
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from NSCLC_native_MFA_interpretation_inputs_v1 import load, quota_gate, sha, VIEWS, PANELS


def ranks(x):
    return pd.Series(x).rank(method='average').to_numpy()


def corr(a,b):
    da=a-a.mean();db=b-b.mean();den=np.linalg.norm(da)*np.linalg.norm(db)
    return float(da@db/den) if den>1e-12 else np.nan


def point(a,b,blocks):
    good=np.isfinite(a)&np.isfinite(b);a=a[good];b=b[good];blocks=blocks[good]
    active=sum(int((blocks==s).sum()) for s in set(blocks) if (blocks==s).sum()>=2)
    if active<30:return np.nan
    ra=ranks(a);rb=ranks(b)
    for s in set(blocks):
        ix=blocks==s;ra[ix]-=ra[ix].mean();rb[ix]-=rb[ix].mean()
    den=np.linalg.norm(ra)*np.linalg.norm(rb)
    return float(ra@rb/den) if den>1e-12 else np.nan


def partial(a,b,sex,blocks):
    good=np.isfinite(a)&np.isfinite(b)&np.isfinite(sex)
    a=a[good];b=b[good];sex=sex[good];blocks=blocks[good]
    design=np.column_stack([*(blocks==s for s in sorted(set(blocks))),sex]).astype(float)
    rank=np.linalg.matrix_rank(design)
    if len(a)<30 or len(a)-rank<10 or rank<design.shape[1]:return np.nan
    q,_=np.linalg.qr(design,mode='reduced');values=np.column_stack([ranks(a),ranks(b)])
    residual=values-q@(q.T@values);den=np.linalg.norm(residual[:,0])*np.linalg.norm(residual[:,1])
    return float(residual[:,0]@residual[:,1]/den) if den>1e-12 else np.nan


def main():
    data=load();root=data['root'];out=root/'NSCLC_native_MFA_interpretation_20261016_v1'
    terminal=out/'private_terminal.json'
    for _ in range(240):
        if terminal.exists():break
        job=root/'HNSCC_receiver_followup_job_control_20261016_v2/measure_NSCLC_native_MFA_interpretation_v1/launch_receipt.json'
        if job.exists():
            receipt=json.loads(job.read_text());proc=Path('/proc')/str(receipt['pid'])/'stat'
            assert proc.exists(),'Producer is absent before terminal'
            fields=proc.read_text().rsplit(')',1)[1].split()
            assert fields[0]!='Z' and fields[19]==receipt['proc_start_ticks'],'Producer stopped before terminal'
        time.sleep(2)
    assert terminal.exists(),'No complete producer terminal'
    proof=out/'private_independent_interpretation_audit.json';assert not proof.exists();quota_gate(root)
    report=json.loads(terminal.read_text())
    for name,record in report['files'].items():assert sha(out/name)==record['sha256']
    plan=json.loads((out/'plan_before_any_new_interpretation_effect.json').read_text())
    assert plan['new_P_or_FDR'] is False and plan['clinical_or_model_refit'] is False
    assert plan['helper_sha256']==sha(Path(__file__).with_name('NSCLC_native_MFA_interpretation_inputs_v1.py'))
    tr=data['train'];ho=~tr;n=len(tr);blocks=data['cases'].holdout_block.to_numpy(str)[ho]
    with np.load(out/'private_all_function_readouts_and_QC_residuals.npz',allow_pickle=False) as z:readouts={k:z[k] for k in z.files}
    with np.load(out/'private_all_mixture_models_and_scores.npz',allow_pickle=False) as z:models={k:z[k] for k in z.files}
    funcs=pd.read_csv(out/'private_all2000_unused_gene_function_correlations.tsv',sep='\t')
    weights=pd.read_csv(out/'private_all2320_weight_function_and_panel_annotations.tsv',sep='\t')
    coverage=pd.read_csv(out/'private_all200_RNA_function_coverage.tsv',sep='\t')
    mixtures=pd.read_csv(out/'private_all40_subtype_mixture_predictions.tsv',sep='\t')
    pairs=pd.read_csv(out/'private_all60_composition_and_sex_sensitivity.tsv',sep='\t')
    covariates=pd.read_csv(out/'private_all_covariate_diagnostics.tsv',sep='\t')
    subtypes=pd.read_csv(out/'private_all_author_subtype_diagnostics.tsv',sep='\t')
    variance=pd.read_csv(out/'private_all40_model_variance_descriptions.tsv',sep='\t')
    assert len(funcs)==2000 and len(weights)==2320 and len(coverage)==200 and len(mixtures)==40 and len(pairs)==60 and len(variance)==40
    assert not funcs.duplicated(['view','factor','reference']).any()
    errors={};checks={}
    def close(name,a,b,atol=1e-8,rtol=1e-8):
        a=np.asarray(a,float);b=np.asarray(b,float)
        assert a.shape==b.shape,(name,a.shape,b.shape)
        assert np.allclose(a,b,atol=atol,rtol=rtol,equal_nan=True),name
        good=np.isfinite(a)&np.isfinite(b)
        error=float(np.max(np.abs(a[good]-b[good]))) if good.any() else 0.
        errors[name]=max(errors.get(name,0.),error);checks[name]=checks.get(name,0)+a.size
    u=data['units'];source=np.load(data['native']/'mean_logCP10K.npy',mmap_mode='r')
    gene=data['genes'];remain=data['remaining'];refs=data['refs'];names=sorted(refs)
    sex=data['demographic'].female.to_numpy()[ho]
    independent_residuals=[];independent_scores=[]
    # Independent pooling uses direct unit additions, not production sparse matrix.
    for vi,view in enumerate(VIEWS):
        p=data['prepared'][vi];scores=data['model']['own_view_scores_'+str(vi)]
        independent_scores.append(scores[ho]);labels=data['labels'][view]
        pooled=np.zeros((n,len(gene)),float);count=np.zeros((n,len(labels)),int)
        for ci,case in enumerate(data['cases'].sampleID):
            selected_units=u[(u.sampleID==case)&u.sub_cell_type.isin(labels)]
            total=int(selected_units.n_cells.sum());assert total>=20
            for item in selected_units.itertuples():
                pooled[ci]+=source[item.unit_index]*(item.n_cells/total)
                count[ci,labels.index(item.sub_cell_type)]=item.n_cells
        x=pooled[:,remain];sums=x.sum(axis=1);ss=np.einsum('ij,ij->i',x,x);g=x.shape[1]
        expected={key:np.full((n,50),np.nan) for key in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
        for fi,name in enumerate(names):
            ids=np.flatnonzero(np.isin(gene[remain],list(refs[name])));k=len(ids)
            c=coverage[(coverage.view==view)&(coverage.reference==name)].iloc[0]
            assert int(c.remaining_targets)==k and int(c.remaining_background)==g-k
            assert bool(c.qualified)==(k>=20 and g-k>=500)
            if not c.qualified:continue
            foreground=x[:,ids];total=foreground.sum(axis=1);tm=total/k;bm=(sums-total)/(g-k)
            residual_ss=ss-total*total/k-(sums-total)**2/(g-k)
            t=(tm-bm)/np.sqrt(residual_ss/(g-2)*(1/k+1/(g-k)))
            for key,value in [('ULM_score',t),('target_mean',tm),('nontarget_mean',bm),('expressed_target_genes',(foreground>0).sum(axis=1))]:expected[key][:,fi]=value
        for key,values in expected.items():close('all_native_'+key,values,readouts[str(vi)+'_'+key],atol=1e-7)
        q,r=np.linalg.qr(p['design'][tr],mode='reduced')
        for key in ['ULM_score','target_mean','nontarget_mean']:
            values=expected[key];good=np.isfinite(values).all(axis=0)
            residual=np.full_like(values,np.nan)
            coeff=solve_triangular(r,q.T@values[tr][:,good])
            residual[:,good]=values[:,good]-p['design']@coeff
            close('all_QR_'+key,residual,readouts[str(vi)+'_'+key+'_QC_residual'],atol=1e-7)
            for row in funcs[funcs.view==view].itertuples():
                k=int(row.factor.rsplit('F',1)[1])-1;fi=names.index(row.reference)
                for tag,values2 in [('raw',values),('QC_residual',residual)]:
                    actual=getattr(row,key+'_'+tag+'_rho');expected_point=point(scores[ho,k],values2[ho,fi],blocks)
                    close('all2000_function_'+key+'_'+tag,actual,expected_point)
        w=data['model']['weights_'+str(vi)];selected=gene[p['selected_native_columns']]
        for row in weights[weights.view==view].itertuples():
            k=int(row.factor.rsplit('F',1)[1])-1;members=refs[row.reference] if row.reference_kind=='Hallmark' else set(PANELS[row.reference])
            mask=np.isin(selected,list(members));v=w[:,k];num=int(mask.sum());total=np.sum(v**2)
            close('all2320_weight_share',row.L2_fraction,np.sum(v[mask]**2)/total)
            if num:close('all2320_weight_enrichment',row.L2_enrichment,(np.sum(v[mask]**2)/total)/(num/len(v)))
            valid=num>=(10 if row.reference_kind=='Hallmark' else 1) and num<len(v) and v.std()>1e-12 and total>1e-12
            expected_contrast=(np.mean(v[mask])-np.mean(v[~mask]))/np.std(v) if valid else np.nan
            close('all2320_weight_contrast',row.signed_weight_contrast,expected_contrast)
        fractions=count/count.sum(axis=1)[:,None];close('all_subtype_fractions',fractions,models['fractions_'+str(vi)])
        mean=fractions[tr].mean(axis=0);sd=fractions[tr].std(axis=0);sd[sd<1e-12]=1.
        design=(fractions-mean)/sd;centre=scores[tr].mean(axis=0)
        augmented=np.vstack([design[tr],np.eye(len(labels))]);outcomes=np.vstack([scores[tr]-centre,np.zeros((len(labels),10))])
        coef=np.linalg.lstsq(augmented,outcomes,rcond=None)[0];prediction=centre+design@coef;residual=scores-prediction
        independent_residuals.append(residual[ho]);close('all40_ridge_coefficients',coef,models['mixture_coefficients_'+str(vi)])
        close('all8320_ridge_predictions',prediction,models['mixture_prediction_'+str(vi)])
        close('all8320_mixture_residuals',residual,models['mixture_residual_'+str(vi)])
        for k in range(10):
            factor='NSCLC_MFA_20261016_F%02d'%(k+1);row=mixtures[(mixtures.view==view)&(mixtures.factor==factor)].iloc[0]
            r2=1-np.sum(residual[ho,k]**2)/np.sum((scores[ho,k]-centre[k])**2)
            close('all40_mixture_R2',row.holdout_R2,r2)
            close('all40_prediction_rho',row.holdout_prediction_rho,corr(ranks(prediction[ho,k]),ranks(scores[ho,k])))
            vr=variance[(variance.view==view)&(variance.factor==factor)].iloc[0];y=p['X'][tr];z=data['model']['discovery_factor_scores']
            close('all40_variance',vr.isolated_discovery_R2,1-np.sum((y-z[:,[k]]*w[:,k])**2)/np.sum(y**2))
            close('all40_total_variance',vr.all_factor_discovery_R2,1-np.sum((y-z@w.T)**2)/np.sum(y**2))
        for row in subtypes[subtypes.view==view].itertuples():
            k=int(row.factor.rsplit('F',1)[1])-1;j=labels.index(row.author_subtype)
            close('all_subtype_correlations',row.block_rank_rho,point(scores[ho,k],fractions[ho,j],blocks))
            assert row.holdout_zero_cases==int((count[ho,j]==0).sum())
        global_values={}
        for case,group in u.groupby('sampleID'):
            global_values[case]=[np.log1p(group.n_cells.sum()),np.log1p(group.native_total_RNA.sum()),np.average(group.mean_logRNA,weights=group.n_cells),np.average(group.mean_loggenes,weights=group.n_cells)]
        gv=np.array([global_values[s] for s in data['cases'].sampleID])
        tech={str(s):p['quality'][:,j] for j,s in enumerate(p['covariate_names'].astype(str))}
        tech.update({s:gv[:,j] for j,s in enumerate(['global_logN','global_logRNA','global_mean_logRNA','global_mean_loggenes'])})
        tech.update(source_age=data['demographic'].age_numeric.to_numpy(),source_female=data['demographic'].female.to_numpy())
        for row in covariates[covariates.view==view].itertuples():
            k=int(row.factor.rsplit('F',1)[1])-1
            if row.covariate in tech:
                values=tech[row.covariate][ho];close('all_covariate_block_rho',row.block_rank_rho,point(scores[ho,k],values,blocks))
                good=np.isfinite(values);v=values[good];s=scores[ho,k][good]
                value=corr(ranks(v),ranks(s)) if len(v)>=3 and np.ptp(v)>=1e-12 and np.ptp(s)>=1e-12 else np.nan
                close('all_covariate_raw_rho',row.raw_Spearman,value)
            else:
                groups=data['cases'][row.covariate].to_numpy(str)[ho];s=scores[ho,k];mean=s.mean()
                between=sum(np.sum(groups==g)*(s[groups==g].mean()-mean)**2 for g in set(groups))
                close('all_covariate_eta2',row.eta_squared,between/np.sum((s-mean)**2))
    for row in pairs.itertuples():
        k=int(row.factor.rsplit('F',1)[1])-1;a=VIEWS.index(row.view_A);b=VIEWS.index(row.view_B)
        raw=[independent_scores[j][:,k] for j in [a,b]];mix=[independent_residuals[j][:,k] for j in [a,b]]
        close('all60_raw_coordination',row.raw_block_rho,point(*raw,blocks))
        close('all60_mixture_coordination',row.mixture_residual_block_rho,point(*mix,blocks))
        for tag,values in [('sex',raw),('mixture_and_sex',mix)]:
            close('all60_'+tag+'_partial',getattr(row,tag+'_partial_rho'),partial(*values,sex,blocks))
            good=np.isfinite(sex);close('all60_'+tag+'_same_subset',getattr(row,tag+'_same_subset_block_rho'),point(values[0][good],values[1][good],blocks[good]))
    # Verify zero/constant and unidentifiable scenarios without changing real data.
    assert not np.isfinite(point(np.ones(63),np.arange(63),blocks))
    assert not np.isfinite(partial(np.arange(63),np.arange(63),np.zeros(63),blocks))
    result=dict(status='PASS_ALL2000_FUNCTION_ROWS_2320_WEIGHT_ANNOTATIONS_ALL40_RIDGE_AND_ALL60_SENSITIVITY_INDEPENDENT_DIRECT_POOL_OLS_QR',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(terminal),
        numerical_checks=checks,maximum_absolute_errors=errors,full_scope_preserved=True,new_P_or_FDR=False,
        raw_matrix_replayed=False,pooling_scope='All original units supporting all832 case-view rows, all native genes from previously audited cache; not all source cells re-extracted.',
        limitations='Same cohort and shared RNA denominator; descriptive diagnostics are not independent biological or clinical validation.')
    proof.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],finished_UTC=result['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
