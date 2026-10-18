"""All registered new-source and clinical tests; no effect-driven selection."""
import json
import os
import hashlib
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import read,wait
from immune_RNA_MCP_new_source_math_v1 import design,estimate,wild_values,bootstrap_indices,bootstrap,leave_one,rank_effect,rank_permutations,adjust,seed


def point(x,y,C,names,kind,marginal=False,binary=False):
    if marginal:
        if len(x)<3 or min(np.ptp(x),np.ptp(y))==0:return np.nan,'INSUFFICIENT_VARIATION'
        return rank_effect(x,y,binary),'DESCRIPTIVE_RANK_POINT'
    m,status=design(x,y,C,names,kind)
    if m is None:return np.nan,status
    if kind=='clinical':value=float(np.linalg.pinv(m['D'],rcond=1e-12)[-1]@m['y'])
    else:
        rx=m['x']-m['Q']@(m['Q'].T@m['x']);ry=m['y']-m['Q']@(m['Q'].T@m['y']);value=float(rx@ry/np.linalg.norm(rx)/np.linalg.norm(ry))
    return value,'ESTIMABLE_POINT'


def rank_boot(x,y,indices,binary=False):
    return np.array([rank_effect(x[ix],y[ix],binary) for ix in indices])


def shared_indices(row,a,cache):
    key=('clinical|'+row.endpoint if row.family.startswith('clinical') else 'source|'+row.context_id)+'|'+row.pair+'|'+hashlib.sha256(a['cases'].tobytes()).hexdigest()
    if key not in cache:cache[key]=bootstrap_indices(a['strata'],key)
    return cache[key],key


def mc_p(values,observed,total_streams):
    k=int((np.abs(values)>=abs(observed)-1e-12).sum());n=len(values);p=(1+k)/(n+1);bound=1. if k==n else float(beta_distribution.ppf(1-.01/total_streams,k+1,n-k));return p,min(1.,(1+n*bound)/(n+1)),k


def summaries(values,loo):
    good=np.isfinite(values);ci=np.quantile(values[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan];finite=np.isfinite(loo)
    return dict(CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstrap=int(good.sum()),LOO_valid=int(finite.sum()),LOO_min=float(loo[finite].min()) if finite.any() else np.nan,LOO_max=float(loo[finite].max()) if finite.any() else np.nan)


def sensitivity_rows(row,a):
    rows=[];marginal=row.family=='source_marginal';clinical=row.family.startswith('clinical');names=a['names'];n=len(a['x'])
    for mode in ['unadjusted_metric_matched','all_captured_own_quality','fixed15_stress_complement','RNA_quality_only','abundance_only','plus_view_lineage_mix','component_A','component_B']:
        effect=np.nan;status='NOT_APPLICABLE';x=a['x'];y=a['y'];C=a['C'];nn=names
        if not row.metadata_preeligible:status='INPUT_NOT_METADATA_PREELIGIBLE'
        elif mode=='unadjusted_metric_matched':
            if clinical and row.kind=='clinical':effect=float(y[x==1].mean()-y[x==0].mean());status='DESCRIPTIVE_GROUP_MEAN_DIFFERENCE'
            else:effect,status=point(x,y,np.empty((n,0)),np.array([],str),row.kind,marginal=True)
        elif mode=='plus_view_lineage_mix' and clinical:effect,status=point(x,y,a['fineC'],a['fine_names'],row.kind)
        elif mode in ['component_A','component_B'] and clinical:
            yy=a[mode]
            if row.kind=='clinical':
                m,st=design(x,y,C,nn,row.kind)
                if m is not None:effect=float(np.linalg.pinv(m['D'],rcond=1e-12)[-1]@yy);status='COMPONENT_OLS_POINT_NOT_INDEPENDENT_CLINICAL_TEST'
                else:status=st
            else:effect,status=point(x,yy,C,nn,row.kind)
        elif mode in ['all_captured_own_quality','fixed15_stress_complement','RNA_quality_only','abundance_only']:
            if mode=='all_captured_own_quality':x=a['allcapture_x'];y=a['allcapture_y'];C=a['allcapture_C']
            elif mode=='fixed15_stress_complement':
                y=a['stress_y']
                if not clinical:x=a['stress_x']
            elif mode in ['RNA_quality_only','abundance_only']:
                columns=np.array([not s.startswith('quality::') or s.split('::')[-1] in (['cellQ','geneQ'] if mode=='RNA_quality_only' else ['fraction','logN']) for s in nn]);C=C[:,columns];nn=nn[columns]
            effect,status=point(x,y,C,nn,row.kind,marginal=marginal)
        rows.append(dict(model_id=row.model_id,mode=mode,effect=effect,status=status,n_units=n,no_new_P=True))
    return rows


def correct(frame,by_family=True):
    result=frame.copy();result['BH']=np.nan;result['BY']=np.nan;result['MC_upper_BH']=np.nan
    groups=result.groupby('family').groups if by_family else {'all':result.index}
    for _,ids in groups.items():
        use=np.array([i for i in ids if result.loc[i,'metadata_preeligible']],int)
        if not len(use):continue
        p=result.loc[use,'P'].fillna(1).to_numpy(float);u=result.loc[use,'MC_upper_P'].fillna(1).to_numpy(float)
        result.loc[use,'BH']=adjust(p);result.loc[use,'BY']=adjust(p,True);result.loc[use,'MC_upper_BH']=adjust(u)
    result.loc[result.P.isna(),['BH','BY','MC_upper_BH']]=np.nan
    result['CI_excludes_zero']=result.CI_low*result.CI_high>0
    result['priority_association']=result.BH.lt(.05)&result.MC_upper_BH.lt(.05)&result.CI_excludes_zero&result.LOO_same_sign
    return result


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,32*1024**2);out=root/'immune_RNA_MCP_new_source_associations_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    ca=root/'immune_RNA_MCP_new_source_family_null_audit_20261019_v1/private_terminal.json';wait(root,out,'audit_immune_RNA_MCP_new_source_null_v1',ca,{'ACCEPTED_ALL_UNIQUE_NULL_SIMULATIONS_RANDOM_STREAMS_AND_PER_DESIGN_FLAGS_NO_TRUE_EFFECTS'})
    _,ia,_,registry,arrays=read(root);calroot=root/'immune_RNA_MCP_new_source_family_null_check_20261019_v1';calmap=pd.read_csv(calroot/'private_complete_original_model_design_mapping.tsv',sep='\t').set_index('model_id');total_streams=int(registry.metadata_preeligible.sum()+registry.loc[registry.family.str.startswith('clinical'),'metadata_preeligible'].sum())
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,input_acceptance_sha256=sha(ia),null_acceptance_sha256=sha(ca),scope='All420registered positions. All metadata-preeligible primary/source/secondary families, including unavailable designs asP1. Complete16unadjusted clinical sensitivity family separately. No source/candidate/contrast/endpoint replacement.',seed=2026101904,draws=19999,bootstrap=2000,CI_min_valid=1900,MC_simultaneous_streams=total_streams,random_storage='Every numerical main/permutation/bootstrap/LOO value,caseorder,shared indices,seeds/reasons; full independent replay before review.',component_interpretation='Binary OLS beta contrast decomposes0.5*(betaA +/-betaB). Ordinal partial-rank effect does NOT decompose additively into component correlations.',limits='Source/context marginal coactivity and within-cohort clinical association, not direct communication or two independently clinical endpoints. Programme learning uncertainty not included; new sources already used elsewhere in project.')
    (out/'plan_before_any_true_new_source_or_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n');rows=[];rawrows=[];sens=[];index_cache={}
    for r in registry.itertuples(index=False):
        a=arrays[r.model_id];row=r._asdict();row.update(effect=np.nan,P=np.nan,MC_upper_P=np.nan,HC3_t=np.nan,CI_low=np.nan,CI_high=np.nan,valid_bootstrap=0,LOO_valid=0,LOO_min=np.nan,LOO_max=np.nan,LOO_same_sign=False,wild_exceedances=np.nan)
        raw=row.copy();raw['family']='unadjusted_clinical';raw['status']='INPUT_NOT_METADATA_PREELIGIBLE';store={};clinical=r.family.startswith('clinical')
        if r.metadata_preeligible:
            ix,key=shared_indices(r,a,index_cache);store.update(indices=ix,cases=a['cases'],bootstrap_key=np.array(key))
            if r.status=='ESTIMABLE' and r.family=='source_marginal':
                e=rank_effect(a['x'],a['y']);vals=rank_permutations(a['x'],a['y'],r.model_id);p,upper,count=mc_p(vals,e,total_streams);b=rank_boot(a['x'],a['y'],ix);loo=np.array([rank_effect(np.delete(a['x'],j),np.delete(a['y'],j)) for j in range(len(a['x']))]);row.update(effect=e,P=p,MC_upper_P=upper,wild_exceedances=count,LOO_same_sign=bool(np.isfinite(loo).all() and (loo*e>0).all()),**summaries(b,loo));store.update(main_statistics=vals,bootstrap_effects=b,bootstrap_reasons=np.where(np.isfinite(b),'ESTIMABLE','CONSTANT_BOOTSTRAP_TARGET'),leave_one=loo)
            elif r.status=='ESTIMABLE':
                assert r.model_id in calmap.index
                if not bool(calmap.loc[r.model_id,'pass_all_null_scenarios']):row['status']='NULL_DIAGNOSTIC_FAILED_NO_ADJUSTED_INFERENCE'
                else:
                    m,st=design(a['x'],a['y'],a['C'],a['names'],r.kind);assert st=='ESTIMABLE'
                    try:g,e=estimate(m)
                    except ValueError as err:row['status']=str(err)
                    else:
                        vals=wild_values(g,19999,seed(r.model_id,'HC3_wild'));p,upper,count=mc_p(vals,g['t'],total_streams);b,reasons=bootstrap(a['x'],a['y'],a['C'],a['names'],r.kind,ix,m['nuisance_rank']);loo,lr=leave_one(a['x'],a['y'],a['C'],a['names'],r.kind,m['nuisance_rank']);row.update(effect=e,P=p,MC_upper_P=upper,HC3_t=g['t'],wild_exceedances=count,LOO_same_sign=bool(np.isfinite(loo).all() and (loo*e>0).all()),**summaries(b,loo));store.update(main_statistics=vals,bootstrap_effects=b,bootstrap_reasons=reasons,leave_one=loo,leave_one_reasons=lr)
            if clinical:
                if min(np.ptp(a['x']),np.ptp(a['y']))>0:
                    e=rank_effect(a['x'],a['y'],r.binary_endpoint);vals=rank_permutations(a['x'],a['y'],r.model_id+'::unadjusted',r.binary_endpoint);p,upper,count=mc_p(vals,e,total_streams);b=rank_boot(a['x'],a['y'],ix,r.binary_endpoint);loo=np.array([rank_effect(np.delete(a['x'],j),np.delete(a['y'],j),r.binary_endpoint) for j in range(len(a['x']))]);raw.update(status='ESTIMABLE',effect=e,P=p,MC_upper_P=upper,wild_exceedances=count,LOO_same_sign=bool(np.isfinite(loo).all() and (loo*e>0).all()),**summaries(b,loo));store.update(unadjusted_statistics=vals,unadjusted_bootstrap_effects=b,unadjusted_leave_one=loo)
                else:raw['status']='CONSTANT_UNADJUSTED_TARGET_OR_READOUT'
            np.savez_compressed(out/('private_'+r.model_id+'_complete_stream.npz'),**store)
        rows.append(row)
        if clinical:rawrows.append(raw)
        sens.extend(sensitivity_rows(r,a))
        if len(rows)%12==0 or len(rows)==len(registry):(out/'private_progress.json').write_text(json.dumps(dict(phase='COMPLETE_FIXED_FAMILIES_NO_PARTIAL_REVIEW',models_completed=len(rows),registered_models=len(registry),partial_effect_review_allowed=False,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    result=correct(pd.DataFrame(rows));unadjusted=correct(pd.DataFrame(rawrows),by_family=False);assert len(result)==420 and len(unadjusted)==48 and len(sens)==420*8
    result.to_csv(out/'private_all420_registered_main_results.tsv',sep='\t',index=False);unadjusted.to_csv(out/'private_all48_unadjusted_clinical_rank_results.tsv',sep='\t',index=False);pd.DataFrame(sens).to_csv(out/'private_all3360_prespecified_point_sensitivities.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<32*1024**2
    terminal=dict(status='COMPLETE_ALL420_NEW_SOURCE_AND_CLINICAL_RESULTS_AND_FULL_RANDOM_STREAMS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,input_acceptance_sha256=sha(ia),null_acceptance_sha256=sha(ca),new_gene_model_fit=False,registered_models=420,independent_patient_or_causal_interaction_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
