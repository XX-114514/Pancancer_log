"""All-factor descriptive interpretation; primary model and 60 P values frozen."""
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import load, pool_native, quota_gate, sha, VIEWS, PANELS, BASIS
from NSCLC_native_MFA_holdout_math_v1 import point
from disjoint_hallmark_ULM_kernel_v1 import score


def rho(a, b, blocks):
    valid = np.isfinite(a) & np.isfinite(b)
    try:
        value, _, _, _, active = point(a[valid], b[valid], blocks[valid])
        return value, int(valid.sum()), active, 'ESTIMABLE'
    except ValueError as error:
        return np.nan, int(valid.sum()), 0, str(error)


def rawrho(a, b):
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3 or np.ptp(a[valid]) < 1e-12 or np.ptp(b[valid]) < 1e-12:
        return np.nan
    return float(np.corrcoef(rankdata(a[valid]), rankdata(b[valid]))[0, 1])


def partial_sex(a, b, sex, blocks):
    good = np.isfinite(sex) & np.isfinite(a) & np.isfinite(b)
    a, b, sex, blocks = a[good], b[good], sex[good], blocks[good]
    base = rho(a, b, blocks)[0]
    design = np.column_stack([*(blocks == k for k in sorted(set(blocks))), sex]).astype(float)
    rank = np.linalg.matrix_rank(design)
    if len(a) < 30 or len(a) - rank < 10 or rank < design.shape[1]:
        return base, np.nan, len(a), len(a)-rank, 'UNIDENTIFIABLE_OR_COVERAGE'
    values = np.column_stack([rankdata(a), rankdata(b)])
    residual = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
    norms = np.linalg.norm(residual, axis=0)
    if (norms < 1e-12).any():
        return base, np.nan, len(a), len(a)-rank, 'CONSTANT_RANK_RESIDUAL'
    return base, float(residual[:, 0] @ residual[:, 1] / np.prod(norms)), len(a), len(a)-rank, 'ESTIMABLE'


def main():
    data = load(); root = data['root']
    out = root / 'NSCLC_native_MFA_interpretation_20261016_v1'
    assert not out.exists(); quota = quota_gate(root); out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__),
        helper_sha256=sha(Path(__file__).with_name('NSCLC_native_MFA_interpretation_inputs_v1.py')),
        ULM_kernel_sha256=sha(Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py')),
        frozen_upstream_evidence=BASIS, scope='All10 factors x all4 views x all50 Hallmarks; all40 technical/composition and all60 pair diagnostics.',
        weight_annotation='Signed target minus complement mean divided by overall weight SD; L2 share and share divided by feature fraction. Hallmark >=10 selected targets, fixed small panels >=1; missing retained.',
        RNA='Cell-weighted native mean(log1pCP10K). Exclude union of all four selected feature sets from target AND background. No renormalization. All50 references retained; >=20 remaining targets and >=500 non-target native genes. Binary-membership ULM plus target/background means and positive-target support.',
        RNA_QC='Train-only OLS on the original complete technical design, applied to raw ULM and both mean readouts; holdout within histology x therapy rank correlations, no P.',
        RNA_low_support_flag='Median number positive targets <10 or median target mean <=0. Retain every row; low support cannot support pathway activation.',
        composition='All original author subtype fractions per view, true zeros retained, discovery-only standardization, constant SD replaced by1 and recorded; ridge alpha1 with free intercept. Holdout R2=1-SSE/SST about discovery score mean; negatives retained.',
        pair_sensitivity='Same63 holdout cases and original blocks, subtract discovery-learned composition prediction for all60. Also rank partial correlations for all60 using full block dummies plus original author female covariate, and both mixture and sex. >=30 finite cases, >=10 residual df, full rank; no silent covariate removal.',
        demographics='Source F/M only, numeric age1..130; unknown remains missing. Associations with own quality, global case quality, age, sex, histology and therapy are descriptive. Source IDs remain author cases.',
        panels=PANELS, quota=quota, new_P_or_FDR=False, clinical_or_model_refit=False,
        limitations=['Gene disjointness is not independent patients; shared full-RNA denominator can couple genes.', 'Stress/lineage RNA panels are diagnostic, not proof of technical artifact or contamination.', 'Composition sensitivity may remove genuine ecology; do not automatically reject it.', 'Weight annotation is not independent function validation.'])
    save('plan_before_any_new_interpretation_effect.json', plan)
    pooled, counts = pool_native(data)
    tr = data['train']; ho = ~tr; cases = data['cases']; blocks = cases.holdout_block.to_numpy(str)[ho]
    genes = data['genes']; remaining = data['remaining']; refs = data['refs']; names = sorted(refs)
    factors = ['NSCLC_MFA_20261016_F%02d' % (k+1) for k in range(10)]
    weight_rows, functional_rows, coverage_rows, mixture_rows, covariate_rows, subtype_rows, model_rows = [], [], [], [], [], [], []
    artifacts = dict(cases=cases.sampleID.to_numpy(str), train=tr, views=np.array(VIEWS), functions=np.array(names))
    readouts = {}
    units = data['units']; global_table = units.groupby('sampleID').agg(n_cells=('n_cells','sum'), RNA=('native_total_RNA','sum'))
    for name in ['mean_logRNA','mean_loggenes']:
        global_table[name] = (units.n_cells * units[name]).groupby(units.sampleID).sum() / global_table.n_cells
    global_table = global_table.reindex(cases.sampleID)
    for vi, view in enumerate(VIEWS):
        p = data['prepared'][vi]; w = data['model']['weights_'+str(vi)]
        scores = data['model']['own_view_scores_'+str(vi)]; artifacts['scores_'+str(vi)] = scores
        selected = genes[p['selected_native_columns']]
        for kind, sets in [('Hallmark', refs), ('fixed_diagnostic_panel', {k:set(v) for k,v in PANELS.items()})]:
            for name, members in sets.items():
                mask = np.isin(selected, list(members)); n = int(mask.sum()); gate = n >= (10 if kind=='Hallmark' else 1) and n < len(mask)
                for k, factor in enumerate(factors):
                    v = w[:,k]; ss = float(v@v); sd = float(v.std()); valid = gate and sd>1e-12 and ss>1e-12
                    share = float(v[mask]@v[mask]/ss) if ss>1e-12 else np.nan
                    weight_rows.append(dict(factor=factor,view=view,reference=name,reference_kind=kind,reference_genes=len(members),native_targets=int(np.isin(genes,list(members)).sum()),selected_targets=n,
                        status='DESCRIPTIVE' if valid else 'INSUFFICIENT_TARGETS_OR_WEIGHT_VARIATION',signed_weight_contrast=float((v[mask].mean()-v[~mask].mean())/sd) if valid else np.nan,L2_fraction=share,L2_enrichment=share/(n/len(mask)) if n and np.isfinite(share) else np.nan))
        target_indices=[]; qualified=[]
        for fi,name in enumerate(names):
            targets=np.flatnonzero(np.isin(genes[remaining],list(refs[name]))); good=len(targets)>=20 and len(remaining)-len(targets)>=500
            coverage_rows.append(dict(view=view,reference=name,reference_genes=len(refs[name]),native_targets=int(np.isin(genes,list(refs[name])).sum()),remaining_targets=len(targets),remaining_background=len(remaining)-len(targets),excluded_union_features=len(data['union']),qualified=good))
            if good: qualified.append(fi); target_indices.append(targets)
        measurements={key:np.full((len(cases),50),np.nan) for key in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
        if target_indices:
            result=score(pooled[vi][:,remaining],target_indices)
            for key in measurements: measurements[key][:,qualified]=result[key]
        for key,value in measurements.items(): readouts[str(vi)+'_'+key]=value
        residuals={}
        for key in ['ULM_score','target_mean','nontarget_mean']:
            residuals[key]=np.full_like(measurements[key],np.nan)
            if qualified:
                b=np.linalg.lstsq(p['design'][tr],measurements[key][tr][:,qualified],rcond=None)[0]
                residuals[key][:,qualified]=measurements[key][:,qualified]-p['design']@b
            readouts[str(vi)+'_'+key+'_QC_residual']=residuals[key]
        for k,factor in enumerate(factors):
            for fi,name in enumerate(names):
                row=dict(view=view,factor=factor,reference=name,qualified=fi in qualified,
                    median_positive_target_genes=float(np.nanmedian(measurements['expressed_target_genes'][ho,fi])) if fi in qualified else np.nan,
                    median_target_mean=float(np.nanmedian(measurements['target_mean'][ho,fi])) if fi in qualified else np.nan)
                row['low_foreground_support']=bool(fi in qualified and (row['median_positive_target_genes']<10 or row['median_target_mean']<=0))
                for key in ['ULM_score','target_mean','nontarget_mean']:
                    for tag,values in [('raw',measurements[key]),('QC_residual',residuals[key])]:
                        r,n,active,status=rho(scores[ho,k],values[ho,fi],blocks)
                        row[key+'_'+tag+'_rho']=r;row[key+'_'+tag+'_status']=status
                functional_rows.append(row)
        fractions=counts[vi]/counts[vi].sum(axis=1,keepdims=True)
        mean=fractions[tr].mean(axis=0);sd=fractions[tr].std(axis=0);constant=sd<1e-12;sd[constant]=1.
        x=(fractions-mean)/sd;centre=scores[tr].mean(axis=0)
        coef=np.linalg.solve(x[tr].T@x[tr]+np.eye(x.shape[1]),x[tr].T@(scores[tr]-centre))
        prediction=centre+x@coef;mix_residual=scores-prediction
        for key,value in [('fractions',fractions),('fraction_mean',mean),('fraction_sd',sd),('fraction_constant',constant),('mixture_coefficients',coef),('mixture_prediction',prediction),('mixture_residual',mix_residual)]:artifacts[key+'_'+str(vi)]=value
        technical={str(n):p['quality'][:,j] for j,n in enumerate(p['covariate_names'].astype(str))}
        technical.update(global_logN=np.log1p(global_table.n_cells.to_numpy()),global_logRNA=np.log1p(global_table.RNA.to_numpy()),global_mean_logRNA=global_table.mean_logRNA.to_numpy(),global_mean_loggenes=global_table.mean_loggenes.to_numpy(),source_age=data['demographic'].age_numeric.to_numpy(),source_female=data['demographic'].female.to_numpy())
        for k,factor in enumerate(factors):
            denom=float(np.sum((scores[ho,k]-centre[k])**2))
            mixture_rows.append(dict(view=view,factor=factor,subtypes=len(mean),constant_discovery_subtypes=int(constant.sum()),holdout_R2=1-float(np.sum(mix_residual[ho,k]**2))/denom if denom>1e-12 else np.nan,holdout_prediction_rho=rawrho(prediction[ho,k],scores[ho,k]),holdout_residual_SD=float(mix_residual[ho,k].std()),holdout_original_SD=float(scores[ho,k].std())))
            for name,values in technical.items():
                r,n,active,status=rho(scores[ho,k],values[ho],blocks)
                covariate_rows.append(dict(view=view,factor=factor,covariate=name,raw_Spearman=rawrho(scores[ho,k],values[ho]),block_rank_rho=r,n_cases=n,active_block_cases=active,status=status))
            for name in ['cancer_type','therapy_group']:
                values=cases[name].to_numpy(str)[ho];overall=float(scores[ho,k].mean());den=float(np.sum((scores[ho,k]-overall)**2))
                between=sum(int((values==s).sum())*(float(scores[ho,k][values==s].mean())-overall)**2 for s in set(values))
                covariate_rows.append(dict(view=view,factor=factor,covariate=name,raw_Spearman=np.nan,block_rank_rho=np.nan,n_cases=int(ho.sum()),active_block_cases=np.nan,status='CATEGORICAL_DESCRIPTIVE',eta_squared=between/den if den>1e-12 else np.nan,groups=len(set(values)),minimum_group_cases=min(int((values==s).sum()) for s in set(values))))
            for j,label in enumerate(data['labels'][view]):
                r,n,active,status=rho(scores[ho,k],fractions[ho,j],blocks)
                subtype_rows.append(dict(view=view,factor=factor,author_subtype=label,block_rank_rho=r,status=status,holdout_zero_cases=int((counts[vi][ho,j]==0).sum())))
            z=data['model']['discovery_factor_scores'][:,k]; y=p['X'][tr]; fit=np.outer(z,w[:,k]);ss=float(np.sum(y*y))
            model_rows.append(dict(view=view,factor=factor,isolated_discovery_R2=1-float(np.sum((y-fit)**2))/ss,all_factor_discovery_R2=1-float(np.sum((y-data['model']['discovery_factor_scores']@w.T)**2))/ss))
    pairs=[];sex=data['demographic'].female.to_numpy()[ho]
    original=pd.read_csv(root/'NSCLC_native_MFA_holdout_20261016_v2/private_all60_internal_holdout_coordination.tsv',sep='\t')
    for k,factor in enumerate(factors):
        for a,b in itertools.combinations(range(4),2):
            s=[artifacts['scores_'+str(j)][ho,k] for j in [a,b]]
            m=[artifacts['mixture_residual_'+str(j)][ho,k] for j in [a,b]]
            raw,n,active,status=rho(*s,blocks);adjusted,_,_,ms=rho(*m,blocks)
            old=original[(original.factor==factor)&(original.view_A==VIEWS[a])&(original.view_B==VIEWS[b])]
            assert len(old)==1 and abs(float(old.iloc[0].block_rank_rho)-raw)<1e-10
            r=dict(factor=factor,view_A=VIEWS[a],view_B=VIEWS[b],n_cases=n,active_block_cases=active,raw_block_rho=raw,mixture_residual_block_rho=adjusted,mixture_status=ms)
            for tag,values in [('sex',s),('mixture_and_sex',m)]:
                same,partial,n,df,status=partial_sex(*values,sex,blocks)
                r[tag+'_same_subset_block_rho']=same;r[tag+'_partial_rho']=partial;r[tag+'_n']=n;r[tag+'_df']=df;r[tag+'_status']=status
            pairs.append(r)
    quota_gate(root)
    tables={'private_all2320_weight_function_and_panel_annotations.tsv':weight_rows,'private_all2000_unused_gene_function_correlations.tsv':functional_rows,'private_all200_RNA_function_coverage.tsv':coverage_rows,'private_all40_subtype_mixture_predictions.tsv':mixture_rows,'private_all_covariate_diagnostics.tsv':covariate_rows,'private_all_author_subtype_diagnostics.tsv':subtype_rows,'private_all40_model_variance_descriptions.tsv':model_rows,'private_all60_composition_and_sex_sensitivity.tsv':pairs}
    for name,rows in tables.items():pd.DataFrame(rows).to_csv(out/name,sep='\t',index=False)
    np.savez_compressed(out/'private_all_function_readouts_and_QC_residuals.npz',**readouts)
    np.savez_compressed(out/'private_all_mixture_models_and_scores.npz',**artifacts)
    data['demographic'].to_csv(out/'private_source_demographics.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values())<256*1024**2
    save('private_terminal.json',dict(status='COMPLETE_ALL_FACTORS_FUNCTION_COMPOSITION_TECHNICAL_DESCRIPTIONS_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,new_P_or_FDR=False,clinical_or_model_refit=False))
    print(json.dumps(dict(status='COMPLETE_INTERPRETATION_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
