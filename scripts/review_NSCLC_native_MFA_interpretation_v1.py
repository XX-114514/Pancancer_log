"""Complete aggregate interpretation; no patient or learned-gene export."""
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load, quota_gate, sha, VIEWS


def main():
    data=load();root=data['root'];src=root/'NSCLC_native_MFA_interpretation_20261016_v1'
    audit=src/'private_independent_interpretation_audit.json'
    assert json.loads(audit.read_text())['status'].startswith('PASS_')
    out=root/'NSCLC_native_MFA_interpretation_review_20261016_v1'
    assert not out.exists();quota_gate(root);out.mkdir(mode=0o700)
    names=['private_all2320_weight_function_and_panel_annotations.tsv','private_all2000_unused_gene_function_correlations.tsv','private_all200_RNA_function_coverage.tsv','private_all40_subtype_mixture_predictions.tsv','private_all_covariate_diagnostics.tsv','private_all60_composition_and_sex_sensitivity.tsv','private_all40_model_variance_descriptions.tsv']
    w,f,c,m,t,p,v=[pd.read_csv(src/name,sep='\t') for name in names]
    factors=sorted(p.factor.unique());pairs=list(dict.fromkeys(zip(p.view_A,p.view_B)))
    labels={'author_CD8':'CD8','author_Treg':'Treg','author_Macrophage':'Macrophage','author_B_Plasma':'B / Plasma'}
    pair_labels=[labels[a]+' - '+labels[b] for a,b in pairs]
    summaries=[]
    for factor in factors:
        pp=p[p.factor==factor];mm=m[m.factor==factor];tt=t[t.factor==factor];ff=f[f.factor==factor];vv=v[v.factor==factor]
        one=dict(factor=factor,raw_pair_min=float(pp.raw_block_rho.min()),raw_pair_max=float(pp.raw_block_rho.max()),
            mixture_pair_min=float(pp.mixture_residual_block_rho.min()),mixture_pair_max=float(pp.mixture_residual_block_rho.max()),
            sex_pair_min=float(pp.sex_partial_rho.min()),sex_pair_max=float(pp.sex_partial_rho.max()),
            mixture_and_sex_pair_min=float(pp.mixture_and_sex_partial_rho.min()),mixture_and_sex_pair_max=float(pp.mixture_and_sex_partial_rho.max()),
            max_abs_female_rho=float(tt[tt.covariate=='source_female'].block_rank_rho.abs().max()),
            max_abs_age_rho=float(tt[tt.covariate=='source_age'].block_rank_rho.abs().max()),
            max_mixture_holdout_R2=float(mm.holdout_R2.max()),min_mixture_holdout_R2=float(mm.holdout_R2.min()),
            min_isolated_discovery_R2=float(vv.isolated_discovery_R2.min()),max_isolated_discovery_R2=float(vv.isolated_discovery_R2.max()))
        for view in VIEWS:
            for tag,values in [('weight',w[(w.factor==factor)&(w.view==view)&(w.reference_kind=='Hallmark')&w.signed_weight_contrast.notna()]),('unused_RNA',ff[(ff.view==view)&ff.qualified&~ff.low_foreground_support&ff.ULM_score_QC_residual_rho.notna()])]:
                col='signed_weight_contrast' if tag=='weight' else 'ULM_score_QC_residual_rho'
                values=values.assign(absolute=values[col].abs()).sort_values(['absolute','reference'],ascending=[False,True]).head(3)
                one[view+'_'+tag+'_top3']='; '.join('%s:%.3f'%(r.reference,float(getattr(r,col))) for r in values.itertuples())
        summaries.append(one)
    summary=pd.DataFrame(summaries);summary.to_csv(out/'private_all10_factor_interpretation_summary.tsv',sep='\t',index=False)
    # Complete60 at each diagnostic condition, not a selected subset.
    fig,axes=plt.subplots(2,2,figsize=(16,12),sharey=True)
    conditions=[('raw_block_rho','Original within-block rank correlation'),('mixture_residual_block_rho','Discovery-trained subtype mixture removed'),('sex_partial_rho','Rank partial correlation controlling author sex'),('mixture_and_sex_partial_rho','Subtype mixture and author sex')]
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d9d9d9')
    for ax,(column,title) in zip(axes.ravel(),conditions):
        values=np.array([[float(p[(p.factor==factor)&(p.view_A==a)&(p.view_B==b)].iloc[0][column]) for a,b in pairs] for factor in factors])
        im=ax.imshow(values,vmin=-1,vmax=1,cmap=cmap,aspect='auto')
        ax.set_title(title,fontsize=11);ax.set_xticks(range(6),pair_labels,rotation=35,ha='right',fontsize=8)
        ax.set_yticks(range(10),[x.rsplit('_',1)[1] for x in factors])
        for i in range(10):
            for j in range(6):ax.text(j,i,'%.2f'%values[i,j] if np.isfinite(values[i,j]) else 'NA',ha='center',va='center',fontsize=8)
    fig.suptitle('All10 native lung factors / all6 cell-pool pairs\nDescriptive sensitivity only: original model and60 inference slots unchanged',fontsize=13)
    fig.tight_layout(rect=[0,0,0.94,0.94]);fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.025,pad=.04)
    fig.savefig(out/'private_all60_sensitivity.png',dpi=140,bbox_inches='tight');plt.close(fig)
    # Full continuous technical/demographic diagnostics and mixture prediction R2.
    covs=['source_female','source_age','logN','mean_logRNA','mean_loggenes','symbol_pct_mito','symbol_pct_ribo','fraction_within_broad_parent','global_logN','global_logRNA','global_mean_logRNA','global_mean_loggenes']
    fig,axes=plt.subplots(2,2,figsize=(17,12),sharey=True)
    for ax,view in zip(axes.ravel(),VIEWS):
        values=np.full((10,len(covs)+1),np.nan)
        for i,factor in enumerate(factors):
            values[i,0]=m[(m.view==view)&(m.factor==factor)].iloc[0].holdout_R2
            for j,cov in enumerate(covs):
                z=t[(t.view==view)&(t.factor==factor)&(t.covariate==cov)]
                if len(z):values[i,j+1]=z.iloc[0].block_rank_rho
        im=ax.imshow(values,cmap=cmap,vmin=-1,vmax=1,aspect='auto');ax.set_title(labels[view])
        ax.set_xticks(range(len(covs)+1),['mixture holdout R2']+covs,rotation=80,fontsize=7)
        ax.set_yticks(range(10),[x.rsplit('_',1)[1] for x in factors])
        for i in range(10):
            for j in range(len(covs)+1):
                if np.isfinite(values[i,j]):ax.text(j,i,'%.2f'%values[i,j],ha='center',va='center',fontsize=6)
    fig.suptitle('All40 own-view scores: composition predictability and technical/demographic diagnostics\nR2 retains negative values; colors clipped to[-1,1]. Gray=not measured/constant; no newP',fontsize=12)
    fig.tight_layout(rect=[0,0,0.94,.94]);fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.025,pad=.04)
    fig.savefig(out/'private_all40_background_diagnostics.png',dpi=140,bbox_inches='tight');plt.close(fig)
    # Every50 function is visible for each view; zero support is marked, not hidden.
    functions=sorted(f.reference.unique())
    for view in VIEWS:
        a=np.full((50,10),np.nan);b=a.copy();low=np.zeros_like(a,bool)
        for i,name in enumerate(functions):
            for j,factor in enumerate(factors):
                wr=w[(w.view==view)&(w.factor==factor)&(w.reference==name)&(w.reference_kind=='Hallmark')].iloc[0]
                fr=f[(f.view==view)&(f.factor==factor)&(f.reference==name)].iloc[0]
                a[i,j]=wr.signed_weight_contrast;b[i,j]=fr.ULM_score_QC_residual_rho;low[i,j]=bool(fr.low_foreground_support)
        fig,axes=plt.subplots(1,2,figsize=(16,15),sharey=True)
        for ax,values,title,limit in zip(axes,[a,b],['Learned weight contrast (annotation)','Unused-gene RNA / technical-residual rank rho'],[3,1]):
            im=ax.imshow(values,aspect='auto',cmap=cmap,vmin=-limit,vmax=limit)
            ax.set_title(title,fontsize=11);ax.set_xticks(range(10),[x.rsplit('_',1)[1] for x in factors],fontsize=8)
            ax.set_yticks(range(50),[x.removeprefix('HALLMARK_') for x in functions],fontsize=8)
            fig.colorbar(im,ax=ax,fraction=.03,pad=.025)
        axes[0].tick_params(labelleft=True);axes[1].tick_params(labelleft=False)
        for i,j in zip(*np.where(low)):axes[1].text(j,i,'x',ha='center',va='center',fontsize=6,color='black')
        fig.suptitle(labels[view]+': all10 factors x all50 reference functions\nGray=coverage/variation unavailable; x=low target RNA support. Same cohort, no newP or independent validation.',fontsize=12)
        fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(out/('private_'+view+'_all50_function_map.png'),dpi=130,bbox_inches='tight');plt.close(fig)
    report=['# Complete native lung MFA interpretation','',
        'All10 factors, all4 views, all50 frozen functions, all40 mixture predictions and all60 sensitivities retained. Audit: '+sha(audit),
        'No newP/FDR, clinical test, model refit, formal state label or oldCM bridge.',
        'Function rankings below are descriptive organization of the complete tables, not new selected hypothesis tests.',
        'Signs follow the frozen factor orientation. RNA functions exclude every selected model feature across every view.',
        'Composition and sex sensitivity diagnose sources of variation; attenuation can still represent meaningful composition ecology.',
        '', '```text', summary.to_string(index=False), '```', '',
        'Interpretation limits: library normalization couples remaining genes; same author cases are not independent external patients; these post-treatment data do not establish baseline prediction. Small fixed lineage/stress panels do not prove ambient RNA or contamination.']
    (out/'PRIVATE_COMPLETE_FACTOR_INTERPRETATION.md').write_text('\n'.join(report)+'\n')
    scope=dict(created_UTC=datetime.now(timezone.utc).isoformat(),basis_audit_sha256=sha(audit),
        question='Before any clinical inference, inspect complete author clinical availability, sampleID role and timing; nominate biological candidates only after reviewing all-factor composition/sex/function evidence. Keep original60 unchanged; do not refit merely to rescue their P values.',
        new_clinical_test_or_gene_list_export=False)
    (out/'PRIVATE_NEXT_QUESTION.json').write_text(json.dumps(scope,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    terminal=dict(status='COMPLETE_ALL_FACTOR_INTERPRETATION_AGGREGATE_REVIEW_VISUAL_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),basis_audit_sha256=sha(audit),files=files)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
