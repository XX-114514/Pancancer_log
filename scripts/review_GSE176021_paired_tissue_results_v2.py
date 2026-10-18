"""Complete aggregate review/figures and descriptive covariance identities."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def dot(ax,row,y,color):
    if np.isfinite(row.CI_low) and np.isfinite(row.CI_high):ax.plot([row.CI_low,row.CI_high],[y,y],color=color,lw=2)
    ax.plot(row.effect,y,'o',color=color,ms=7)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);audit=root/'GSE176021_paired_tissue_association_audit_20261019_v1/private_terminal.json';assert sha(audit)=='abcb336145d958a729b6eb47d4065f6fe70f4e5c9b20ad664d2b260e099deec4';src=root/'GSE176021_paired_tissue_associations_20261019_v1';t=json.loads((src/'private_terminal.json').read_text())
    for n,v in t['files'].items():assert sha(src/n)==v['sha256']
    out=root/'GSE176021_paired_tissue_complete_review_20261019_v2';assert not out.exists();out.mkdir(mode=0o700);main=pd.read_csv(src/'private_all6_primary_paired_tissue_results.tsv',sep='\t');sens=pd.read_csv(src/'private_all54_point_sensitivity_results.tsv',sep='\t');co=pd.read_csv(src/'private_all30_coactivity_context_points_only.tsv',sep='\t');colors={'MCP1':'#777777','MCP2':'#146e9d'};plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(15,7));a=main[main.outcome.ne('RE_tissue_cochange')]
    for y,r in zip([3,2,1,0],a.itertuples(index=False)):
        dot(axes[0],r,y,colors[r.program]);axes[0].text(1.03,y,f'{r.BH:.3g} / {r.sign_BH:.3g}',transform=axes[0].get_yaxis_transform(),va='center',fontsize=10)
    axes[0].set_yticks([3,2,1,0],[r.program+' '+('author Treg' if r.outcome.startswith('R_') else 'author CD8 union') for r in a.itertuples(index=False)]);axes[0].set_ylim(-.5,3.6);axes[0].axvline(0,color='black',lw=.8);axes[0].set_xlabel('Median tumour-minus-normal difference\n(original discovery case-SD units)');axes[0].set_title('Separate lineage shifts; paired Wilcoxon',loc='left',fontweight='bold',pad=30);axes[0].text(1.03,3.58,'BH6 / sign BH4',transform=axes[0].get_yaxis_transform(),fontsize=10)
    b=main[main.outcome.eq('RE_tissue_cochange')]
    for y,r in zip([1,0],b.itertuples(index=False)):
        dot(axes[1],r,y,colors[r.program]);axes[1].text(1.03,y,f'{r.BH:.3g} / {r.MC_upper_BH:.3g}',transform=axes[1].get_yaxis_transform(),va='center',fontsize=10)
    axes[1].set_yticks([1,0],b.program);axes[1].set_ylim(-.3,1.3);axes[1].set_xlim(-1.04,1.04);axes[1].axvline(0,color='black',lw=.8);axes[1].set_xlabel('Spearman correlation of paired tissue differences');axes[1].set_title('Treg/CD8 tissue cochange; case permutations',loc='left',fontweight='bold',pad=30);axes[1].text(1.03,1.29,'BH6 / MC upper BH6',transform=axes[1].get_yaxis_transform(),fontsize=10)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE176021: all six frozen RNA-MCP paired-tissue questions (10 cases)',x=.03,ha='left',fontweight='bold',fontsize=15);fig.text(.03,.04,'All six primary tests share one BH/BY family; intervals: 2,000 whole-case-pair bootstrap draws. Four direction-only sign tests are sensitivity only.\nA strong cochange does not imply both lineage scores rise on average. Author reference unions are not new formal state annotations.\nFrozen model quality correction differs from complete case-level confounder adjustment. Stage metadata remains unknown; no MPR or direct-communication claim.',fontsize=10);fig.subplots_adjust(left=.16,right=.82,top=.77,bottom=.28,wspace=1.25)
    for ext in ['png','svg']:fig.savefig(out/('private_all6_paired_tissue_forest.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    scopes=['legacy_capture_anchor','author_confirmed_capture','all_author_lineage'];scnames={'legacy_capture_anchor':'original capture','author_confirmed_capture':'author-filtered capture','all_author_lineage':'all author lineage'};cs={'legacy_capture_anchor':'#777777','author_confirmed_capture':'#d58a19','all_author_lineage':'#146e9d'};tissues=['tumor','adjacent normal lung','LN','brain metastases','blood'];fig,axes=plt.subplots(1,2,figsize=(17,13))
    for ax,pr in zip(axes,['MCP1','MCP2']):
        labels=[];ys=[]
        for ti,tissue in enumerate(tissues):
            for si,scope in enumerate(scopes):
                r=co[co.scope.eq(scope)&co.tissue.eq(tissue)&co.program.eq(pr)].iloc[0];y=14-ti*3-si;ys.append(y);labels.append(tissue+' | '+scnames[scope]+' | n='+str(r.case_units))
                if np.isfinite(r.rho):
                    if np.isfinite(r.LOO_min):ax.plot([r.LOO_min,r.LOO_max],[y,y],color=cs[scope],lw=1.5,alpha=.65)
                    ax.plot(r.rho,y,'o',color=cs[scope],ms=6)
                else:ax.text(.4,y,'unavailable',transform=ax.get_yaxis_transform(),color='#999999',va='center',fontsize=9)
        ax.axvline(0,color='black',lw=.7);ax.set_xlim(-1.04,1.04);ax.set_ylim(-.6,14.6);ax.set_yticks(ys,labels,fontsize=9);ax.set_xlabel('R/E coactivity point\nline: leave-one-case range');ax.set_title(pr,loc='left',fontweight='bold');ax.spines[['top','right']].set_visible(False)
    fig.suptitle('GSE176021: all 30 scope/tissue coactivity positions, descriptive only',x=.03,ha='left',fontsize=15,fontweight='bold');fig.text(.03,.035,'No new P values or confidence intervals. LOO lines are sensitivity ranges, not 95% confidence intervals. n<10 is sparse description only.\nTissue-only pooling differs from the earlier source x cancer x tissue x stage strata; these points do not replace the closed earlier tests.\nTreg/CD8 author unions and the original exclusive RNA captures are distinct measurement definitions. Missing author or insufficient-case positions remain unavailable.',fontsize=10);fig.subplots_adjust(left=.27,right=.97,top=.91,bottom=.14,wspace=1.8)
    for ext in ['png','svg']:fig.savefig(out/('private_all30_scope_coactivity_points.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    # Covariance identities are descriptive, and do not decompose Spearman rho.
    case=pd.read_csv(root/'GSE176021_paired_tissue_inputs_20261019_v1/private_all_scope_case_tissue_values.tsv',sep='\t');reg=pd.read_csv(root/'GSE176021_paired_tissue_inputs_20261019_v1/private_all6_registered_primary_tests.tsv',sep='\t');sens_in=pd.read_csv(root/'GSE176021_paired_tissue_inputs_20261019_v1/private_all_prespecified_case_sensitivity_inputs.tsv',sep='\t');sd=pd.read_csv(root/'GSE176021_paired_tissue_inputs_20261019_v1/private_original_discovery_case_scales.tsv',sep='\t').set_index(['view','program']).discovery_case_sd;decomp=[]
    for pr in ['MCP1','MCP2']:
        ids=reg[reg.program.eq(pr)].model_id;first=np.load(root/'GSE176021_paired_tissue_inputs_20261019_v1'/('private_'+ids.iloc[0]+'_inputs.npz'));cases=first['cases'];d=case[case.scope.eq('all_author_lineage')&case['mode'].eq('primary')&case.program.eq(pr)].set_index(['case_block','tissue','view']);v=np.array([[d.loc[(c,t,view),'value']/sd.loc[(view,pr)] for t in ['adjacent normal lung','tumor'] for view in ['R_FOXP3_capture','E_effector_capture']] for c in cases]);cov=np.cov(v,rowvar=False,ddof=1);delta=v[:,2:]-v[:,:2];delta_cov=float(np.cov(delta,rowvar=False,ddof=1)[0,1]);terms={'tumour_cross_lineage_covariance':float(cov[2,3]),'normal_cross_lineage_covariance':float(cov[0,1]),'minus_tumour_R_normal_E':float(-cov[2,1]),'minus_normal_R_tumour_E':float(-cov[0,3])};assert abs(sum(terms.values())-delta_cov)<1e-10
        shape=[]
        for outc in ['R_tissue_shift','E_tissue_shift']:
            mid=reg[(reg.program==pr)&reg.outcome.eq(outc)].model_id.iloc[0];s=sens_in[sens_in.model_id.eq(mid)].pivot(index='case_block',columns='mode',values='y').reindex(cases);raw=s.unresidualized_weighted_RNA.to_numpy();RNA=-(s.unresidualized_weighted_RNA-s.RNA_quality_only_background).to_numpy();composition=-(s.RNA_quality_only_background-s.primary).to_numpy();parts=np.c_[raw,RNA,composition];assert np.max(np.abs(parts.sum(axis=1)-s.primary.to_numpy()))<1e-10;shape.append(parts)
        cross=(shape[0]-shape[0].mean(axis=0)).T@(shape[1]-shape[1].mean(axis=0))/(len(cases)-1);assert abs(cross.sum()-delta_cov)<1e-10;decomp.append({'program':pr,'case_units':10,'Pearson_delta_rho':float(np.corrcoef(delta.T)[0,1]),'delta_covariance':delta_cov,'paired_reference_covariance_terms':terms,'signal_RNA_correction_composition_cross_covariance_private_matrix':cross.tolist(),'note':'Algebraic Pearson/covariance description, not Spearman decomposition or proof of causal/shared-upstream origin.'})
    review={'reviewed_UTC':datetime.now(timezone.utc).isoformat(),'acceptance_sha256':sha(audit),'all6_primary':main.to_dict(orient='records'),'all54_predefined_point_sensitivities':sens.to_dict(orient='records'),'all30_context_points_no_new_P':co.to_dict(orient='records'),'private_covariance_identities':decomp,'interpretation':'MCP2 shows strong paired Treg/CD8 tissue cochange, CI and leave-one direction retained; independent lineage shifts are not corrected-significant. Positive points under full-cell/stress/own/background/unresidualized checks support a source-context coadaptation candidate, not all-patient upregulation, MPR prediction or direct receptor signaling. Different legacy pooled tissue point vs earlier strict context is a denominator/context difference, not silently replaced estimate.','next':'Close six-question family. Explain shared upstream/functional programmes and baseline-reference/nuisance contribution using accepted whole-author RNA; directional LR/receiver response requires separate fixed full family and appropriate cells. Missing myeloid in a T-enriched source is not evidence of no myeloid involvement. MPR remains unmapped.'}
    (out/'private_complete_scientific_review_and_covariance_diagnostics.json').write_text(json.dumps(review,indent=2)+'\n');files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2;d={'status':'COMPLETE_ALL6_54_30_PRIVATE_REVIEW_COVARIANCE_IDENTITIES_FIGURES_VISUAL_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'files':files,'new_P_or_gene_fit':False}
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json'),'PNG':{k:v for k,v in files.items() if k.endswith('.png')}}))


if __name__=='__main__':main()
