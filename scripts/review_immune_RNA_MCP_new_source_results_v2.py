"""Complete accepted aggregate review and all-source/clinical forest figures."""
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


def dot(ax,r,y,color,marker,label=None):
    if not np.isfinite(r.effect):return
    if np.isfinite(r.CI_low) and np.isfinite(r.CI_high):ax.plot([r.CI_low,r.CI_high],[y,y],color=color,lw=1.8)
    ax.plot(r.effect,y,marker=marker,color=color,ms=6,label=label,linestyle='none',markerfacecolor=color if np.isfinite(r.CI_low) else 'white')


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);audit=root/'immune_RNA_MCP_new_source_association_audit_20261019_v3/private_terminal.json';assert sha(audit)=='681cf45ecbfa98f4febb0519309241c698dbb1ce8ac503a2ffdcbf8f472ad1f7';src=root/'immune_RNA_MCP_new_source_associations_20261019_v2';term=json.loads((src/'private_terminal.json').read_text());assert sha(src/'private_terminal.json')==json.loads(audit.read_text())['producer_terminal_sha256']
    for name,v in term['files'].items():assert sha(src/name)==v['sha256']
    m=pd.read_csv(src/'private_all420_registered_main_results.tsv',sep='\t');raw=pd.read_csv(src/'private_all48_unadjusted_clinical_rank_results.tsv',sep='\t');s=pd.read_csv(src/'private_all3360_prespecified_point_sensitivities.tsv',sep='\t');contexts=pd.read_csv(root/'immune_RNA_MCP_new_source_family_inputs_20261019_v2/private_all_new_source_contexts.tsv',sep='\t',keep_default_na=False)
    out=root/'immune_RNA_MCP_new_source_complete_review_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','svg.fonttype':'none'})
    sources=m[m.family.eq('source_marginal')&m.metadata_preeligible].copy();adjusted=m[m.family.eq('source_adjusted')&m.metadata_preeligible].set_index(['context_id','pair','program']);assert len(sources)==12 and len(adjusted)==4
    sources=sources.merge(contexts[['context_id','cancer_label','tissue','phase_class','material_class']],on='context_id',validate='many_to_one');fig,ax=plt.subplots(figsize=(15,9));y=np.arange(len(sources))[::-1];labels=[]
    for j,r in enumerate(sources.itertuples(index=False)):
        pair=r.pair.replace('R_FOXP3_capture','R').replace('E_effector_capture','E').replace('Myeloid','My').replace('|',' / ');source='CAR infusion product' if r.dataset.startswith('CXG_') else r.dataset+' '+r.tissue+' '+r.phase_class;labels.append(source+' | '+pair+' '+r.program+' | n='+str(r.n_units));dot(ax,r,y[j]+.12,'#777777','o','Marginal Spearman (BH12)' if j==0 else None);ax.text(1.04,y[j]+.12,f'{r.BH:.3g}',color='#555555',va='center')
        key=(r.context_id,r.pair,r.program)
        if key in adjusted.index:
            z=adjusted.loc[key];dot(ax,z,y[j]-.12,'#16699b','s','Quality/abundance adjusted (BH4)' if j==0 else None);ax.text(1.25,y[j]-.12,f'{z.BH:.3g}',color='#16699b',va='center')
        else:ax.text(1.25,y[j]-.12,'--',color='#999999',va='center')
    ax.axvline(0,color='black',lw=.8);ax.set_xlim(-1.02,1.02);ax.set_ylim(-.65,len(sources)-.15);ax.set_yticks(y,labels);ax.set_xlabel('Spearman / partial-Spearman correlation; case-bootstrap 95% interval');ax.set_title('2026-10-18 RNA MCP v1: all new-source estimable contexts',loc='left',fontweight='bold',pad=32);ax.text(1.04,len(sources)-.02,'q(marg)',fontsize=9);ax.text(1.25,len(sources)-.02,'q(adj)',fontsize=9);ax.legend(loc='lower left',bbox_to_anchor=(0,-.20),frameon=False,ncol=2);ax.spines[['top','right']].set_visible(False)
    fig.text(.03,.025,'R/E are RNA capture views, not formal Treg/CD8 states. CAR products are not tumour ICI replication. Unknown stage/processing retained.\nAdjusted tests require >=20 units / >=10 residual df. Open markers: interval unavailable (<1900 valid bootstrap draws); --: no adjusted inference.\nRepeated materials collapsed within source/context; contexts sharing patients are not independent replications. Programme learning uncertainty is not included.',fontsize=9)
    fig.subplots_adjust(left=.43,right=.80,top=.88,bottom=.22);fig.savefig(out/'private_all_source_context_forest.png',dpi=180,bbox_inches='tight');fig.savefig(out/'private_all_source_context_forest.svg',bbox_inches='tight');plt.close(fig)
    clinical=m[m.family.str.startswith('clinical')&m.metadata_preeligible].copy();ridx=raw.set_index('model_id');assert len(clinical)==16
    fig,axes=plt.subplots(2,2,figsize=(15,10));epnames={'response3m':'Primary: author 3-month response (R=CR)','ICANS_group':'Secondary: ICANS grade3-4 vs0-2','CRS_max_grade':'Secondary: author maximum CRS grade','prolonged_cytopenia':'Secondary: author cytopenia label'}
    for ax,(ep,title) in zip(axes.flat,epnames.items()):
        part=clinical[clinical.endpoint.eq(ep)];yy=np.arange(4)[::-1];lab=[]
        for j,r in enumerate(part.itertuples(index=False)):
            lab.append(r.program+' '+('shared mean' if r.contrast=='shared_mean' else 'R minus E'));dot(ax,r,yy[j],'#16699b','o');rr=ridx.loc[r.model_id];ax.text(1.03,yy[j],f'{r.BH:.3g} / {rr.BH:.3g}',transform=ax.get_yaxis_transform(),va='center',fontsize=9)
        ax.axvline(0,color='black',lw=.8);ax.set_ylim(-.5,3.6);ax.set_yticks(yy,lab);ax.set_title(title+'; n='+str(int(part.n_units.iloc[0])),loc='left',fontweight='bold',fontsize=11);ax.spines[['top','right']].set_visible(False);ax.text(1.03,3.55,'q adjusted / q raw ranks',transform=ax.get_yaxis_transform(),fontsize=9)
        ax.set_xlabel('Partial Spearman rho' if ep=='CRS_max_grade' else 'Adjusted difference in 0.5(R/SD_R +/- E/SD_E) units')
        if ep=='CRS_max_grade':ax.set_xlim(-1.02,1.02)
    fig.suptitle('CAR infusion products: all fixed RNA-MCP clinical contrasts',x=.03,ha='left',fontweight='bold',fontsize=15)
    fig.text(.03,.025,'Major family: all4 response contrasts. Separate secondary family: all12 toxicity/author-label contrasts. Raw Wilcoxon/Spearman family: all16.\nFull adjustment: both8 original RNA/abundance variables, sex, parent CD4/CD8 proportions. Source is previously used in project; no causal or predictive validation.\nUnknown outcomes stay missing; response is not ORR. Cytopenia timing remains unverified. Rank effects for ordinal CRS do not decompose additively.',fontsize=9)
    fig.subplots_adjust(left=.13,right=.83,top=.89,bottom=.18,hspace=.65,wspace=1.05);fig.savefig(out/'private_all16_clinical_forest.png',dpi=180,bbox_inches='tight');fig.savefig(out/'private_all16_clinical_forest.svg',bbox_inches='tight');plt.close(fig)
    review=dict(reviewed_UTC=datetime.now(timezone.utc).isoformat(),acceptance_sha256=sha(audit),source_main_all_metadata_eligible=m[m.metadata_preeligible].to_dict(orient='records'),all_unadjusted_clinical_metadata_eligible=raw[raw.metadata_preeligible].to_dict(orient='records'),all_predefined_sensitivities=s.to_dict(orient='records'),all_new_source_contexts=contexts.to_dict(orient='records'),
        interpretation='New NSCLC tumour MCP2 R/E marginal coactivity is strong and held under fullcapture/stress points; insufficient cases for full8QC inference, adjacentnormal same-source direction is also positive but does not surviveBH12. CAR both programmes show quality-conditioned R/E coactivity, and MCP2 is weak marginally: do not call it unconditioned tumour-specific replication. None of4 primary or12 secondary clinical contrasts passes its fixed BH family; complete raw Wilcoxon/Spearman family also provides no corrected clinical support. Nominal ICANS hints remain exploratory, not resistance/toxicity mechanisms.',
        next='Close these finite families. Prioritise biological explanation of R/E shared programmes and qualify all previously frozen directional hypotheses in genuinely new sources using existing all-gene sufficient statistics; no clinical endpoint/covariate retuning. First diagnose source-context/quality suppression and receiver readout fidelity; spatial bridge only if measurement fidelity and genuine processing context support it.',goal_complete=False)
    (out/'private_complete_scientific_review_and_next_decision.json').write_text(json.dumps(review,indent=2)+'\n');files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    terminal=dict(status='COMPLETE_ALL_ACCEPTED_SOURCE_CLINICAL_RESULTS_PRIVATE_REVIEW_AND_FOREST_FIGURES_VISUAL_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,files=files,new_P_or_gene_fit=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),files={k:v for k,v in files.items() if k.endswith('.png')})))


if __name__=='__main__':main()
