"""Private full-family figures and programme interpretation after complete audit."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);out=root/'NSCLC_communication_tensor_family_review_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_family_review.json').write_text(json.dumps({'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'scope':'All20main slots, all60background sensitivities, all20secondary Wilcoxon, all15model objectives and allrestart stability. No positive-only plots. Learned LR identities/weights and biological interpretation remain server private. Figures require actual visual review after creation, not claimed by this producer.'},indent=2)+'\n')
    accepted=wait(root,out,'audit_NSCLC_communication_tensor_family_v2',root/'NSCLC_communication_tensor_family_audit_20261019_v2/private_terminal.json',{'ACCEPTED_COMPLETE_LR_SUPPORT_ALL15_OBJECTIVES_ALL20_HELDOUT_PROJECTIONS_INPUTS_AND_INFERENCE_REPLAY'})
    src=root/'NSCLC_communication_tensor_family_20261019_v1';assert sha(src/'private_terminal.json')==accepted['producer_terminal_sha256']
    read=lambda n:pd.read_csv(src/n,sep='\t')
    results=read('private_all20_heldout_clinical_results.tsv').set_index('model_id');sensitivity=read('private_all60_background_sensitivity_effects.tsv').set_index(['model_id','sensitivity']);ranks=read('private_all20_secondary_Wilcoxon_results.tsv').set_index('model_id');stability=read('private_all_restart_factor_stability.tsv');models=read('private_all15_CP_models.tsv');registry=read('private_all20_factor_registry.tsv').set_index('model_id');project=read('private_all1740_factor_context_projections.tsv');assert len(results)==len(ranks)==len(registry)==20 and len(models)==15
    axes=np.load(src/'private_tensor_axes.npz');states=axes['states'].tolist();vocab=json.loads((src/'private_selected_LR_vocabulary.json').read_text());slots=results.index.tolist();full=[];private_top=[]
    for slot in slots:
        rr=registry.loc[slot];rank=int(rr['rank']);j=int(rr.factor)-1;f=dict(np.load(src/('private_selected_K%d_factors_and_projection.npz'%rank)));B=f['LR'][:,j];U=f['sender'][:,j];V=f['receiver'][:,j];b=B/B.sum();u=U/U.sum();v=V/V.sum();s=stability[(stability['rank']==rank)&(stability.factor==j+1)&(stability.restart!=rr.best_restart)].mean_four_axis_cosine
        r=results.loc[slot];w=ranks.loc[slot];q=project[project.model_id.eq(slot)];row=dict(model_id=slot,source='GSE243013',rank=rank,factor=j+1,author_records=int(r.author_records),reviewed_independent_patients=0,status=r.status,effect=r.effect,CI_low=r.CI_low,CI_high=r.CI_high,BH20=r.BH20,BY20=r.BY20,wild_P=r.wild_P,bootstrap_valid=int(r.bootstrap_valid),median_other_restart_cosine=float(s.median()),min_other_restart_cosine=float(s.min()),median_heldout_reconstruction_error=float(q.relative_observed_projection_error.median()),LR_effective_number=float(1/np.sum(b*b)),top1_LR_mass_share=float(np.max(b)),top5_LR_mass_share=float(np.sort(b)[-5:].sum()),secondary_Wilcoxon_rank_biserial=w.rank_biserial,secondary_Wilcoxon_BH20=w.secondary_BH20)
        for name in ['unadjusted','clinical_only','clinical_RNA_no_factor_abundance']:row[name+'_effect']=sensitivity.loc[(slot,name)].effect
        full.append(row)
        # All20programmes receive identical descriptive top-loading summaries.
        for kind,weights,labels in [('LR',b,[d['canonical_resource_LR_id'] for d in vocab]),('sender',u,states),('receiver',v,states)]:
            for pos,index in enumerate(np.argsort(-weights,kind='stable')[:10]):
                item=dict(model_id=slot,axis=kind,position=pos+1,label=labels[index],normalized_mass=float(weights[index]),interpretation='descriptive_factor_loading_not_causal_dominance')
                if kind=='LR':item.update(ligand_complex='_'.join(vocab[index]['ligand_components']),receptor_complex='_'.join(vocab[index]['receptor_components']))
                private_top.append(item)
    full=pd.DataFrame(full);full.to_csv(out/'private_all20_source_aggregate_review.tsv',sep='\t',index=False);pd.DataFrame(private_top).to_csv(out/'private_all20_descriptive_top_loading_interpretation.tsv',sep='\t',index=False)
    os.environ['MPLCONFIGDIR']=str(out/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    figures=[]
    def finish(fig,name):
        fig.tight_layout()
        for suffix in ['png','svg']:
            p=out/(name+'.'+suffix);fig.savefig(p,dpi=150,bbox_inches='tight');figures.append({'name':p.name,'sha256':sha(p),'bytes':p.stat().st_size})
        plt.close(fig)
    y=np.arange(20);fig,ax=plt.subplots(figsize=(10,11));ax.axvline(0,color='.55',lw=.8)
    for col,color,offset,label in [('unadjusted_effect','.72',-.17,'Unadjusted'),('clinical_only_effect','#2c7fb8',-.06,'Clinical background'),('clinical_RNA_no_factor_abundance_effect','#e0902e',.06,'Clinical + RNA background')]:ax.scatter(full[col],y+offset,s=17,c=color,label=label)
    for i,r in full.iterrows():
        if np.isfinite(r.effect):
            color='#b2182b' if r.BH20<.05 else '#222222';ax.plot([r.CI_low,r.CI_high],[i+.19,i+.19],c=color,lw=1.3);ax.scatter([r.effect],[i+.19],s=28,c=color)
        else:ax.text(0,i+.19,'not estimable',ha='center',fontsize=8)
    ax.set_yticks(y,slots);ax.invert_yaxis();ax.set_xlabel('Activity difference: pCR/MPR minus nonMPR (frozen discovery SD)');ax.set_title('GSE243013 | all20 heldout communication-factor questions\nMain: clinical + RNA + factor-weighted abundance; 95% bootstrap CI');ax.legend(loc='best',fontsize=8);finish(fig,'all20_main_and_background_forest')
    fig,ax=plt.subplots(figsize=(9,10));ax.axvline(0,c='.55',lw=.8);ax.scatter(full.secondary_Wilcoxon_rank_biserial,y,c=['#b2182b' if q<.05 else '.35' for q in full.secondary_Wilcoxon_BH20],s=28);ax.set_yticks(y,slots);ax.invert_yaxis();ax.set_xlim(-1.05,1.05);ax.set_xlabel('Unadjusted rank-biserial effect (positive = higher pCR/MPR)');ax.set_title('All20 secondary Wilcoxon rank-sum comparisons\nSeparate BH20 family; no abundance or RNA adjustment');finish(fig,'all20_secondary_Wilcoxon')
    matrix=np.column_stack([full.author_records/87,full.bootstrap_valid/2000,full.median_other_restart_cosine,1-full.median_heldout_reconstruction_error]);fig,ax=plt.subplots(figsize=(8,11));im=ax.imshow(matrix,aspect='auto',vmin=0,vmax=1,cmap='viridis');ax.set_yticks(y,slots);ax.set_xticks(np.arange(4),['Heldout N /87','Valid bootstrap /2000','Other-start cosine','1 - heldout error'],rotation=25,ha='right')
    for i in range(20):
        for j in range(4):ax.text(j,i,'%.2f'%matrix[i,j],ha='center',va='center',color='white' if matrix[i,j]<.5 else 'black',fontsize=8)
    ax.set_title('All20 coverage and model diagnostics\nAuthor records are not independently reviewed patients');fig.colorbar(im,ax=ax,shrink=.5);finish(fig,'all20_coverage_stability_and_projection')
    fig,ax=plt.subplots(figsize=(8,5))
    for start,g in models.groupby('restart'):ax.plot(g['rank'],g.observed_relative_error,'o-',label='Fixed start%d'%start)
    ax.set_xticks(range(2,7));ax.set_xlabel('Fixed CP rank');ax.set_ylabel('Observed discovery relative reconstruction error');ax.set_title('All15 discovery fits | no clinical rank selection');ax.legend();finish(fig,'all15_fixed_model_reconstruction')
    report={'status':'PRIVATE_COMPLETE_FAMILY_REVIEW_PRODUCED_ACTUAL_VISUAL_INSPECTION_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'complete_audit_sha256':sha(root/'NSCLC_communication_tensor_family_audit_20261019_v2/private_terminal.json'),'models':15,'clinical_slots':20,'figures':figures,'complete_aggregate_table_sha256':sha(out/'private_all20_source_aggregate_review.tsv'),'private_loading_interpretation_sha256':sha(out/'private_all20_descriptive_top_loading_interpretation.tsv'),'limits':'Post-neoadjuvant association, source previously explored, author records not independently reviewed patients. Factor loading share is descriptive, not causal communication dominance. Original9CM and RNA MCP identities are distinct. All learned names/weights and effect tables remain private; source aggregate review only under current user authorization.'}
    (out/'private_complete_review.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<8*1024**2
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':report['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
