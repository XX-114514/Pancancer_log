"""Private all20 report: distinguish original inference from post-result checks."""
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
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);out=root/'NSCLC_communication_tensor_family_review_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_revised_review.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All20original slots plus complete independently audited post-result intervals/backgrounds/descriptive redundancy. V1visual review identified misleading bootstrap-CI title because all original percentile intervals unavailable; preserveV1, explicitly label supplementary fixed-design wild-t intervals, unavailable originals, secondary Wilcoxon and post-result checks. No underlying main inference changed. No individual or learned gene vectors displayed.'),indent=2)+'\n')
    accepted=wait(root,out,'audit_NSCLC_tensor_complete_family_diagnostics_v1',root/'NSCLC_tensor_complete_family_diagnostics_audit_20261019_v1/private_terminal.json',{'ACCEPTED_ALL20_SUPPLEMENTARY_WILD_T_INTERVALS_ALL190_PAIRS_60_BACKGROUNDS_AND_BOOTSTRAP_FAILURE_DIAGNOSTICS'})
    supp=root/'NSCLC_tensor_complete_family_diagnostics_20261019_v1';assert sha(supp/'private_terminal.json')==accepted['producer_terminal_sha256'];original=root/'NSCLC_communication_tensor_family_20261019_v1';old=root/'NSCLC_communication_tensor_family_review_20261019_v1';oldreceipt=json.loads((old/'private_terminal.json').read_text());assert sha(old/'private_all20_source_aggregate_review.tsv')==oldreceipt['complete_aggregate_table_sha256']
    read=lambda p:pd.read_csv(p,sep='\t');main=read(old/'private_all20_source_aggregate_review.tsv').set_index('model_id');iv=read(supp/'private_all20_supplementary_fixed_design_intervals.tsv').set_index('model_id');roles=read(supp/'private_all20_original_author_role_and_redundancy.tsv').set_index('model_id');sens=read(supp/'private_all60_total_signal_and_endpoint_RNA_sensitivities.tsv').set_index(['model_id','sensitivity']);models=read(original/'private_all15_CP_models.tsv');assert len(main)==len(iv)==len(roles)==20 and len(sens)==60
    slots=main.index.tolist();np.testing.assert_allclose(main.effect,iv.loc[slots].effect,atol=1e-12);full=main.copy()
    for field in iv.columns:
        if field not in full.columns:full[field]=iv.loc[slots,field]
    for field in roles.columns:full[field]=roles.loc[slots,field]
    for label in ['total_LR_and_capture','endpoint_RNA_quality','both']:
        for field in ['effect','HC3_t_CI_low','HC3_t_CI_high','status']:full[label+'__'+field]=[sens.loc[(slot,label),field] for slot in slots]
    full.to_csv(out/'private_all20_complete_source_review_v2.tsv',sep='\t');os.environ['MPLCONFIGDIR']=str(out/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'});y=np.arange(20);figures=[]
    def finish(fig,name):
        fig.tight_layout()
        for suffix in ['png','svg']:
            p=out/(name+'.'+suffix);fig.savefig(p,dpi=150,bbox_inches='tight');figures.append(dict(name=p.name,sha256=sha(p),bytes=p.stat().st_size))
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,12));ax.axvline(0,color='.55',lw=.8)
    for col,color,offset,label in [('unadjusted_effect','.72',-.2,'Unadjusted'),('clinical_only_effect','#2c7fb8',-.07,'Clinical background'),('clinical_RNA_no_factor_abundance_effect','#e0902e',.06,'Clinical + RNA background')]:ax.scatter(full[col],y+offset,s=16,c=color,label=label)
    for i,slot in enumerate(slots):
        r=full.loc[slot];color='#b2182b' if r.BH20<.05 else '#222222';ax.plot([r.supplementary_fixed_design_wild_t_CI_low,r.supplementary_fixed_design_wild_t_CI_high],[i+.2,i+.2],color=color,lw=1.2);ax.scatter(r.effect,i+.2,c=color,s=28)
    ax.set_yticks(y,slots);ax.invert_yaxis();ax.set_xlabel('Activity difference: pCR/MPR minus nonMPR (discovery SD)');ax.set_title('GSE243013 | all20 original main effects\nBars: supplementary pointwise fixed-design wild-t 95% CI\nOriginal case-bootstrap CIs unavailable: valid draws below predeclared 95% gate')
    handles,labels=ax.get_legend_handles_labels();handles+=[Line2D([0],[0],marker='o',color='#b2182b',label='Original main BH20 < .05'),Line2D([0],[0],marker='o',color='#222222',label='Original main BH20 >= .05')];ax.legend(handles=handles,loc='best',fontsize=8);finish(fig,'all20_main_effects_supplementary_fixed_design_CI')
    fig,ax=plt.subplots(figsize=(9,10));ax.axvline(0,c='.55',lw=.8);ax.scatter(full.secondary_Wilcoxon_rank_biserial,y,c=['#b2182b' if q<.05 else '.35' for q in full.secondary_Wilcoxon_BH20],s=28);ax.set_yticks(y,slots);ax.invert_yaxis();ax.set_xlim(-1.05,1.05);ax.set_xlabel('Unadjusted rank-biserial effect (positive = higher pCR/MPR)');ax.set_title('All20 secondary Wilcoxon comparisons\nUnadjusted; red = secondary BH20 < .05, gray = other');finish(fig,'all20_secondary_Wilcoxon_explicit_legend')
    matrix=np.column_stack([full.author_records/87,full.bootstrap_valid/2000,full.median_other_restart_cosine,1-full.median_heldout_reconstruction_error]);fig,ax=plt.subplots(figsize=(8,11));im=ax.imshow(matrix,aspect='auto',vmin=0,vmax=1,cmap='viridis');ax.set_yticks(y,slots);ax.set_xticks(range(4),['Heldout N /87','Case-bootstrap valid /2000','Other-start cosine','1 - heldout relative error'],rotation=25,ha='right')
    for i in range(20):
        for j in range(4):ax.text(j,i,'%.2f'%matrix[i,j],ha='center',va='center',color='#a50026' if j==1 and matrix[i,j]<.95 else 'white' if matrix[i,j]<.5 else 'black',fontsize=8)
    ax.set_title('All20 coverage and model diagnostics\nAll case-bootstrap valid fractions <.95: original CIs withheld\nAuthor records are not independently reviewed patients');fig.colorbar(im,ax=ax,shrink=.5);finish(fig,'all20_coverage_with_original_CI_failure_gate')
    fig,ax=plt.subplots(figsize=(8,5))
    for start,g in models.groupby('restart'):ax.plot(g['rank'],g.observed_relative_error,'o-',label='Fixed start%d'%start)
    ax.set_xticks(range(2,7));ax.set_xlabel('Fixed CP rank');ax.set_ylabel('Observed discovery relative reconstruction error');ax.set_title('All15 discovery fits | no clinical rank selection');ax.legend();finish(fig,'all15_fixed_model_reconstruction')
    fig,ax=plt.subplots(figsize=(11,12));ax.axvline(0,c='.55',lw=.8);ax.scatter(full.effect,y-.22,c='#222222',s=22,label='Original main effect')
    for label,color,offset,title in [('total_LR_and_capture','#e0902e',-.07,'+ total LR score / captured types'),('endpoint_RNA_quality','#2c7fb8',.08,'+ sender/receiver RNA quality'),('both','#238b45',.23,'+ both (HC3 t 95% CI)')]:
        ax.scatter(full[label+'__effect'],y+offset,c=color,s=19,label=title)
        if label=='both':
            for i,slot in enumerate(slots):ax.plot([full.loc[slot,label+'__HC3_t_CI_low'],full.loc[slot,label+'__HC3_t_CI_high']],[i+offset,i+offset],c=color,lw=.9)
    ax.set_yticks(y,[slot+' ['+full.loc[slot,'descriptive_redundancy_group']+']' for slot in slots]);ax.invert_yaxis();ax.set_xlabel('Activity difference (discovery SD)');ax.set_title('All20 POST-RESULT background diagnostics\nConditioning on total communication asks a relative-pattern question\nDescriptive groups deduplicate similar rank variants; no group-level FDR');ax.legend(loc='best',fontsize=8);finish(fig,'all20_post_result_global_and_endpoint_backgrounds')
    receipt=dict(status='PRIVATE_REVISED_ALL20_REVIEW_PRODUCED_ACTUAL_VISUAL_INSPECTION_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),original_audit_sha256=sha(root/'NSCLC_communication_tensor_family_audit_20261019_v2/private_terminal.json'),supplement_audit_sha256=sha(root/'NSCLC_tensor_complete_family_diagnostics_audit_20261019_v1/private_terminal.json'),figures=figures,complete_source_review_sha256=sha(out/'private_all20_complete_source_review_v2.tsv'),original_main_or_original_case_CI_changed=False,scope='All20; source-level only; descriptive group/no group FDR; supplementary intervals explicitly post-result/fixed-design/pointwise. Learned names and weights private.')
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');assert sum(p.stat().st_size for p in out.rglob('*') if p.is_file())<8*1024**2;print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
