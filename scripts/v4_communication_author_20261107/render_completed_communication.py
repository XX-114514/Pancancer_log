"""Reviewed static figures of aggregate real outputs; no patient/spot rows."""
import os, json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(os.environ['CCI_REVIEW_ROOT'])
tables=root/'completed_global_response_v1';out=root/'figures';out.mkdir(exist_ok=True)
proof=json.loads((tables/'complete.json').read_text());assert proof['status'].startswith('PASS_')
cms=['CM%02d'%i for i in range(1,15)]
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
g=pd.read_csv(tables/'global_all14_CM_by_CM_candidate_CCI.tsv',sep='\t')
matrix=g.pivot(index='source_CM',columns='target_CM',values='mean_candidate_CCI').reindex(index=cms,columns=cms)
fig,ax=plt.subplots(figsize=(10.3,9));fig.subplots_adjust(left=.12,right=.87,bottom=.23,top=.88)
cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#e4e4e4')
im=ax.imshow(np.ma.masked_invalid(matrix.to_numpy()),cmap=cmap,aspect='equal',vmin=0,vmax=float(np.nanmax(matrix)))
ax.set_xticks(range(14),cms,rotation=45,ha='right');ax.set_yticks(range(14),cms);ax.set_xlabel('CM of identity B in the ordered pair A|B');ax.set_ylabel('CM of identity A in the ordered pair A|B')
for i in range(14):
    for j in range(14):
        v=matrix.iloc[i,j];ax.text(j,i,'NA' if pd.isna(v) else f'{v:.0f}',ha='center',va='center',fontsize=7.5,color='black' if pd.isna(v) or v>np.nanmax(matrix)*.65 else 'white')
cax=fig.add_axes([.89,.26,.025,.48]);fig.colorbar(im,cax=cax,label='Mean candidate LR count per ordered identity pair')
fig.suptitle('2026-11-06 V4 | global CellPhoneDB candidate landscape\n38,711 balanced cells / 45 labels / 1,000 permutations',fontsize=13,y=.97)
fig.text(.08,.045,'All 14 original CMs retained; original empty CM01/02/05 = NA.\nMemberships overlap and include self pairs; matrix cells are not additive.\nA|B is an ordered identity pair, not a biological ligand-to-receptor direction.\nGlobal expression-specificity candidates, not High/Low effects, patient evidence or contact.',fontsize=9,ha='left')
fig.savefig(out/'global_CM_candidate_CCI.png',dpi=170);plt.close(fig)

response=pd.read_csv(tables/'response_full_family.tsv',sep='\t');elig=pd.read_csv(tables/'CM_member_signature_eligibility.tsv',sep='\t')
eligible=elig[elig.status.eq('READY_CELL_LEVEL_ONLY')].groupby('cm').size().reindex(cms,fill_value=0)
fig,axs=plt.subplots(2,1,figsize=(12,8.2),sharex=True);fig.subplots_adjust(left=.08,right=.98,bottom=.23,top=.85,hspace=.42)
labels=['Author-style unsigned signatures / annotation universe','Directional signatures / measured-gene universe (separate diagnostic)'];branches=['author_unsigned_annotation_universe','directional_measured_universe']
for ax,branch,title in zip(axs,branches,labels):
    a=response[response.branch.eq(branch)];one=a[a.q_within_CM_member.lt(.05)].groupby('cm').size().reindex(cms,fill_value=0);full=a[a.q_full_family.lt(.05)].groupby('cm').size().reindex(cms,fill_value=0);x=np.arange(14)
    ax.bar(x-.18,one,width=.36,color='#73889d',label='q < .05 within CM-member')
    ax.bar(x+.18,full,width=.36,color='#16798b',label='BH < .05 across complete branch')
    for i,cm in enumerate(cms):
        if eligible[cm]==0:ax.text(i,.5,'NA',ha='center',transform=ax.get_xaxis_transform(),color='#555555',fontsize=8)
        elif full[cm]>0:ax.text(i+.18,full[cm]+.7,str(full[cm]),ha='center',fontsize=8)
    ax.set_ylim(0,max(5,float(one.max())*1.18));ax.set_ylabel('CM-member / cytokine / direction\nassociations (not unique cytokines)',fontsize=9);ax.set_title(title,fontsize=11,loc='left')
handles,legendlabels=axs[0].get_legend_handles_labels();fig.legend(handles,legendlabels,loc='upper right',bbox_to_anchor=(.98,.94),ncol=2,frameon=False,fontsize=9)
axs[-1].set_xticks(range(14),[f'{cm}\n({eligible[cm]})' for cm in cms],fontsize=9);axs[-1].set_xlabel('Original V4 CM (number of eligible immune member comparisons)',fontsize=9)
fig.suptitle('2026-11-06 V4 | cytokine-response sensitivity to analysis definition',fontsize=14,y=.97)
fig.text(.08,.055,'CM03 has no author-style within-member q < .05 response; diagnostic positives do not replace this result.\nUp/down denote signature-gene directions in cytokine perturbation; down-signature matches do not imply lower cytokine activity.\nSame-label pooled cells; source/patient confounding and reference-identity limits remain. No intervention was performed.',fontsize=9)
fig.savefig(out/'CM_response_branch_comparison.png',dpi=170);plt.close(fig)
print(json.dumps({'status':'REAL_AGGREGATE_FIGURES_GENERATED_PENDING_VISUAL_REVIEW','n_figures':2}))
