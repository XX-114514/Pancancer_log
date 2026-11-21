#!/usr/bin/env python3
"""Recreate manuscript panels from supplied plot arrays and aggregate tables.

No expression matrix, patient row, learned model coefficient or fitted model is
required. Statistical values are read from frozen external tables.
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

from umap_labels import add_umap_labels
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#147D92';PURPLE='#8F849D';GREY='#9AA8B3';ORANGE='#B8664D'
def panel(ax,k,title):ax.set_title(title,loc='left',pad=12);ax.text(-.1,1.09,k,transform=ax.transAxes,fontsize=13,fontweight='bold')
def save(fig,out,name):
 for ext in ['png','pdf','svg']:fig.savefig(out/(name+'.'+ext),dpi=330,bbox_inches='tight',pad_inches=.14)
 plt.close(fig);print('RENDERED',name,flush=True)
def run(inputs,out):
 out.mkdir(parents=True,exist_ok=False)
 j=json.loads((inputs/'integrated_evidence/integrated_aggregate_evidence_v3.json').read_text()) if (inputs/'integrated_evidence/integrated_aggregate_evidence_v3.json').exists() else json.loads((inputs/'integrated_evidence/integrated_aggregate_evidence_v2.json').read_text())
 cfg=json.loads((inputs/'figure_configuration.json').read_text());ev=inputs/'integrated_evidence';immune=inputs/'immune'
 core=cfg['immune_core'];means=pd.read_csv(immune/'T_node_state_means.tsv',sep='\t',index_col=0);states=means.columns.tolist();z=(means-means.mean())/means.std(ddof=0)
 bal=pd.DataFrame(j['balance']['component_comparisons']);b=bal[(bal.reference_module=='CM02')&(bal['mode']=='equal_source_weight')]
 fig=plt.figure(figsize=(9.2,7.8));gs=fig.add_gridspec(2,2,hspace=.62,wspace=.65)
 ax=fig.add_subplot(gs[0,0]);ax.axis('off');panel(ax,'a','A seven-member T-cell combination')
 for i,n in enumerate(core):ax.text(.02,.9-i*.12,n,fontsize=10,color=BLUE)
 ax=fig.add_subplot(gs[0,1]);a=z.reindex(core);assert not a.isna().any().any();im=ax.imshow(a,aspect='auto',vmin=-2,vmax=2,cmap='PuOr');ax.set_yticks(range(len(core)),[n.replace('T_','') for n in core],fontsize=8);ax.set_xticks(range(5),[s.replace('state_score_','').replace('_',' ') for s in states],rotation=45,ha='right',fontsize=7);fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label='Across-node state z score');panel(ax,'b','The combination spans distinct RNA states')
 ax=fig.add_subplot(gs[1,0]);ax.scatter(np.arange(len(b)),b.effective_A_cosine,s=40,color=BLUE);ax.set_ylim(.9,1.01);ax.set_xticks(range(len(b)),['Seed '+str(i+1) for i in range(len(b))]);ax.set_ylabel('Matched effective-basis cosine');panel(ax,'c','The core persists under source balancing');ax.text(0,-.3,'Equal total objective weight per source; fixed-rank sensitivity.',transform=ax.transAxes,fontsize=8)
 ax=fig.add_subplot(gs[1,1]);m=j['MCP_TNBC_heldout'];rho=m['effect_after_original8_controls_and_both_fixed_node_mixture_scores'];lo,hi=m['conditional95CI'];ax.plot([lo,hi],[0,0],color=BLUE,lw=2);ax.scatter(rho,0,color=BLUE,s=45);ax.axvline(0,color=GREY,ls='--');ax.set_xlim(-.05,1);ax.set_ylim(-.8,.8);ax.set_yticks([0],['Effector T - myeloid'],fontsize=8);ax.set_xlabel('Adjusted RNA-program correlation');panel(ax,'d','State coordination extends to myeloid cells');ax.text(0,-.29,f"{m['units']} TNBC sample units | BH6={m['followup_BH6']:.4f}\nFixed node-mixture and RNA-quality controls",transform=ax.transAxes,fontsize=8)
 save(fig,out,'Figure3_immune_state_combination')
 clinical=pd.read_csv(ev/'immune_clinical_all_model_slots.tsv',sep='\t');c=clinical[clinical.axis=='pCM02_regulatory_cytotoxic'];fig=plt.figure(figsize=(9,4.5));gs=fig.add_gridspec(1,3,width_ratios=[1,1.2,1.3],wspace=.75)
 ax=fig.add_subplot(gs[0,0]);ax.axis('off');panel(ax,'a','A fixed two-arm RNA readout')
 for yy,label in [(.78,'Regulatory gene panel'),(.53,'Cytotoxic gene panel')]:
  ax.add_patch(FancyBboxPatch((.03,yy-.06),.94,.15,boxstyle='round,pad=.02',fc='#EDF3F7',ec=GREY));ax.text(.5,yy,label,ha='center',va='center',fontsize=8)
 ax.text(.5,.26,'Joint score = min(two arm means)\nPan-T-adjusted | per SD',ha='center',fontsize=8,linespacing=1.6)
 ax=fig.add_subplot(gs[0,1]);q=c[(c.cohort=='IMvigor210')&(c.family=='ICI_response')].iloc[0];ax.plot([q.CI_low,q.CI_high],[0,0],color=BLUE,lw=2);ax.scatter(q.effect,0,color=BLUE,s=40);ax.axvline(1,color=GREY,ls='--');ax.set_xscale('log');ax.set_xlim(.8,9);ax.xaxis.set_major_locator(FixedLocator([1,2,4,8]));ax.xaxis.set_major_formatter(FixedFormatter(['1','2','4','8']));ax.xaxis.set_minor_locator(NullLocator());ax.set_ylim(-.8,.8);ax.set_yticks([0],['IMvigor210'],fontsize=8);ax.set_xlabel('Response OR (95% interval)');panel(ax,'b','Association with ICI response');ax.text(0,-.22,f"n={int(q.n)} | responders={int(q.events)}\nBH={q.BH_q:.4g}",transform=ax.transAxes,fontsize=8)
 ax=fig.add_subplot(gs[0,2]);q=c[(c.cohort=='SU2C_MARK')&(c.family=='SU2C_survival')].set_index('endpoint').loc[['OS','PFS']]
 for i,r in enumerate(q.itertuples()):ax.plot([r.CI_low,r.CI_high],[i,i],color=BLUE,lw=2);ax.scatter(r.effect,i,color=BLUE,s=40);ax.text(.04,i-.22,f'n={int(r.n)}; events={int(r.events)}; BH={r.BH_q:.3g}',transform=ax.get_yaxis_transform(),fontsize=6.8)
 ax.axvline(1,color=GREY,ls='--');ax.set_xscale('log');ax.set_xlim(.25,1.2);ax.xaxis.set_major_locator(FixedLocator([.25,.5,1]));ax.xaxis.set_major_formatter(FixedFormatter(['0.25','0.5','1']));ax.xaxis.set_minor_locator(NullLocator());ax.set_ylim(1.6,-.7);ax.set_yticks([0,1],['SU2C OS','SU2C PFS']);ax.set_xlabel('Hazard ratio (95% interval)');panel(ax,'c','Association with survival')
 save(fig,out,'Figure4_immune_clinical_readouts')
 fig,ax=plt.subplots(figsize=(8.8,4.6));cv=pd.DataFrame(j['CV_contrasts']);cv=cv[(cv.module_proxy=='CM02_proxy')&cv.pooled_paired_logloss_improvement.notna()].reset_index(drop=True)
 for i,r in enumerate(cv.itertuples()):ax.plot([r.conditional_CI_low,r.conditional_CI_high],[i,i],color=PURPLE,lw=2);ax.scatter(r.pooled_paired_logloss_improvement,i,color=PURPLE,s=35)
 ax.set_yticks(range(len(cv)),[r.cohort+' | '+r.candidate+' vs '+r.reference for r in cv.itertuples()],fontsize=8);ax.invert_yaxis();ax.axvline(0,color=GREY,ls='--');ax.set_xlabel('OOF log-loss improvement (conditional 95% interval)');ax.set_title('Added prediction information for the fixed immune readout',loc='left');save(fig,out,'FigureS09_prediction_information')
 r=pd.read_csv(inputs/'R4_figure_tables/Figure7_all15_immune_associations.tsv',sep='\t');q=r[(r.branch=='original_network_members')&(r.model=='purity_age_sex')&(r.policy=='LF_single')];names=[f'CM{i:02d}' for i in range(1,16)];can=sorted(q.cancer.unique());mat=q.pivot(index='cm',columns='cancer',values='partial_rho').reindex(index=names,columns=can);qq=q.pivot(index='cm',columns='cancer',values='q_joint_with_unadjusted').reindex_like(mat)
 fig=plt.figure(figsize=(10,7.4));gs=fig.add_gridspec(2,1,height_ratios=[1.8,1],hspace=.72);ax=fig.add_subplot(gs[0]);cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#D9DDE1');im=ax.imshow(mat,aspect='auto',cmap=cm,vmin=-1,vmax=1);ax.set_xticks(range(len(can)),can,rotation=90,fontsize=7);ax.set_yticks(range(15),names,fontsize=7)
 for i,k in zip(*np.where(qq.to_numpy()<.05)):ax.text(k,i,'*',ha='center',va='center',fontsize=8)
 fig.colorbar(im,ax=ax,fraction=.018,pad=.02,label='Adjusted rank correlation');panel(ax,'a','Cancer context differentiates multicellular RNA readouts')
 ax=fig.add_subplot(gs[1]);v=mat.loc['CM12'];sig=qq.loc['CM12']<.05;ax.bar(range(len(can)),v,color=[BLUE if s else '#C6CED4' for s in sig]);ax.set_xticks(range(len(can)),can,rotation=90,fontsize=7);ax.axhline(0,color=GREY,lw=.7);ax.set_ylabel('Adjusted rank correlation');panel(ax,'b','The vascular-associated readout has cancer-specific directions')
 save(fig,out,'Figure7_cancer_context')
 # Recreate both immune archive displays from identifier-free plotting inputs.
 with np.load(immune/'immune_embedding.npz') as zz:xy=zz['xy'];major=zz['major'];ds=zz['dataset'];majors=zz['majors'].tolist();datasets=zz['datasets'].tolist()
 palette=plt.get_cmap('tab20')(np.linspace(0,1,len(majors)));rng=np.random.default_rng(20261014);order=rng.permutation(len(xy));fig,axes=plt.subplots(2,2,figsize=(9,8));fig.subplots_adjust(hspace=.5,wspace=.35,bottom=.16)
 for ax,codes,pal,title in [(axes[0,0],major,palette,'Major cell identities'),(axes[0,1],ds,plt.get_cmap('tab10')(np.arange(len(datasets))),'Source composition')]:
  ax.scatter(xy[order,0],xy[order,1],c=pal[codes[order]],s=.07,linewidths=0,rasterized=True);ax.set_title(title,loc='left');ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('Frozen UMAP 1');ax.set_ylabel('Frozen UMAP 2')
 legend_only = {'Unknown', 'Ambiguous', 'Unresolved', 'T_NK_boundary'}
 add_umap_labels(axes[0,0], xy, major, [
  None if name in legend_only else name.replace('_', ' ') for name in majors
 ], fontsize=6)
 s=pd.read_csv(immune/'source_coverage.tsv',sep='\t');display_source=s.dataset.map(lambda x:'CXG CAR collection' if str(x).startswith('CXG_') else x);axes[1,0].barh(range(len(s)),s.cells/1000,color=BLUE);axes[1,0].set_yticks(range(len(s)),display_source,fontsize=7);axes[1,0].set_xlabel('Archived cells (thousands)');axes[1,1].barh(range(len(s)),s.archived_native_endpoint_eligible/s.cells*100,color=BLUE);axes[1,1].set_yticks(range(len(s)),display_source,fontsize=7);axes[1,1].set_xlabel('Archived native endpoint eligibility (%)')
 fig.legend(handles=[plt.Line2D([],[],ls='',marker='o',color=palette[i],label=m,markersize=4) for i,m in enumerate(majors)],loc='lower center',ncol=6,frameon=False,fontsize=7);save(fig,out,'FigureS07_immune_atlas')
 t=pd.read_csv(immune/'T_state_display.tsv.gz',sep='\t');lims=json.loads((immune/'state_display_scale.json').read_text());fig,axes=plt.subplots(2,3,figsize=(10,7.5));fig.subplots_adjust(hspace=.5,wspace=.65,bottom=.18)
 for ax,key in zip(axes.flat,states):
  im=ax.scatter(t.x,t.y,c=t[key],s=.2,cmap='cividis',vmin=0,vmax=lims[key],rasterized=True,linewidths=0);ax.set_title(key.replace('state_score_','').replace('_',' '),loc='left');ax.set_xlabel('Frozen UMAP 1');ax.set_ylabel('Frozen UMAP 2');fig.colorbar(im,ax=ax,fraction=.04,pad=.02)
 ax=axes[1,2];im=ax.imshow(z,vmin=-2,vmax=2,cmap='PuOr',aspect='auto');ax.set_xticks(range(5),[s.replace('state_score_','').replace('_',' ') for s in states],rotation=90,fontsize=7);ax.set_yticks(range(len(z)),[x.replace('T_','') for x in z.index],fontsize=7);ax.set_title('T-node state means',loc='left');fig.colorbar(im,ax=ax,fraction=.04,pad=.02);save(fig,out,'FigureS08_T_state_dimensions')
 (out/'render_receipt.json').write_text(json.dumps({'status':'RENDERED','new_statistical_tests':False,'figures':[p.stem for p in out.glob('*.png')]},indent=2)+'\n')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--inputs',type=Path,required=True);a.add_argument('--output',type=Path,required=True);v=a.parse_args();run(v.inputs,v.output)
