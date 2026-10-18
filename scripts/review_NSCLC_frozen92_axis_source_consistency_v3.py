"""All92 source-only figures; no learned molecular names in display tables."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);out=root/'NSCLC_frozen92_axis_source_review_20261019_v3';assert not out.exists();out.mkdir(mode=0o700);src=root/'NSCLC_frozen92_axis_source_consistency_20261019_v1';proof=root/'NSCLC_frozen92_axis_source_result_audit_20261019_v1/private_terminal.json';assert sha(proof)=='79e8aa29aec12db05d437f6e4ebc24dec7e97772bf62ca594434efe36522b9b7';p=json.loads(proof.read_text());assert p['producer_terminal_sha256']==sha(src/'private_terminal.json');t=json.loads((src/'private_terminal.json').read_text())
    for name,d in t['files'].items():assert sha(src/name)==d['sha256']
    r=pd.read_csv(src/'private_all92_source_consistency_results.tsv',sep='\t');s=pd.read_csv(src/'private_all276_fixed_sensitivity_effects.tsv',sep='\t');h=pd.read_csv(src/'private_all184_histology_descriptors.tsv',sep='\t').set_index(['axis_id','cancer_type']);reg=pd.read_csv(root/'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_all92_frozen_axis_registry.tsv',sep='\t');assert len(r)==92 and r.axis_id.tolist()==reg.axis_id.tolist();background=s[s.sensitivity.eq('main_plus_total_communication')].set_index('axis_id');z=r.merge(reg[['axis_id','sender_original_author_type','receiver_original_author_type']],on='axis_id',validate='one_to_one');z['evidence_class']=np.select([~z.status.eq('ESTIMABLE'),z.BY92.lt(.05),z.BH92.lt(.05),z.wild_P.lt(.05)],['COVERAGE_OR_DESIGN_INSUFFICIENT','BH_AND_BY_SOURCE_CANDIDATE','BH_ONLY_SOURCE_CANDIDATE','NOMINAL_ONLY_NO_MULTIPLICITY_SUPPORT'],default='NO_CORRECTED_SUPPORT');z.to_csv(out/'private_all92_source_only_review.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False});panels=[r.iloc[i:i+23] for i in range(0,92,23)]
    def canvas(title,subtitle):
        fig,axes=plt.subplots(2,2,figsize=(15,13));fig.suptitle(title+'\n'+subtitle,fontsize=12,y=.99);fig.subplots_adjust(top=.91,bottom=.09,hspace=.25,wspace=.40);return fig,axes.ravel()
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=150);fig.savefig(out/(name+'.pdf'));plt.close(fig)
    def axis_rows(ax,part):
        y=np.arange(len(part));ax.set_yticks(y);ax.set_yticklabels(['']*len(part))
        for label in list(ax.texts):
            if label.get_gid()=='manual_axis_id':label.remove()
        for position,identifier in zip(y,part.axis_id.astype(str).tolist()):ax.text(-.025,position,identifier,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=9,clip_on=False,gid='manual_axis_id')
        ax.set_ylim(len(part)-.4,-.8);ax.grid(axis='x',alpha=.15);return y
    fig,axes=canvas('GSE243013: complete frozen 92-axis clinical family','Separate-label-set within-cohort check; 129 author records, not verified independent patients; pointwise WCU-B 95% intervals')
    limits=[float(r.CI_low.min())-.12,float(r.CI_high.max())+.12]
    for ax,part in zip(axes,panels):
        axis_rows(ax,part);ax.axvline(0,color='black',lw=.7)
        for y,row in enumerate(part.itertuples()):
            colour='#a12b37' if row.BY92<.05 else '#266ab2' if row.BH92<.05 else '#777777'
            if np.isfinite(row.effect):ax.errorbar(row.effect,y,xerr=[[max(0,row.effect-row.CI_low)],[max(0,row.CI_high-row.effect)]],fmt='o',color=colour,ms=3,capsize=2,lw=.8)
            else:ax.text(.04,y,'Not estimable',transform=ax.get_yaxis_transform(),fontsize=8,color='#777777')
            ax.text(1.02,y,'n=%d  q=%.3g'%(row.author_records,row.BH92),transform=ax.get_yaxis_transform(),va='center',fontsize=7)
        ax.set_xlim(*limits);ax.set_xlabel('Adjusted difference: pCR/MPR minus non-MPR (axis SD)')
    fig.legend(handles=[Line2D([0],[0],color='#266ab2',marker='o',label='Primary BH92 < .05'),Line2D([0],[0],color='#a12b37',marker='o',label='Primary BY92 < .05'),Line2D([0],[0],color='#777777',marker='o',label='No corrected support')],loc='lower center',ncol=3,frameon=False);save(fig,'all92_primary_source_effects')
    fig,axes=canvas('GSE243013: complete per-axis clinical and expression coverage','Bars: complete clinical records. Black tick: both endpoints >=20 cells and every complex gene >=10% detection; detection never filters main records.')
    for ax,part in zip(axes,panels):
        y=axis_rows(ax,part);ax.barh(y,part.nonMPR,color='#647b98',label='non-MPR');ax.barh(y,part.pCR_or_MPR,left=part.nonMPR,color='#d6aa69',label='pCR/MPR');ax.plot(part.strict_detection_supported,y,'k|',ms=7,label='Strict expression support');ax.axvline(50,color='#777777',ls=':',lw=.8);ax.set_xlim(0,130);ax.set_xlabel('Author record count; not an independent-patient count');axis_rows(ax,part)
    fig.legend(handles=[Patch(color='#647b98',label='non-MPR clinical records'),Patch(color='#d6aa69',label='pCR/MPR clinical records'),Line2D([0],[0],color='black',marker='|',linestyle='',label='Expression-support records')],loc='lower center',ncol=3,frameon=False);save(fig,'all92_clinical_expression_coverage')
    fig,axes=canvas('GSE243013: additional total-communication background','Blue: primary model, WCU-B CI. Orange: total LR / capture additionally controlled, CV3 t-reference CI. Orange has no new multiplicity-adjusted test.')
    lower=min(float(r.CI_low.min()),float(background.CV3_t_CI_low.min()))-.12;upper=max(float(r.CI_high.max()),float(background.CV3_t_CI_high.max()))+.12
    for ax,part in zip(axes,panels):
        axis_rows(ax,part);ax.axvline(0,color='black',lw=.7)
        for y,row in enumerate(part.itertuples()):
            for offset,value,lo,hi,colour,marker in [(-.12,row.effect,row.CI_low,row.CI_high,'#266ab2','o'),(.12,background.loc[row.axis_id].effect,background.loc[row.axis_id].CV3_t_CI_low,background.loc[row.axis_id].CV3_t_CI_high,'#c67b36','^')]:
                if np.isfinite(value):ax.errorbar(value,y+offset,xerr=[[max(0,value-lo)],[max(0,hi-value)]],fmt=marker,color=colour,ms=3,lw=.6,capsize=1)
        ax.set_xlim(lower,upper);ax.set_xlabel('Adjusted difference (frozen axis SD)');axis_rows(ax,part)
    fig.legend(handles=[Line2D([0],[0],color='#266ab2',marker='o',label='Primary; WCU-B interval'),Line2D([0],[0],color='#c67b36',marker='^',label='Additional background; CV3 t interval')],loc='lower center',ncol=2,frameon=False);save(fig,'all92_main_and_total_background')
    fig,axes=canvas('GSE243013: histology descriptions for all 92 axes','Unadjusted differences only; no histology P values or independent-cancer replication claim. Gray: fewer than10 records in either response class.')
    for ax,part in zip(axes,panels):
        values=np.full((len(part),2),np.nan)
        for y,row in enumerate(part.itertuples()):
            for c,cancer in enumerate(['LUAD','LUSC']):
                rr=h.loc[(row.axis_id,cancer)]
                if min(rr.nonMPR,rr.pCR_or_MPR)>=10:values[y,c]=rr.unadjusted_effect
        cmap=plt.get_cmap('coolwarm').copy();cmap.set_bad('#eeeeee');im=ax.imshow(values,aspect='auto',cmap=cmap,vmin=-1.5,vmax=1.5);ax.set_xticks([0,1]);ax.set_xticklabels(['LUAD','LUSC']);ax.set_yticks(range(len(part)));ax.set_yticklabels(part.axis_id)
        for y,row in enumerate(part.itertuples()):
            for c,cancer in enumerate(['LUAD','LUSC']):
                rr=h.loc[(row.axis_id,cancer)];label='%.2f'%values[y,c] if np.isfinite(values[y,c]) else 'low coverage';ax.text(c,y,label+'  (%d/%d)'%(rr.nonMPR,rr.pCR_or_MPR),ha='center',va='center',fontsize=7)
    fig.text(.5,.045,'Cell labels show unadjusted difference and (non-MPR / pCR-or-MPR counts). Color clipped at +/-1.5 SD.',ha='center',fontsize=9);save(fig,'all92_histology_descriptions')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<4*1024**2;receipt=dict(status='COMPLETE_ALL92_SOURCE_ONLY_FOUR_PNG_PDF_SETS_ACTUAL_VISUAL_INSPECTION_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),result_audit_sha256=sha(proof),resource_guard=guard,files=files,learned_molecular_names_exported=False,individual_records_exported=False,scientific_outputs_public=False)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
