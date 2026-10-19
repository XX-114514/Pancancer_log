"""Complete validation figures, including every held study and failed gate."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(a):
    r=a.root/'model_validation_v1'
    audit=a.root/'model_independent_audit_v2/independent_model_audit_technical_public.json'
    assert json.loads(audit.read_text())['producer_terminal_sha256']==hashlib.sha256((r/'model_validation_technical_public.json').read_bytes()).hexdigest()
    d=pd.read_csv(r/'private_all_model_metrics.tsv',sep='\t');c=pd.read_csv(r/'private_all_class_precision_recall.tsv',sep='\t')
    g=pd.read_csv(r/'private_all_frozen_engineering_gates.tsv',sep='\t');missing=pd.read_csv(r/'private_missing_spatial_marker_recovery_diagnostics.tsv',sep='\t')
    out=a.root/'validation_figures_v2';assert not out.exists();out.mkdir(mode=0o700)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    sources=sorted(d.holdout.unique());modes=['native','matched_depth_stress'];colors=['#1679a5','#df7b24'];names=[]
    def save(fig,name,caption):
        fig.text(.03,.01,caption,fontsize=8,ha='left',va='bottom');fig.subplots_adjust(bottom=.18)
        for ext in ['png','pdf','svg']:fig.savefig(out/(name+'.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
        fig.savefig(out/(name+'_preview.jpg'),dpi=130,bbox_inches='tight',facecolor='white');plt.close(fig);names.append(name)
    fig,ax=plt.subplots(figsize=(8,4.2));yy=np.arange(len(sources))
    for j,mode in enumerate(modes):
        vals=[d[(d.kind=='lineage')&(d.holdout==s)&(d['mode']==mode)].macroF1.iloc[0] for s in sources]
        ax.scatter(vals,yy+(j-.5)*.15,s=60,color=colors[j],label=['Native reference','Fixed target-depth stress'][j])
    ax.axvline(.8,color='#555',ls='--',lw=1);ax.set(yticks=yy,yticklabels=sources,xlim=(0,1),xlabel='Patient-balanced macro F1',title='CM02 panel bridge v20261020: all held-study lineage models');ax.legend(frameon=False,loc='lower right')
    save(fig,'All3_lineage_validation','Five reference-proposal classes; every study shown. Fixed engineering threshold = 0.80.\nReference annotation agreement is not independent biological cell-identity validation.')
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
    for ax,target in zip(axes,['regulatory','cytotoxic']):
        for j,mode in enumerate(modes):
            values=[d[(d.kind=='programme')&(d.target==target)&(d.holdout==s)&(d['mode']==mode)].rho.iloc[0] for s in sources]
            ax.scatter(values,yy+(j-.5)*.15,s=55,color=colors[j],label=['Native','Target-depth stress'][j])
        ax.axvline(.3,color='#777',ls=':',lw=1);ax.axvline(.5,color='#555',ls='--',lw=1);ax.set(xlim=(-1,1),yticks=yy,yticklabels=sources,xlabel='Case-weighted Spearman rho',title=target.capitalize()+' RNA programme')
    axes[1].legend(frameon=False,loc='lower right');fig.suptitle('All three held studies: both RNA programmes and both assay conditions')
    save(fig,'All3_two_programme_validation','No original 14 core marker or tested LR gene is a model feature or normalization-denominator feature.\nTwo continuous programmes can coexist. Per-study rho >=0.30 and median >=0.50 are engineering gates; no P values.')
    classes=['B_Plasma','Myeloid','NK','Nonimmune','T'];columns=[(s,m) for s in sources for m in modes]
    matrix=np.full((5,len(columns)),np.nan)
    for i,label in enumerate(classes):
        for j,(source,mode) in enumerate(columns):
            row=c[(c.lineage==label)&(c.holdout==source)&(c['mode']==mode)].iloc[0]
            if row.coverage_eligible:matrix[i,j]=row.precision
    fig,ax=plt.subplots(figsize=(10,4.6));im=ax.imshow(matrix,vmin=0,vmax=1,cmap='viridis',aspect='auto')
    ax.set(yticks=range(5),yticklabels=classes,xticks=range(6),xticklabels=[s+'\n'+('native' if m=='native' else 'stress') for s,m in columns],title='Every held-study class: patient-balanced precision')
    for i in range(5):
        for j in range(6):ax.text(j,i,'NA' if np.isnan(matrix[i,j]) else f'{matrix[i,j]:.2f}',ha='center',va='center',color='white' if np.nan_to_num(matrix[i,j])<.65 else '#111')
    fig.colorbar(im,ax=ax,label='Precision')
    save(fig,'All_class_precision_coverage','NA = fewer than100 cells or5 cases. Adopted-class precision gate=0.90, in addition to complete macro-F1 gate.\nNo failed class/source is removed and no cell-identity label is overwritten.')
    genes=['TNFRSF18','CTSW','GZMB','NKG7'];matrix=np.full((4,6),np.nan)
    for i,gene in enumerate(genes):
        for j,(source,mode) in enumerate(columns):matrix[i,j]=missing[(missing.missing_spatial_gene==gene)&(missing.holdout==source)&(missing['mode']==mode)].programme_prediction_vs_native_gene_rho.iloc[0]
    fig,ax=plt.subplots(figsize=(10,4.4));im=ax.imshow(matrix,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto');ax.set(yticks=range(4),yticklabels=genes,xticks=range(6),xticklabels=[s+'\n'+('native' if m=='native' else 'stress') for s,m in columns],title='All missing spatial-core genes: reference recovery diagnostics')
    for i in range(4):
        for j in range(6):ax.text(j,i,'NA' if np.isnan(matrix[i,j]) else f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=9)
    fig.colorbar(im,ax=ax,label='Case-weighted rank correlation')
    save(fig,'All_missing_marker_recovery','Programme prediction versus held-out native logCP10K gene expression; these genes are absent from the MERFISH panel.\nDiagnostic association is not recovered spatial gene counts or a gene-level imputation accuracy claim.')
    targets=classes+['regulatory','cytotoxic'];matrix=np.array([[bool(g[(g.target==t)&(g['mode']==m)].passed.iloc[0]) for m in modes] for t in targets])
    fig,ax=plt.subplots(figsize=(7,5));im=ax.imshow(matrix,cmap=matplotlib.colors.ListedColormap(['#d7dce0','#187c64']),vmin=0,vmax=1,aspect='auto');ax.set(yticks=range(7),yticklabels=targets,xticks=[0,1],xticklabels=['Native reference','Target-depth stress'],title='All frozen engineering gates: complete admission ledger')
    for i in range(7):
        for j in range(2):ax.text(j,i,'PASS' if matrix[i,j] else 'NOT PASSED',ha='center',va='center',color='white' if matrix[i,j] else '#333')
    save(fig,'Complete_engineering_admission','Both required T identity and both RNA programmes must pass before dual-arm spatial inference.\nA technical/numerical audit passing does not override failed scientific-use gates. No new spatial or clinical P values.')
    links=''.join('<section><h2>'+n.replace('_',' ')+'</h2><img style="width:100%" src="'+n+'.png"><p><a href="'+n+'.pdf">PDF</a> | <a href="'+n+'.svg">SVG</a></p></section>' for n in names)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>CM02 bridge validation</title><main style="max-width:1100px;margin:auto;font:16px system-ui"><h1>CM02 panel bridge: complete validation</h1>'+links+'</main>')
    receipt=dict(status='ALL5_COMPLETE_MODEL_VALIDATION_FIGURES_RENDERED_VISUAL_REVIEW_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),figures=names,
                 model_audit_sha256=hashlib.sha256(audit.read_bytes()).hexdigest(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),no_new_statistics=True)
    (out/'render_technical_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args())
