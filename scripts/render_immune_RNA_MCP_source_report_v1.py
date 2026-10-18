"""Render all source-level effects; no case records or learned genes/weights."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,16*1024**2)
    audit=root/'immune_RNA_MCP_stability_coordination_audit_20261018_v1/private_terminal.json'
    assert sha(audit)=='27efdf31dfee1a3937de6790a03947ff13ab7af594f8a74ddd3e5ffeeaf11a6a'
    source=root/'immune_RNA_MCP_source_coordination_20261018_v1'
    main_path=source/'private_all6_primary_source_tests.tsv'
    discovery_path=source/'private_all12_fitted_source_descriptions.tsv'
    main=pd.read_csv(main_path,sep='\t').sort_values('slot');discovery=pd.read_csv(discovery_path,sep='\t')
    assert len(main)==6 and len(discovery)==12
    output=root/'immune_RNA_MCP_source_report_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    (output/'plan_before_complete_source_render.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        main_source_table_sha256=sha(main_path),discovery_description_sha256=sha(discovery_path),
        numerical_acceptance_sha256=sha(audit),scope='All6 primary retained-source effects and all18 source descriptions, including unsupported uncertainty. No case/sample IDs, learned genes/coefficients or feature matrix.',
        visual_acceptance='PENDING actual pixel inspection; figures and scientific effects stay private, not public Git.'),indent=2)+'\n')
    short={'R_FOXP3_capture':'R','E_effector_capture':'E','Myeloid':'M'}
    labels=[f'{r.program}   {short[r.view_A]} - {short[r.view_B]}' for r in main.itertuples()]
    y=np.arange(6)[::-1];green='#008676';amber='#C88720';grey='#777777'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(13.5,7.6));fig.subplots_adjust(left=.19,right=.68,top=.80,bottom=.26)
    for row,yy in zip(main.itertuples(index=False),y):
        supported=row.BH6<.05 and (row.lower>0 or row.upper<0)
        color=green if supported else amber if row.BH6<.05 else grey
        ax.plot([row.lower,row.upper],[yy,yy],color=color,lw=2)
        ax.plot([row.lower,row.lower],[yy-.07,yy+.07],color=color,lw=1)
        ax.plot([row.upper,row.upper],[yy-.07,yy+.07],color=color,lw=1)
        ax.scatter(row.adjusted_rho,yy,s=65,color=color,zorder=4)
        ax.scatter(row.unadjusted_rho,yy,marker='D',facecolors='white',edgecolors='black',s=35,zorder=5)
        for xpos,value in [(1.08,str(row.complete_units)),(1.31,f'{row.BH6:.4f}'),(1.65,f'{row.BY6:.4f}')]:
            ax.text(xpos,yy,value,transform=ax.get_yaxis_transform(),va='center',ha='left',fontsize=10)
    for xpos,text in [(1.08,'N'),(1.31,'BH / 6'),(1.65,'BY / 6')]:
        ax.text(xpos,5.70,text,transform=ax.get_yaxis_transform(),fontweight='bold',ha='left',fontsize=11)
    ax.axvline(0,color='#999999',ls='--',lw=1)
    ax.set_xlim(-1,1);ax.set_ylim(-.6,5.55);ax.set_yticks(y,labels)
    ax.set_xlabel('Within-source rank correlation')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
    fig.suptitle('2026-10-18 RNA MCPs (v1, K2, source-adjusted)',y=.965,fontweight='bold',fontsize=16)
    fig.text(.19,.885,'GSE305078 / TNBC: 28 retained sample units; exact patient mapping pending',fontsize=12)
    handles=[Line2D([],[],marker='o',color=green,label='BH < .05 and bootstrap CI excludes zero'),
        Line2D([],[],marker='o',color=amber,label='BH < .05, bootstrap CI includes zero'),
        Line2D([],[],marker='D',color='black',mfc='white',ls='',label='Unadjusted Spearman (no CI shown)')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.18,.035),frameon=False,fontsize=10)
    fig.text(.19,.20,'R: FOXP3-supported T; E: effector-supported T; M: broad myeloid. RNA capture views, not formal subtypes.',fontsize=9.5)
    fig.text(.19,.16,'Adjusted: both view fractions/cell counts and native RNA/gene quality. 19,999 residual permutations; BH/BY over all six tests.',fontsize=9)
    fig.text(.19,.13,'Intervals: 2,000 planned sample-block bootstraps (1,999-2,000 valid). Evidence concerns coordinated RNA states, not LR binding.',fontsize=9)
    for ext in ['png','svg']:fig.savefig(output/('private_all6_retained_source_forest.'+ext),dpi=400,facecolor='white')
    plt.close(fig)
    order=['GSE286827','GSE294109','GSE305078'];values=np.full((6,3),np.nan)
    for i,row in enumerate(main.itertuples(index=False)):
        values[i,2]=row.unadjusted_rho
        for j,s in enumerate(order[:2]):
            d=discovery[discovery.slot.eq(row.slot)&discovery.dataset.eq(s)];assert len(d)==1
            values[i,j]=d.unadjusted_rho.iloc[0]
    fig,ax=plt.subplots(figsize=(10.5,7.3));fig.subplots_adjust(left=.22,right=.84,top=.75,bottom=.24)
    im=ax.imshow(values,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
    for i,j in np.ndindex(values.shape):
        ax.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',color='white' if abs(values[i,j])>.65 else 'black',fontsize=13)
    ax.set_yticks(np.arange(6),labels)
    ax.set_xticks(np.arange(3),['HNSCC / GSE286827\n28 discovery cases','ccRCC / GSE294109\n10 discovery cases','TNBC / GSE305078\n28 retained samples'],fontsize=9)
    ax.axvline(1.5,color='black',lw=2)
    fig.colorbar(im,ax=ax,fraction=.05,pad=.04,label='Unadjusted within-source Spearman')
    fig.suptitle('2026-10-18 RNA MCP source descriptions (v1)',y=.95,fontweight='bold',fontsize=15)
    fig.text(.22,.845,'Discovery descriptions helped fit the model; they are not independent replication.',fontsize=10)
    fig.text(.22,.135,'All six pairs shown. Retained-source adjusted effects and uncertainty are in the forest plot.',fontsize=9.5)
    fig.text(.22,.095,'Discovery units: reviewed cases, mixed sampling phases. Retained TNBC units: pretreatment samples, patient/treatment mapping pending.',fontsize=8.5)
    for ext in ['png','svg']:fig.savefig(output/('private_all18_source_descriptions.'+ext),dpi=400,facecolor='white')
    plt.close(fig)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(status='ALL_SOURCE_FIGURES_RENDERED_ACTUAL_VISUAL_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        biological_clinical_or_direct_interaction_validation=False),indent=2)+'\n')
    print('COMPLETE_SOURCE_FIGURES_RENDERED_VISUAL_ACCEPTANCE_PENDING')


if __name__=='__main__':main()
