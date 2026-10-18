"""Approved source-only follow-up forest; no learned node coefficients."""
import json
import os
from datetime import datetime, timezone
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
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    audit=root/'immune_RNA_MCP_original_node_decomposition_audit_20261018_v1/private_terminal.json'
    assert sha(audit)=='f109cecf33be6ad8a937bce8c50d567aa1426f1a1e1ea05a0c20335ce4e123d7'
    table=root/'immune_RNA_MCP_original_node_decomposition_20261018_v1/private_all6_composition_followup_tests.tsv'
    assert sha(table)=='c1913e150ff748dda1b7ebde2468aa50113359412988b2faa348fb888c2c1cdd'
    d=pd.read_csv(table,sep='\t').sort_values('slot');assert len(d)==6
    out=root/'immune_RNA_MCP_node_composition_report_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_source_forest.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        source_table_sha256=sha(table),numeric_acceptance_sha256=sha(audit),
        scope='All six additional source-level composition follow-up effects; no case records, genes, learned node coefficients or matrix. Actual pixel review pending.'),indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(13.5,7.6));fig.subplots_adjust(left=.19,right=.68,top=.80,bottom=.37)
    short={'R_FOXP3_capture':'R','E_effector_capture':'E','Myeloid':'M'}
    labels=[f'{r.program}   {short[r.view_A]} - {short[r.view_B]}' for r in d.itertuples()]
    for r,y in zip(d.itertuples(index=False),np.arange(6)[::-1]):
        color='#008676' if r.BH6<.05 and (r.lower>0 or r.upper<0) else '#C88720' if r.BH6<.05 else '#777777'
        ax.errorbar(r.adjusted_for_frozen_node_composition,y,xerr=[[r.adjusted_for_frozen_node_composition-r.lower],[r.upper-r.adjusted_for_frozen_node_composition]],fmt='o',color=color,capsize=3,lw=2,ms=7)
        ax.scatter(r.original_adjusted_same_cases,y-.13,marker='D',facecolors='white',edgecolors='black',s=35,zorder=5)
        for xpos,value in [(1.08,str(r.complete_units)),(1.28,f'{r.BH6:.4f}'),(1.48,f'{r.BY6:.4f}')]:
            ax.text(xpos,y,value,transform=ax.get_yaxis_transform(),va='center',ha='left',fontsize=10)
    for xpos,text in [(1.08,'N'),(1.28,'BH / 6'),(1.48,'BY / 6')]:
        ax.text(xpos,5.7,text,transform=ax.get_yaxis_transform(),fontweight='bold',fontsize=11)
    ax.axvline(0,color='#999999',ls='--',lw=1);ax.set_xlim(-1,1);ax.set_ylim(-.6,5.55)
    ax.set_yticks(np.arange(6)[::-1],labels);ax.set_xlabel('Within-source partial rank correlation')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
    fig.suptitle('2026-10-18 RNA MCPs: original-node composition follow-up (v1)',y=.965,fontweight='bold',fontsize=15)
    fig.text(.19,.885,'GSE305078 / TNBC: source retained from MCP training; exact patient mapping pending',fontsize=11.5)
    fig.text(.19,.24,'R: FOXP3-supported T; E: effector-supported T; M: broad myeloid. Capture views, not formal subtypes.',fontsize=9.5)
    fig.text(.19,.205,'Original abundance/RNA controls + both node-mixture predictions fixed from 38 discovery cases.',fontsize=9)
    fig.text(.19,.17,'N = 25 for R pairs (3 samples have unseen training nodes); N = 28 for E-M. Missing is not zero.',fontsize=9)
    fig.text(.19,.135,'19,999 permutations; 2,000 planned bootstraps; separate BH/BY6 follow-up family, same retained source.',fontsize=9)
    fig.text(.19,.10,'Intervals condition on fixed programmes/node effects; these are not ligand-receptor or clinical tests.',fontsize=9)
    handles=[Line2D([],[],marker='o',color='#008676',label='BH < .05 and bootstrap CI excludes zero'),
        Line2D([],[],marker='o',color='#C88720',label='BH < .05, bootstrap CI includes zero'),
        Line2D([],[],marker='o',color='#777777',label='Follow-up evidence insufficient'),
        Line2D([],[],marker='D',color='black',mfc='white',ls='',label='Original eight controls, same samples (no CI)')]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.18,.012),frameon=False,fontsize=9.5,ncol=2)
    for ext in ['png','svg']:fig.savefig(out/('private_all6_node_adjusted_forest.'+ext),dpi=400,facecolor='white')
    plt.close(fig)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ALL6_FOLLOWUP_FOREST_RENDERED_ACTUAL_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files),indent=2)+'\n')
    print('ALL6_FOLLOWUP_SOURCE_FOREST_RENDERED_ACTUAL_REVIEW_PENDING')


if __name__=='__main__':main()
