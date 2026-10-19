"""Source-only all-slot plot; no identities, learned genes or patient rows."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def groups(frame):
    out=np.full(len(frame),'Non-T/NK',object)
    out[frame.major_cell_type.eq('T/NK cell')]='Other T/NK'
    out[frame.sub_cell_type.str.startswith('NK_')]='Author NK'
    out[frame.sub_cell_type.eq('T_gdT_TRDV1')]='Author Vdelta1'
    out[frame.sub_cell_type.eq('T_gdT_TRDV2')]='Author Vdelta2'
    return out


def main():
    root=Path(os.environ['PANCANCER_NORMALIZED_GAMMA_SOURCE_REVIEW'])
    receipt=json.loads((root/'private_selected_summary_copy_receipt.json').read_text())
    assert receipt['producer_terminal_sha256']=='ef64454f3d5e9867f904c5cf50599a421670910a6bf6c39e1cb904a6abca5209'
    assert receipt['audit_terminal_sha256']=='666e90720403e7b764c3e3dd21cea7ac76eb54b97fbcd6d942412fcbf4d2fa0f'
    for name,d in receipt['selected_files'].items():assert sha(root/name)==d['sha256']
    f=pd.read_csv(root/'audit/private_complete51_by4_author_reference_counts.tsv',sep='\t')
    joint=pd.read_csv(root/'audit/private_complete2448_author_AB_NK_rule_model_descriptions.tsv.gz',sep='\t')
    assert len(f)==204 and len(joint)==2448
    f['reference_group']=groups(f);joint['reference_group']=groups(joint)
    summary=f.groupby(['source','model','reference_group'],sort=True)[['pilot_cells','pure_cells']].sum().reset_index()
    conflict=joint.groupby(['source','model','AB_annotation_status','NK_model_pure','any_original_raw_rule'],sort=True)[['pilot_cells','model_pure_cells']].sum().reset_index()
    assert summary.groupby('model').pilot_cells.sum().eq(20000).all()
    assert conflict.groupby('model').pilot_cells.sum().eq(20000).all()
    summary.to_csv(root/'private_complete_source_reference_counts.tsv',sep='\t',index=False)
    conflict.to_csv(root/'private_complete_source_annotation_NK_rule_overlaps.tsv',sep='\t',index=False)
    models=['HiTME_full_symbol_normalized','HiTME_native_subset','PBMC_exact','NK_HiTME_exact']
    labels=['HiTME full: public-symbol normalized','HiTME original measured subset','PBMC original model','NK control (different target)']
    categories=['Author Vdelta1','Author Vdelta2','Author NK','Other T/NK','Non-T/NK']
    counts=summary.pivot(index='model',columns='reference_group',values='pure_cells').reindex(index=models,columns=categories).to_numpy(float)
    totals=summary.pivot(index='model',columns='reference_group',values='pilot_cells').reindex(index=models,columns=categories).to_numpy(float)
    assert np.array_equal(totals,np.repeat(totals[:1],4,axis=0))
    coverage=np.divide(counts,totals,out=np.full_like(counts,np.nan),where=totals>0)
    composition=np.divide(counts,counts.sum(1)[:,None],out=np.full_like(counts,np.nan),where=counts.sum(1)[:,None]>0)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,2,figsize=(16,6.8))
    for ax,values,title in zip(axs,[coverage,composition],['Coverage of each author reference group','Composition among model calls']):
        im=ax.imshow(values*100,vmin=0,vmax=100,cmap='Blues',aspect='auto')
        ax.set_xticks(range(5),[c+'\n(n=%d)'%n for c,n in zip(categories,totals[0])],rotation=23,ha='right')
        ax.set_yticks(range(4),labels)
        ax.set_title(title,pad=12)
        for i in range(4):
            for j in range(5):
                v=values[i,j]
                ax.text(j,i,('%.1f%%\nn=%d'%(100*v,counts[i,j])) if np.isfinite(v) else 'undefined',ha='center',va='center',fontsize=9,color='white' if np.isfinite(v) and v>.6 else 'black')
        fig.colorbar(im,ax=ax,shrink=.7,label='%')
    fig.suptitle('GSE243013 | same fixed 20,000-cell pilot | all four executed model slots',fontsize=14,y=.98)
    fig.text(.025,.025,'Author reference labels are not TCR ground truth. Coverage is not validated sensitivity.\nNK control has a different target; overlapping calls are not independent confirmations.\nOriginal literal-failure slot is retained in the parent run. No clinical test or formal identity replacement.',fontsize=9)
    fig.tight_layout(rect=(0,.15,1,.94))
    for ext in ['png','svg']:fig.savefig(root/('source_complete_normalized_model_reference_comparison.'+ext),dpi=170,bbox_inches='tight')
    plt.close(fig)
    report=dict(status='COMPLETE_ALL_SOURCE_REFERENCE_COUNTS_OVERLAPS_AND_FIGURE_VISUAL_REVIEW_PENDING',created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),models=4,reference_groups=5,author_model_positions=204,joint_metadata_positions=2448,individual_cells_or_samples_displayed=False,formal_identity_or_clinical_inference_changed=False,scientific_outputs_public=False,files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in root.iterdir() if p.is_file() and p.name!='private_source_report_terminal.json'})
    p=root/'private_source_report_terminal.json';assert not p.exists();p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],terminal_sha256=sha(p))))


if __name__=='__main__':
    main()
