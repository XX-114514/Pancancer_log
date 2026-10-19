"""All eight reference contexts, two fixed models and conditional loss intervals."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    out = root / 'private_source_figures_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    t = json.loads((root / 'private_terminal.json').read_text())
    path = root / 'private_all8_reference_contexts.json'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == t['files'][path.name]['SHA256']
    records = json.loads(path.read_text()); assert len(records) == 8
    labels = [x['heldout_source'] + ' / ' + ('4-class lineage' if x['task'] == 'coarse_lineage' else '7-class T RNA capture') for x in records]
    fig, axes = plt.subplots(1,2,figsize=(14,7.6),sharey=True)
    for i, x in enumerate(records):
        if x['status'] != 'COMPLETE_FIXED_SOURCE_HELDOUT_DIAGNOSTIC':
            for ax in axes: ax.text(.02,i,'Not estimable',va='center',color='0.5',fontsize=9)
            continue
        if x['all_heldout_classes_support_full_macro']:
            q = x['metrics']['QC_only']['macroF1']; r = x['metrics']['RNA_PCs_plus_QC']['macroF1']
            axes[0].plot([q,r],[i,i],color='0.65',lw=1)
            axes[0].plot(q,i,'o',color='#888888',label='QC only' if i==0 else None)
            axes[0].plot(r,i,'o',color='#236c94',label='RNA PCs + QC' if i==0 else None)
        else: axes[0].text(.02,i,'Full-class support insufficient',va='center',fontsize=8)
        delta = x['paired_class_material_balanced_logloss_improvement']; b = x['conditional_library_bootstrap']
        if 'CI_low' in b:
            axes[1].errorbar(delta,i,xerr=np.array([[delta-b['CI_low']],[b['CI_high']-delta]]),fmt='o',color='#236c94',capsize=3,lw=1.2)
        else: axes[1].plot(delta,i,'o',color='#236c94')
    axes[0].set_yticks(range(8),labels,fontsize=9); axes[0].set_xlim(0,1)
    axes[0].set_xlabel('Class/material-balanced macro F1')
    axes[0].set_title('Recovery of frozen operational labels',fontsize=11)
    axes[0].legend(loc='lower right',frameon=False,fontsize=9)
    axes[1].axvline(0,color='0.35',lw=.8,ls='--'); axes[1].set_xlim(-.02,1.02)
    axes[1].set_title('Paired heldout loss improvement over QC',fontsize=11)
    axes[1].set_xlabel('Balanced logloss: QC minus RNA + QC; positive = better')
    axes[0].set_ylim(7.6,-.7)
    for ax in axes:
        ax.axhline(3.5,color='0.7',lw=.7); ax.grid(axis='x',color='0.91'); ax.spines[['top','right']].set_visible(False)
    fig.subplots_adjust(left=.29,right=.98,top=.84,bottom=.22,wspace=.16)
    fig.suptitle('2026-10-14 immune atlas CM (v1): marker-disjoint reference qualification\nEvery source held out completely; all 8 contexts / two fixed models',fontsize=13)
    fig.text(.025,.045,'Labels are operational RNA/annotation classes; recovery is not functional, receptor, clinical or deconvolution validation.\nIntervals resample original material libraries with fixed heldout predictions; not retraining uncertainty or reviewed-patient inference.\nClass/material balancing describes this diagnostic sample, not tissue prevalence; grouping genes are excluded from all training features.',fontsize=9)
    for ext in ['png','svg']: fig.savefig(out/('private_all8_reference_contexts_and_loss.'+ext),dpi=200)
    plt.close(fig)
    files={p.name:dict(bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in out.iterdir()}
    assert sum(x['bytes'] for x in files.values())<4*1024**2
    result=dict(status='COMPLETE_SOURCE_DIAGNOSTIC_FIGURE_VISUAL_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parent_terminal_SHA256=hashlib.sha256((root/'private_terminal.json').read_bytes()).hexdigest(),files=files)
    (out/'private_figure_terminal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],files=files)))


if __name__=='__main__': main()
