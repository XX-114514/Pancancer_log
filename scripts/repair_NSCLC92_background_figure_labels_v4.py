"""Independent figure-coordinate identifier labels; same audited numerical data."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);out=root/'NSCLC_frozen92_axis_source_background_labels_20261019_v4';assert not out.exists();out.mkdir(mode=0o700);proof=root/'NSCLC_frozen92_axis_source_result_audit_20261019_v1/private_terminal.json';assert sha(proof)=='79e8aa29aec12db05d437f6e4ebc24dec7e97772bf62ca594434efe36522b9b7';src=root/'NSCLC_frozen92_axis_source_consistency_20261019_v1';t=json.loads((src/'private_terminal.json').read_text());names=['private_all92_source_consistency_results.tsv','private_all276_fixed_sensitivity_effects.tsv']
    for n in names:assert sha(src/n)==t['files'][n]['sha256']
    r=pd.read_csv(src/names[0],sep='\t');b=pd.read_csv(src/names[1],sep='\t');b=b[b.sensitivity.eq('main_plus_total_communication')].set_index('axis_id');assert len(r)==92
    plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False});fig,axes=plt.subplots(2,2,figsize=(15,13));fig.subplots_adjust(left=.13,right=.93,top=.91,bottom=.09,hspace=.25,wspace=.35);fig.suptitle('GSE243013: all 92 primary and total-background effects\nBlue: WCU-B pointwise intervals. Orange: additional-background CV3 t intervals, without new FDR testing.',fontsize=12,y=.99);panels=[r.iloc[i:i+23] for i in range(0,92,23)];lo=min(r.CI_low.min(),b.CV3_t_CI_low.min())-.12;hi=max(r.CI_high.max(),b.CV3_t_CI_high.max())+.12
    for ax,part in zip(axes.ravel(),panels):
        ax.set_yticks(np.arange(len(part)));ax.set_yticklabels(['']*len(part));ax.set_ylim(len(part)-.4,-.8);ax.set_xlim(lo,hi);ax.axvline(0,color='black',lw=.7);ax.grid(axis='x',alpha=.15);ax.set_xlabel('Adjusted pCR/MPR minus non-MPR difference (axis SD)')
        for y,row in enumerate(part.itertuples()):
            rr=b.loc[row.axis_id]
            for offset,value,lower,upper,colour,marker in [(-.12,row.effect,row.CI_low,row.CI_high,'#266ab2','o'),(.12,rr.effect,rr.CV3_t_CI_low,rr.CV3_t_CI_high,'#c67b36','^')]:
                if np.isfinite(value):ax.errorbar(value,y+offset,xerr=[[value-lower],[upper-value]],fmt=marker,color=colour,ms=3,lw=.7,capsize=1)
    fig.canvas.draw();labels=[]
    for ax,part in zip(axes.ravel(),panels):
        for y,identifier in enumerate(part.axis_id.astype(str)):
            xy=fig.transFigure.inverted().transform(ax.transData.transform((0.,y)));text=fig.text(ax.get_position().x0-.014,xy[1],identifier,ha='right',va='center',fontsize=9,color='black',gid='frozen_axis_id');labels.append(text.get_text())
    assert labels==r.axis_id.tolist() and len(set(labels))==92
    fig.legend(handles=[Line2D([0],[0],color='#266ab2',marker='o',label='Primary; WCU-B interval'),Line2D([0],[0],color='#c67b36',marker='^',label='Additional background; CV3 t interval')],loc='lower center',ncol=2,frameon=False)
    fig.savefig(out/'all92_main_and_total_background_explicit_figure_labels.png',dpi=150);fig.savefig(out/'all92_main_and_total_background_explicit_figure_labels.pdf');plt.close(fig);files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(d['bytes'] for d in files.values())<2*1024**2
    receipt=dict(status='COMPLETE_ALL92_BACKGROUND_FIGURE_LABELS_SEMANTICALLY_VERIFIED_ACTUAL_VISUAL_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,input_result_audit_sha256=sha(proof),complete92_frozen_identifier_text_objects_verified=True,files=files,scientific_data_changed=False,scientific_outputs_public=False);(out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
