"""Scientific figure from complete source/author-role depth aggregates only."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib,json,os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


def main():
    root=Path(os.environ['PANCANCER_RAM_OUTPUT']);out=root/'private_source_figures_v2'
    assert not out.exists()
    p=root/'private_all448_author_group_depth_descriptions.tsv'
    t=json.loads((root/'private_terminal.json').read_text())
    assert hashlib.sha256(p.read_bytes()).hexdigest()==t['files'][p.name]['sha256']
    d=pd.read_csv(p,sep='\t',float_precision='round_trip')
    d=d[d.scope.eq('five_disjoint_author_roles')]
    roles=['Author_Vdelta1','Author_Vdelta2','Author_NK','Other_T_NK','Non_T_NK']
    names=['Author Vdelta1','Author Vdelta2','Author NK','Other author T/NK','Author non-T/NK']
    colors=['#005b96','#b75000','#7a5195','#2a8f68','#777777']
    fig,axes=plt.subplots(1,2,figsize=(11,5.6),sharey=True)
    fig.subplots_adjust(left=.08,right=.98,top=.81,bottom=.40,wspace=.15)
    for ax,axis,title in zip(axes,['native_total_RNA','native_detected_genes'],['Total native RNA','Detected genes']):
        for role,name,color in zip(roles,names,colors):
            s=d[d.author_group.eq(role)&d.axis.eq(axis)].sort_values('depth_bin')
            assert len(s)==4 and s.depth_bin.tolist()==[1,2,3,4]
            ax.plot(s.depth_bin,s.fixed_RNA_candidate_union_fraction,marker='o',color=color,label=name,lw=1.8,ms=5)
        ax.set_title(title,fontsize=11);ax.set_xticks([1,2,3,4],['1: lowest','2','3','4: highest'])
        ax.set_xlabel('Detection-depth interval');ax.grid(axis='y',alpha=.25)
        ax.yaxis.set_major_formatter(PercentFormatter(1));ax.set_ylim(0,float(d.fixed_RNA_candidate_union_fraction.max())+.04)
    axes[0].set_ylabel('Fixed RNA-rule capture fraction')
    fig.suptitle('GSE243013: RNA candidate capture and detection depth',fontsize=14,y=.95)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.20),ncol=3,frameon=False,fontsize=9)
    fig.text(.08,.045,'All five disjoint author roles are shown. Intervals use quartiles in the original author T/NK pool; ties stay together.\nAuthor labels are reference labels, not productive receptor truth. Rule capture is not validated sensitivity.\nSource-level descriptions; no independent-patient or causal dropout inference.',fontsize=9)
    out.mkdir(mode=0o700)
    for suffix in ['png','svg']:fig.savefig(out/('private_gamma_depth_capture_all_author_roles.'+suffix),dpi=180)
    plt.close(fig)
    files={p.name:dict(SHA256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in out.iterdir()}
    x=dict(status='RENDERED_SOURCE_DEPTH_FIGURE_VISUAL_REVIEW_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_SHA256=t['files']['private_all448_author_group_depth_descriptions.tsv']['sha256'],files=files)
    (out/'private_figure_terminal.json').write_text(json.dumps(x,indent=2)+'\n')
    print(json.dumps(dict(status=x['status'])))


if __name__=='__main__':main()
