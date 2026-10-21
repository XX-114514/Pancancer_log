"""All source-level state-count contrasts; pointwise intervals, no clinical claims."""
from pathlib import Path
import hashlib,json,os
from datetime import datetime,timezone
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(os.environ['PANCANCER_PARENT_OUTPUT']);out=root/'private_source_figures_v1';assert not out.exists()
    p=root/'private_all56_source_count_combinations.json';t=json.loads((root/'private_terminal.json').read_text())
    assert hashlib.sha256(p.read_bytes()).hexdigest()==t['files'][p.name]['SHA256']
    rows=json.loads(p.read_text());assert len(rows)==56
    reg=list(dict.fromkeys(r['reg_author_state'] for r in rows));other=list(dict.fromkeys(r['other_author_T_state'] for r in rows));assert len(reg)==4 and len(other)==7
    fig,axes=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(left=.13,right=.97,bottom=.21,top=.87,wspace=.42,hspace=.4)
    for ax,a in zip(axes.flat,reg):
        for phase,offset,color,name in [('pre',-.12,'#2066a8','Before treatment'),('post',.12,'#bd5800','After treatment')]:
            subset=[r for r in rows if r['reg_author_state']==a and r['phase']==phase]
            lookup={r['other_author_T_state']:r for r in subset}
            for j,b in enumerate(other):
                r=lookup[b];x=r['partial_Spearman'];lo,hi=r['pointwise_bootstrap_CI']
                if x is not None:
                    if lo is not None:ax.hlines(j+offset,lo,hi,color=color,lw=1.5)
                    ax.plot(x,j+offset,'o',color=color,ms=4,label=name if j==0 else None)
                else:ax.text(0,j+offset,'NA',fontsize=8)
        ax.set_title(a,fontsize=12);ax.set_yticks(range(7),other,fontsize=9);ax.invert_yaxis()
        ax.set_xlim(-1.03,1.03);ax.set_xticks([-1,-.5,0,.5,1]);ax.axvline(0,color='#888888',lw=.8,ls='--');ax.grid(axis='x',alpha=.15)
        ax.set_xlabel('Partial Spearman correlation',fontsize=9)
    fig.suptitle('GSE286827: all author T-state count combinations',fontsize=15,y=.96)
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.12),ncol=2,frameon=False)
    fig.text(.04,.025,'Separate denominators: state share within author Treg versus within other author T cells. Adjusted for treatment, HPV and both parent counts.\nIntervals: pointwise stratified patient bootstrap, without multiplicity-controlled pair declarations. The two phases are not independent cohorts.\nAuthor Teff is not automatically confirmed CD8; exact CM02 cell mapping, spatial structure and clinical increment are not established.',fontsize=9)
    out.mkdir(mode=0o700)
    for ext in ['png','svg']:fig.savefig(out/('private_all56_author_count_structure.'+ext),dpi=160)
    plt.close(fig)
    result=dict(status='RENDERED_ALL56_SOURCE_COUNT_COMBINATIONS_VISUAL_REVIEW_PENDING',created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_SHA256=hashlib.sha256(p.read_bytes()).hexdigest(),files={p.name:dict(SHA256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in out.iterdir()})
    (out/'private_figure_terminal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'])))


if __name__=='__main__':main()
