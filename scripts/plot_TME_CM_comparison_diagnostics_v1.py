"""Descriptive plots from accepted internal CM comparison summaries only."""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)

def main(R):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    cfg=load(R/'config/spec_v1.json');P=Path(cfg['comparison'])
    assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
    p=P/'outputs/accepted_membership_module_comparison_v1.json'
    assert sha(p)==load(P/'provenance/terminal_v1.json')['acceptance_SHA256'];m=load(p)
    schemes=[x['scheme'] for x in m['crosswalk_checks']];n=len(schemes)
    ranks={}
    for r in m['unmatched']:ranks[r['left_scheme']]=r['left_rank'];ranks[r['right_scheme']]=r['right_rank']
    names=cfg['display_names'];labels=[names[s]+'\nK='+str(ranks[s]) for s in schemes]
    pairs={}
    for r in m['matches']:
        key=(r['basis'],r['left_scheme'],r['right_scheme']);pairs.setdefault(key,[]).append(r)
    figures=[];fig,axes=plt.subplots(1,3,figsize=(15,5.7),layout='constrained')
    for ax,(basis,field,title,vmin,vmax,cmap) in zip(axes,[
        ('W','cosine','Raw W: matched loading similarity',0,1,'Blues'),
        ('A','cosine','Effective W S: matched loading similarity',0,1,'Blues'),
        ('W','activity_Spearman','Same-unit activity: W-matched factors',-1,1,'RdBu_r')]):
        matrix=np.full((n,n),np.nan)
        for i,left in enumerate(schemes):
            for j,right in enumerate(schemes):
                rows=pairs.get((basis,left,right),pairs.get((basis,right,left),[]))
                values=[r[field] for r in rows if r[field] is not None]
                if values:matrix[i,j]=np.median(values)
        cm=plt.get_cmap(cmap).copy();cm.set_bad('#e5e7eb')
        im=ax.imshow(matrix,cmap=cm,norm=Normalize(vmin,vmax))
        ax.set_xticks(range(n),labels,rotation=38,ha='right',fontsize=8)
        ax.set_yticks(range(n),labels,fontsize=8);ax.set_title(title,fontsize=10,pad=13)
        for i in range(n):
            for j in range(n):
                value=matrix[i,j];text='unmatched' if not np.isfinite(value) else f'{value:.2f}'
                dark=np.isfinite(value) and (value>.75 if vmin==0 else abs(value)>.65)
                ax.text(j,i,text,ha='center',va='center',fontsize=8 if np.isfinite(value) else 6,color='white' if dark else '#111827')
        fig.colorbar(im,ax=ax,shrink=.65,pad=.02)
    fig.suptitle('First T/NK–myeloid panel: annotation-workflow comparison',fontsize=14)
    fig.supxlabel('Medians of recorded matches; two final seeds. Diagonal = seed repeat. Gray = no W match above cosine 0.60.\nOwn ranks differ; fixed-K and source checks are separate. Normalized feature loadings are not cell counts.',fontsize=9)
    for ext in ['png','pdf']:
        out=R/'outputs'/('first_panel_CM_similarity_and_activity_v1.'+ext);fig.savefig(out,dpi=180,bbox_inches='tight');figures.append(dict(path=str(out),SHA256=sha(out)))
    plt.close(fig)
    # Seed stability is an internal numerical reproducibility diagnostic.
    fig,ax=plt.subplots(figsize=(8.5,4.6),layout='constrained');x=np.arange(n)
    for field,label,color,offset in [('cosine','Raw W loading similarity','#2563eb',-.15),('activity_Spearman','Same-unit activity Spearman','#059669',.15)]:
        values=[np.median([r[field] for r in pairs['W',s,s] if r[field] is not None]) for s in schemes]
        ax.bar(x+offset,values,.28,label=label,color=color)
    ax.set_xticks(x,labels,fontsize=9);ax.set_ylim(0,1.08);ax.set_ylabel('Median matched factor agreement')
    ax.legend(loc='lower left',fontsize=9);ax.set_title('Within-workflow repeatability across two fixed seeds')
    ax.spines[['top','right']].set_visible(False);ax.text(0,-.26,'Training-domain description; not external or independent-patient validation.',transform=ax.transAxes,fontsize=9)
    for ext in ['png','pdf']:
        out=R/'outputs'/('first_panel_CM_seed_repeatability_v1.'+ext);fig.savefig(out,dpi=180,bbox_inches='tight');figures.append(dict(path=str(out),SHA256=sha(out)))
    plt.close(fig)
    save(R/'outputs/figure_manifest_v1.json',dict(input_comparison_SHA256=sha(p),figures=figures,
        matplotlib_version=matplotlib.__version__,rank_and_matching_scope_retained=True,private_analysis_figures_not_publication_release=True))
    save(R/'provenance/terminal_v1.json',dict(status='ACCEPTED_CM_COMPARISON_DESCRIPTIVE_FIGURES_RENDERED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),manifest_SHA256=sha(R/'outputs/figure_manifest_v1.json')))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
    try:main(R)
    except BaseException:
        import traceback
        save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
