"""Private layout-only v2; preserve original figures and every scientific value."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);src=r/'NSCLC_paired_clinical_complete_evidence_review_20261016_v1';proof=src/'private_independent_display_audit.json'
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL572')
    source=src/'private_all572_complete_same_context_evidence.tsv';terminal=json.loads((src/'private_terminal.json').read_text());assert sha(source)==terminal['files'][source.name]
    t=pd.read_csv(source,sep='\t');out=r/'NSCLC_paired_clinical_forest_layout_20261016_v2';assert not out.exists();out.mkdir(mode=0o700)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8});records=[]
    for (context,program),f in t[t.pair_IUT_P.notna()].groupby(['context','program']):
        f=f.sort_values(['author_subtype_A','author_subtype_B']);height=max(5.8,len(f)*.24+2.2);y=np.arange(len(f))
        fig,axes=plt.subplots(1,2,figsize=(15,height),sharey=True)
        for j,key in enumerate(['A','B']):
            ax=axes[j]
            for i,row in enumerate(f.itertuples()):
                color='#206b93' if row.both_conditional_clinical_supported else '#888888';lo=getattr(row,key+'_CI_low');hi=getattr(row,key+'_CI_high');point=getattr(row,key+'_beta')
                if np.isfinite(lo) and np.isfinite(hi):ax.plot([lo,hi],[i,i],lw=.8,color=color)
                ax.scatter(point,i,s=14,color=color)
            ax.axvline(0,color='#555555',lw=.7,ls='--');ax.set(title=key+' programme conditional on the other end',xlabel='Conditional MPR association per matched-case programme SD')
        axes[0].set_yticks(y,[f'{x.record} | {x.author_subtype_A} / {x.author_subtype_B}' for x in f.itertuples()],fontsize=7);axes[0].invert_yaxis()
        title=fig.suptitle(f'Context {context}, {program}: all estimable pairs; one post-treatment NSCLC source',fontsize=11,y=1-.12/height)
        fig.subplots_adjust(left=.34,right=.98,bottom=.68/height,top=1-.80/height,wspace=.15)
        fig.canvas.draw();renderer=fig.canvas.get_renderer();mainbox=title.get_window_extent(renderer)
        for ax in axes:
            assert not mainbox.overlaps(ax.title.get_window_extent(renderer)),'Title overlap remains'
            for text in [ax.title,ax.xaxis.label]+ax.get_xticklabels()+ax.get_yticklabels():
                if not text.get_visible():continue
                box=text.get_window_extent(renderer)
                assert box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.x1+1 and box.y1<=fig.bbox.y1+1,'Text clipping remains'
        p=out/f'private_complete_forest_C{context}_{program}_v2.png';fig.savefig(p,dpi=130);plt.close(fig)
        with Image.open(p) as im:im.verify()
        data=out/f'private_plot_rows_C{context}_{program}.tsv';f.to_csv(data,sep='\t',index=False,na_rep='NA')
        replay=pd.read_csv(data,sep='\t');pd.testing.assert_frame_equal(f.reset_index(drop=True),replay,check_dtype=False)
        records.append(dict(context=int(context),program=program,rows=len(f),png=p.name,png_sha256=sha(p),plot_data_sha256=sha(data),title_overlap_and_all_text_bounds_checked=True))
    assert len(records)==8 and sum(x['rows'] for x in records)==360
    result=dict(status='PASS_ALL8_LAYOUT_V2_NUMERIC_ROWS_UNCHANGED_ALL_TITLES_AND_TEXT_BOUNDS_CHECKED_VISUAL_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
                code_sha256=sha(__file__),input_source_sha256=sha(source),input_audit_sha256=sha(proof),records=records,original_v1_figures_preserved=True,new_inference=False)
    (out/'private_layout_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
