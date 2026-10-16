"""V3: check rendered in-view ticks and reserve fixed physical title/caption margins."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rendered_text(ax):
    # Axis draws ticks only in its current interval. get_xticklabels() also returns
    # visible-flag labels outside that interval, which are not actually rendered.
    text=[ax.title,ax.xaxis.label,ax.yaxis.label,ax.xaxis.get_offset_text(),ax.yaxis.get_offset_text()]
    for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
        lo,hi=sorted(limits)
        for tick in axis.get_major_ticks():
            if lo-1e-10<=tick.get_loc()<=hi+1e-10:text.extend([tick.label1,tick.label2])
    return [t for t in text if t.get_visible() and t.get_text()]


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);src=r/'NSCLC_paired_clinical_complete_evidence_review_20261016_v1';proof=src/'private_independent_display_audit.json'
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL572')
    source=src/'private_all572_complete_same_context_evidence.tsv';original=json.loads((src/'private_terminal.json').read_text());assert sha(source)==original['files'][source.name]
    t=pd.read_csv(source,sep='\t');out=r/'NSCLC_paired_clinical_forest_layout_20261016_v3';assert not out.exists();out.mkdir(mode=0o700)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8});records=[];passed=int(t.BH_q_pair572.lt(.05).sum())
    for (context,program),f in t[t.pair_IUT_P.notna()].groupby(['context','program']):
        f=f.sort_values(['author_subtype_A','author_subtype_B']);height=max(6.2,len(f)*.24+2.2);y=np.arange(len(f))
        fig,axes=plt.subplots(1,2,figsize=(15,height),sharey=True)
        for key,ax in zip(['A','B'],axes):
            for i,row in enumerate(f.itertuples()):
                color='#206b93' if row.both_conditional_clinical_supported else '#888888';lo=getattr(row,key+'_CI_low');hi=getattr(row,key+'_CI_high')
                if np.isfinite(lo) and np.isfinite(hi):ax.plot([lo,hi],[i,i],lw=.8,color=color)
                ax.scatter(getattr(row,key+'_beta'),i,s=14,color=color)
            ax.axvline(0,color='#555555',lw=.7,ls='--');ax.set(title=key+' programme conditional on the other end',xlabel='Conditional MPR association per matched-case programme SD')
        axes[0].set_yticks(y,[f'{x.record} | {x.author_subtype_A} / {x.author_subtype_B}' for x in f.itertuples()],fontsize=7);axes[0].invert_yaxis()
        title=fig.suptitle(f'Context {context}, {program}: all estimable pairs; one post-treatment NSCLC source',fontsize=11,y=1-.12/height)
        caption=fig.text(.34,.11/height,f'Pointwise 95% bootstrap intervals; {passed}/{len(t)} planned pairs passed the primary joint clinical BH threshold (q<0.05).',fontsize=8)
        fig.subplots_adjust(left=.34,right=.98,bottom=.90/height,top=1-.80/height,wspace=.15)
        fig.canvas.draw();renderer=fig.canvas.get_renderer();mainbox=title.get_window_extent(renderer);checked=0;offview_outside_figure=0
        for ax in axes:
            assert not mainbox.overlaps(ax.title.get_window_extent(renderer)),'Title overlap remains'
            actual=rendered_text(ax)
            for text in actual+[title,caption]:
                box=text.get_window_extent(renderer);inside=box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.x1+1 and box.y1<=fig.bbox.y1+1
                assert inside,(int(context),program,text.get_text(),list(box.bounds),list(fig.bbox.bounds));checked+=1
            actual_ids={id(v) for v in actual}
            for text in ax.get_xticklabels()+ax.get_yticklabels():
                box=text.get_window_extent(renderer)
                if id(text) not in actual_ids and text.get_visible() and (box.x0<0 or box.y0<0 or box.x1>fig.bbox.x1 or box.y1>fig.bbox.y1):offview_outside_figure+=1
        p=out/f'private_complete_forest_C{context}_{program}_v3.png';fig.savefig(p,dpi=130);plt.close(fig)
        with Image.open(p) as im:im.verify()
        data=out/f'private_plot_rows_C{context}_{program}.tsv';f.to_csv(data,sep='\t',index=False,na_rep='NA');pd.testing.assert_frame_equal(f.reset_index(drop=True),pd.read_csv(data,sep='\t'),check_dtype=False)
        records.append(dict(context=int(context),program=program,rows=len(f),png=p.name,png_sha256=sha(p),plot_data_sha256=sha(data),rendered_text_bounds_checked=checked,offview_not_drawn_labels_outside_figure=offview_outside_figure))
    assert len(records)==8 and sum(x['rows'] for x in records)==360
    result=dict(status='PASS_ALL8_LAYOUT_V3_NUMERIC_ROWS_UNCHANGED_RENDERED_TITLES_TICKS_AND_CAPTIONS_CHECKED_VISUAL_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
                input_source_sha256=sha(source),input_audit_sha256=sha(proof),records=records,original_v1_and_failed_v2_preserved=True,new_inference=False)
    (out/'private_layout_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
