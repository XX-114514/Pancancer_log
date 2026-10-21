"""Render all nine source-level primary contrasts; no individual points."""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(a):
    assert sha(a.input)==a.input_sha256
    all_rows=json.loads(a.input.read_text())
    def selected(null):
        return [r for r in all_rows if r['scope']=='all_HNSCC'
            and r['radius_um']==50 and r['minimum_human_gene_molecules']==20
            and r['regulatory_parent']=='dual_regulatory_evidence' and r['null']==null]
    main_rows=selected('local_matched');global_rows=selected('parent_only')
    assert len(main_rows)==len(global_rows)==9
    lookup={(r['R_observed_state'],r['E_observed_state']):r for r in global_rows}
    labels=[];y=np.arange(9)
    names={'TNFRSF9_positive':'TNFRSF9+','LAG3_positive':'LAG3+',
           'cycling_marker_positive':'cycling markers',
           'dysfunction_marker_positive':'dysfunction markers',
           'memory_marker_positive':'TCF7/IL7R markers'}
    fig,ax=plt.subplots(figsize=(13.0,7.8))
    for k,row in enumerate(main_rows):
        key=(row['R_observed_state'],row['E_observed_state'])
        labels.append('R: '+names[key[0]]+'  /  CD8: '+names[key[1]])
        for record,offset,color,marker in [(lookup[key],-.12,'#969DA6','o'),(row,.12,'#2066A8','s')]:
            if 'conditional_excess' not in record:continue
            point=100*record['conditional_excess'];lo,hi=np.array(record['CI95_case_bootstrap'])*100
            ax.plot([lo,hi],[k+offset,k+offset],color=color,lw=1.7)
            ax.plot(point,k+offset,marker,color=color,ms=5)
        q=row.get('BH_nine_primary_opportunities')
        annotation='n='+str(row['supporting_cases'])+'; q='+('NA' if q is None else f'{q:.3f}')
        ax.text(1.02,k,annotation,transform=ax.get_yaxis_transform(),va='center',fontsize=9)
    ax.axvline(0,color='#333333',ls='--',lw=1)
    ax.set_yticks(y,labels,fontsize=9.4);ax.invert_yaxis()
    ax.set_xlabel('Standardized excess cross-parent state pairs (%)',fontsize=11)
    ax.grid(axis='x',alpha=.18)
    ax.spines[['top','right','left']].set_visible(False)
    ax.plot([],[],'o-',color='#969DA6',label='Within-parent shuffling')
    ax.plot([],[],'s-',color='#2066A8',label='Local matched shuffling (primary)')
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.03),frameon=False,ncol=2,fontsize=10)
    fig.suptitle('GSE300147 HNSCC: state proximity at the fixed 50 µm radius',fontsize=14,y=.97)
    fig.text(.5,.917,'All nine prespecified combinations; source cases equally weighted',ha='center',fontsize=10)
    fig.text(.035,.10,'R: T RNA supported by FOXP3 + IL2RA. Patterns are not confirmed functional subtypes.',fontsize=9)
    fig.text(.035,.075,'Pointwise 95% intervals: source-case bootstrap. q: nine primary conditional randomization tests.',fontsize=9)
    fig.text(.035,.05,'Excess = (observed − expected)/(expected + 0.5). ICI timing is not established for this source.',fontsize=9)
    fig.text(.035,.025,a.review_label,fontsize=8,color='#555555')
    fig.subplots_adjust(left=.355,right=.82,top=.83,bottom=.19)
    assert not a.output.exists();a.output.mkdir(mode=0o700)
    for suffix in ['png','pdf']:fig.savefig(a.output/('all9_primary_source_forest.'+suffix),dpi=170,facecolor='white')
    plt.close(fig)
    proof=dict(source_SHA256=sha(__file__),input_SHA256=sha(a.input),
        source_only_no_cell_case_points_or_weights=True,all_nine_primary_slots_rendered=True,
        files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in a.output.iterdir()})
    (a.output/'private_render_receipt.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(proof))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--input-sha256',required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--review-label',default='Independent full primary random-stream replay pending; no clinical validation.')
    main(p.parse_args())
