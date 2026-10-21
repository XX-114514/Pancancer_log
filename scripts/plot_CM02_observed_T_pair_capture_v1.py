"""Private source-level complete primary forest, without significance labels."""
from datetime import datetime,timezone
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R=['R_TNFRSF9','R_LAG3','R_cycling']
E=['E_dysfunction','E_cycling','E_memory']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('source',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();root=args.source;out=args.output
    assert not out.exists();out.mkdir(mode=0o700)
    t=json.loads((root/'private_terminal.json').read_text())
    assert t['status']=='COMPLETE_ALL288_OBSERVED_T_PAIR_CAPTURE_DESCRIPTIONS'
    name='private_all288_observed_T_pair_capture_descriptions.json'
    source=root/name;assert sha(source)==t['files'][name]['SHA256']
    rows=json.loads(source.read_text())
    chosen=[r for r in rows if r['context']=='full_current_source' and r['parent']=='dual_regulatory_evidence' and r['capture_scope']=='all_native_filtered']
    assert len(chosen)==18
    lookup={(r['phase'],r['regulatory_state'],r['effector_state']):r for r in chosen}
    pairs=list(product(R,E));fig,axes=plt.subplots(1,2,figsize=(14,7),sharey=True)
    styles=[('unadjusted_Spearman',None,.22,'#727272','o','Unadjusted'),
            ('abundance_treatment_HPV_partial_Spearman','base',0,'#27648c','s','Counts + treatment + HPV'),
            ('capture_adjusted_partial_Spearman','capture',-.22,'#be6524','D','Add RNA capture proxies')]
    for ax,phase,title in zip(axes,['pre','post'],['Pre-treatment: primary description','Post-treatment: same cohort, secondary']):
        counts=set()
        for j,(r,e) in enumerate(pairs):
            row=lookup[phase,r,e];counts.add(row['reviewed_cases']);y=len(pairs)-1-j
            if j%2==0:ax.axhspan(y-.45,y+.45,color='#f4f5f6',zorder=0)
            for key,prefix,offset,color,marker,label in styles:
                point=row[key]
                if point is None:
                    ax.text(1.025,y+offset,'NE',color=color,fontsize=7,va='center')
                    continue
                interval=row[prefix+'_pointwise_CI'] if prefix else None
                if interval is not None:
                    ax.hlines(y+offset,*interval,color=color,linewidth=1.4,zorder=2)
                ax.scatter(point,y+offset,color=color,marker=marker,s=25,zorder=3)
        assert len(counts)==1
        ax.axvline(0,color='#888888',linestyle='--',linewidth=.8)
        ax.set_xlim(-1.08,1.13);ax.set_xticks([-1,-.5,0,.5,1]);ax.set_ylim(-.6,8.6)
        ax.set_title(title+'\nEligible reviewed patients: '+str(next(iter(counts))),fontsize=10)
        ax.set_xlabel('Rank correlation / partial rank correlation')
        ax.grid(axis='x',alpha=.15)
    axes[0].set_yticks(range(9))
    axes[0].set_yticklabels([r+' / '+e for r,e in pairs][::-1],fontsize=9)
    legend=[Line2D([0],[0],marker=m,color=c,linestyle='-',markersize=5,label=label) for _,_,_,c,m,label in styles]
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.58,.055),ncol=3,frameon=False,fontsize=9)
    fig.suptitle('GSE286827: all nine observed regulatory / effector T-state pairs',fontsize=13)
    fig.text(.02,.014,'Pointwise 95% patient bootstrap intervals; no multiplicity-controlled claim. NE: not estimable.\nRNA proxies are endogenous; these descriptions do not establish direct communication or independent validation.',fontsize=8,color='#444444')
    fig.subplots_adjust(left=.30,right=.98,top=.86,bottom=.20,wspace=.08)
    image=out/'all9_pre_post_adjustment_forest.png';fig.savefig(image,dpi=200)
    pdf=out/'all9_pre_post_adjustment_forest.pdf';fig.savefig(pdf)
    plt.close(fig)
    receipt=dict(status='COMPLETE_PRIVATE_ALL9_PRIMARY_AND_SECONDARY_SOURCE_FOREST',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(source),plotter_SHA256=sha(Path(__file__)),
                 no_individual_records_or_learned_weights=True,all18_primary_secondary_source_positions_shown=True,
                 files={p.name:dict(SHA256=sha(p),bytes=p.stat().st_size) for p in [image,pdf]})
    (out/'private_plot_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status=receipt['status'])))


if __name__=='__main__':main()
