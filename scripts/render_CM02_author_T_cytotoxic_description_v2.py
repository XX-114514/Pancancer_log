"""All four original samples; conditional RNA description and complete exclusions."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main(a):
    root=a.root;d=root/'cytotoxic_author_T_description_v1'
    audit=root/'cytotoxic_description_independent_audit_v2/independent_description_audit_technical_public.json'
    assert json.loads(audit.read_text())['producer_receipt_sha256']==hashlib.sha256((d/'description_technical_public.json').read_bytes()).hexdigest()
    table=pd.read_csv(d/'private_author_T_only_cytotoxic_RNA_description.tsv.gz',sep='\t',dtype={'NAME':str})
    out=root/'conditional_cytotoxic_figures_v2';assert not out.exists();out.mkdir(mode=0o700)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    states=['not author T','pan-T <2','no predictive RNA','reference domain outlier','prediction outside [0,1]','display eligible']
    names=[];summaries=[]
    for donor in ['Pt38','Pt43','Pt67','Pt73']:
        xy=pd.read_csv(a.source/'single_cell_portal'/('spatial_'+donor+'.txt.gz'),sep='\t',skiprows=[1],dtype={'NAME':str}).set_index('NAME')
        t=table[table.donor_id.eq(donor)].set_index('NAME');assert t.index.isin(xy.index).all()
        t=xy[['X','Y']].join(t[['cytotoxic_RNA_description','description_status_code']],how='left')
        codes=t.description_status_code.fillna(0).to_numpy(int);valid=codes==5
        assert valid.sum()==t.cytotoxic_RNA_description.notna().sum()
        fig,axes=plt.subplots(1,2,figsize=(12,5.5),gridspec_kw={'width_ratios':[3,1]})
        ax=axes[0]
        ax.scatter(t.X,t.Y,s=.25,color='#d7d9dc',rasterized=True,linewidths=0)
        im=ax.scatter(t.X[valid],t.Y[valid],c=t.cytotoxic_RNA_description[valid],s=3,cmap='viridis',vmin=0,vmax=1,linewidths=0,rasterized=True)
        ax.set_aspect('equal');ax.set(xlabel='Published X coordinate (physical unit unverified)',ylabel='Published Y coordinate',title=donor+': original-author-T conditional RNA description')
        fig.colorbar(im,ax=ax,shrink=.7,label='Predicted cytotoxic-associated RNA (new bridge)')
        n=np.bincount(codes,minlength=6);summaries.append(dict(donor_id=donor,all_original_cells=len(t),**{states[i]:int(n[i]) for i in range(6)}))
        axes[1].barh(range(6),n,color=['#b7bdc4','#e3b370','#e3b370','#e3b370','#e3b370','#187c64'])
        axes[1].set(yticks=range(6),yticklabels=states,xlabel='Original cells (log scale)',xscale='log',title='All cells accounted for')
        axes[1].set_xlim(.8,max(n.max()*4,10));axes[1].invert_yaxis()
        for i,v in enumerate(n):axes[1].text(max(v,1),i,' '+str(v),va='center',fontsize=8)
        fig.text(.03,.02,'All original positions shown. Grey = outside conditional display scope or uncertain; never a negative cytotoxic label.\nAuthor T identity is conditional, not independently validated. No colocalization, communication or clinical P; no contact claim.',fontsize=8)
        fig.subplots_adjust(wspace=.7,bottom=.2);name=donor+'_conditional_cytotoxic_RNA';names.append(name)
        for ext in ['png','pdf','svg']:fig.savefig(out/(name+'.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
        fig.savefig(out/(name+'_preview.jpg'),dpi=120,bbox_inches='tight',facecolor='white');plt.close(fig)
    pd.DataFrame(summaries).to_csv(out/'private_all4_complete_description_and_unknown_coverage.tsv',sep='\t',index=False)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Conditional single-arm description</title><main style="max-width:1200px;margin:auto;font:16px system-ui"><h1>All four samples: conditional RNA description</h1><p>Only the passed cytotoxic RNA arm; all unknowns retained. Not full CM02 spatial validation.</p>'+''.join('<section><h2>'+n+'</h2><img style="width:100%" src="'+n+'.png"><p><a href="'+n+'.pdf">PDF</a> | <a href="'+n+'.svg">SVG</a></p></section>' for n in names)+'</main>')
    proof=dict(status='ALL4_CONDITIONAL_MAPS_AND_ALL_UNKNOWN_STATES_RENDERED_VISUAL_REVIEW_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),all4_original_samples=True,
               independent_audit_sha256=hashlib.sha256(audit.read_bytes()).hexdigest(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),no_new_spatial_P=True)
    (out/'render_technical_public.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--source',type=Path,required=True);main(p.parse_args())
