"""Redesign an atlas opening figure from frozen, externally supplied data.

No embedding refit, inferred trajectories, model fitting or outcome testing.
Rasterized scatter layers; editable vector labels, diagrams and dotplot.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import FancyBboxPatch, Wedge
import numpy as np
import pandas as pd
from umap_labels import add_umap_labels

PARENTS=['B','Plasma','CD4T','CD8T','ILC','Myeloid','Endothelial','Stromal']
DISPLAY=['B cells','Plasma','CD4 T','CD8 T','NK / ILC','Myeloid','Endothelial','Stromal']
PALETTE=['#0072B2','#56B4E9','#009E73','#117733','#CC79A7','#E69F00','#D55E00','#7C6DAD']
INK='#17232D';MUTED='#52606B';BLUE='#365E78';LIGHT='#EDF3F6'
REGIONS={
    'Brain':['Diffuse glioma'],
    'Head / neck':['Head and neck squamous cell carcinoma','Nasopharyngeal carcinoma'],
    'Lung':['NSCLC'],
    'Breast':['Breast cancer'],
    'Hepatobiliary':['Hepatocellular carcinoma','Intrahepatic cholangiocarcinoma'],
    'Upper GI':['Gastric cancer','Esophageal cancer'],
    'Pancreas':['Pancreatic ductal adenocarcinoma','Pancreatic neuroendocrine tumor'],
    'Colorectal':['Colorectal cancer'],
    'Kidney':['Renal cell carcinoma','Renal medullary carcinoma'],
    'Bladder / prostate':['Bladder cancer','Prostate cancer'],
    'Gynecologic':['Ovarian cancer','Cervical cancer','Endometrial cancer'],
    'Skin / eye':['Cutaneous melanoma','Acral melanoma','Basal cell carcinoma','Uveal melanoma'],
    'Thyroid':['Papillary thyroid carcinoma'],
    'Mesenchymal':['Osteosarcoma','Soft tissue sarcoma'],
    'Other':['Diffuse large B-cell lymphoma','Mesothelioma'],
}
GENES=['CD79A','MS4A1','MZB1','JCHAIN','CD3D','TRAC','IL7R','FOXP3','GZMK','GZMH','XCL1','FCGR3A','LYZ','LST1','VWF','CA4','ESM1','COL1A1','RGS5','MYH11']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,
    'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':7,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'figure.facecolor':'white','savefig.facecolor':'white'})

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def panel(fig,letter,title,x,y):
    fig.text(x,y,letter,fontsize=11,fontweight='bold',color=INK,va='top')
    fig.text(x+.028,y-.001,title,fontsize=8.2,fontweight='bold',color=INK,va='top')

def off(ax):
    ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values():s.set_visible(False)

def coverage(ax,cancer_counts):
    rows=[]
    for region,names in REGIONS.items():
        for name in names:
            rows.append({'cancer':name,'display_group':region,'cells':cancer_counts[name]})
    assert set(cancer_counts)=={x['cancer'] for x in rows}
    cmap=plt.get_cmap('tab20c')
    cols=[cmap(i) for i in [0,1,4,5,8,9,12,13,16,17,2,6,10,14,18]]
    total=sum(cancer_counts.values());theta=90;anchors=[]
    for (region,names),col in zip(REGIONS.items(),cols):
        value=sum(cancer_counts[n] for n in names);delta=360*value/total
        ax.add_patch(Wedge((0,.12),.91,theta-delta,theta,width=.27,facecolor=col,edgecolor='white',lw=.7))
        # Fine radial divisions distinguish all cancer labels within each display group.
        sub=theta
        for n in names:
            sub-=360*cancer_counts[n]/total
            rad=np.deg2rad(sub)
            ax.plot([.64*np.cos(rad),.91*np.cos(rad)],[.12+.64*np.sin(rad),.12+.91*np.sin(rad)],color='white',lw=.35)
        mid=np.deg2rad(theta-delta/2)
        anchors.append({'label':region,'x':np.cos(mid),'y':.12+np.sin(mid)})
        theta-=delta
    # Two-column callouts keep even the rare categories readable.
    for side in [-1,1]:
        items=sorted([a for a in anchors if (1 if a['x']>=0 else -1)==side],key=lambda a:a['y'])
        ys=np.array([a['y'] for a in items]);gap=.20
        for k in range(1,len(ys)):ys[k]=max(ys[k],ys[k-1]+gap)
        if len(ys) and ys[-1]>1.22:ys-=ys[-1]-1.22
        for a,y in zip(items,ys):
            ax.plot([a['x']*.94,side*1.02,side*1.11],[.12+(a['y']-.12)*.94,y,y],color='#A3AEB7',lw=.45)
            label=a['label'].replace('Bladder / prostate','Bladder /\nprostate')
            ax.text(side*1.14,y,label,fontsize=5.9,ha='left' if side>0 else 'right',va='center',color=INK,linespacing=1.05)
    ax.text(0,.18,'4.72M',ha='center',va='center',fontsize=14,fontweight='bold',color=INK)
    ax.text(0,-.10,'cells',ha='center',va='center',fontsize=7.6,color=MUTED)
    ax.set_xlim(-2.1,2.1);ax.set_ylim(-1.85,1.45);ax.set_aspect('equal');off(ax)
    ax.text(0,-1.31,'72 sources · 28 cancer types',fontsize=6.4,ha='center',color=INK)
    ax.text(0,-1.60,'1,258 units · 1,014 registry patients',fontsize=6.1,ha='center',color=MUTED)
    return pd.DataFrame(rows)

def resources(ax,cfg):
    ax.set_xlim(-.02,1.02);ax.set_ylim(-.03,1.04);off(ax)
    def card(x,y,title,number,detail,accent=False):
        ax.add_patch(FancyBboxPatch((x,y),.28,.33,boxstyle='round,pad=.008,rounding_size=.015',
            fc='#E5EEF4' if accent else '#F5F7F8',ec='#93A7B7' if accent else '#C8D1D7',lw=.65))
        ax.text(x+.14,y+.255,title,ha='center',va='center',fontsize=5.7,fontweight='bold',color=INK,linespacing=1.12)
        ax.text(x+.14,y+.135,f'{number:,}',ha='center',va='center',fontsize=9.1,fontweight='bold',color=BLUE if accent else INK)
        ax.text(x+.14,y+.048,detail,ha='center',va='center',fontsize=5.0,color=MUTED)
    card(.0,.61,'Harmonised\npan-cancer atlas',cfg['harmonised']['cells'],'72 sources | 28 cancer types',True)
    card(.36,.61,'Non-malignant\nTME view',cfg['non_malignant']['cells'],'Analysis view',True)
    card(.72,.61,'Resolved TME\nidentity reference',cfg['reference']['cells'],'44 identities | 70 sources',True)
    for a,b in [(.28,.36),(.64,.72)]:
        ax.annotate('',xy=(b-.012,.775),xytext=(a+.012,.775),arrowprops=dict(arrowstyle='->',lw=.85,color=BLUE))
    card(.0,.10,'Discovery\ncollection',cfg['discovery']['cells'],'54 series | 1,322 samples')
    card(.36,.10,'Immunotherapy\ncompendium',cfg['immune']['cells'],'9 collections')
    card(.72,.10,'Curated\nT-cell states',cfg['curated']['cells'],'6 sources')
    ax.annotate('',xy=(.708,.265),xytext=(.652,.265),arrowprops=dict(arrowstyle='->',lw=.8,color=MUTED))
    ax.annotate('',xy=(.14,.592),xytext=(.14,.452),arrowprops=dict(arrowstyle='->',lw=.75,color=MUTED,linestyle='--'))
    ax.text(.167,.517,'re-ingestion',fontsize=5.1,va='center',color=MUTED)
    ax.text(.5,-.015,'Linked, frozen resources; cell totals are not additive.',ha='center',fontsize=5.8,color=MUTED)

def scatter(ax,xy,code,palette,shuffle):
    ax.scatter(xy[shuffle,0],xy[shuffle,1],c=palette[code[shuffle]],s=.07,lw=0,alpha=.72,rasterized=True)
    ax.set_aspect('equal');off(ax)
    # Small orientation axes retain most of the panel for data.
    for end in [(.18,.04),(.04,.18)]:
        ax.annotate('',xy=end,xytext=(.04,.04),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.65,color=MUTED))
    ax.text(.19,.015,'UMAP 1',transform=ax.transAxes,fontsize=5.5,ha='left',color=MUTED)
    ax.text(.005,.20,'UMAP 2',transform=ax.transAxes,fontsize=5.5,ha='center',rotation=90,color=MUTED)

def main(a):
    out=a.output;out.mkdir(parents=True,exist_ok=False);tables=out/'source_tables';tables.mkdir()
    cfg=json.loads(a.resources_json.read_text())
    registry=pd.read_csv(a.registry,sep='\t',keep_default_na=False)
    cancer=Counter()
    for text in registry.analysis_main:
        for item in text.split(';'):
            name,n=item.strip().rsplit(':',1);cancer[name]+=int(n)
    assert len(registry)==cfg['harmonised']['sources']==72
    assert sum(cancer.values())==registry.cells_harmonised_atlas.astype(int).sum()==cfg['harmonised']['cells']
    assert len(cancer)==28 and registry.units.astype(int).sum()==1258 and registry.patients.astype(int).sum()==1014
    complete=json.loads((a.atlas_dir/'complete.json').read_text())
    assert complete['R4_44_identity_cells']==cfg['reference']['cells']==1974309
    mapping=pd.read_csv(a.node_map,sep='\t').set_index('node').majorCluster.to_dict()
    order=[n for p in PARENTS for n,c in mapping.items() if c==p]
    embed=a.atlas_dir/'atlas_embedding_PRIVATE.npz'
    with np.load(embed) as z:
        xy=z['xy'];node=z['node'];nodes=z['nodes'].tolist();src=z['source'];canc=z['cancer'];sources=z['sources'].tolist();cancers=z['cancers'].tolist()
    assert nodes==order and len(nodes)==44 and len(xy)==1974309 and len(sources)==70 and len(cancers)==27
    assert np.isfinite(xy).all() and node.min()==0 and node.max()==43
    pcode=np.array([PARENTS.index(mapping[n]) for n in nodes])[node]
    parent_colors=np.array([colors.to_rgba(c) for c in PALETTE])
    fine=[]
    for p,col in zip(PARENTS,PALETTE):
        ns=[n for n in nodes if mapping[n]==p];base=np.array(colors.to_rgb(col))
        for k in range(len(ns)):
            fine.append(base*(.6+.4*k/max(1,len(ns)-1))+np.ones(3)*(.4-.4*k/max(1,len(ns)-1)))
    fine=np.clip(np.array(fine),0,1)
    shuffle=np.random.default_rng(20261120).permutation(len(xy))
    key=pd.DataFrame({'number':np.arange(1,45),'identity':nodes,'major':[mapping[n] for n in nodes],
        'color_hex':[colors.to_hex(c) for c in fine],'cells':np.bincount(node,minlength=44)})
    key.to_csv(tables/'Figure1d_identity_key.tsv',sep='\t',index=False)
    pc=pd.DataFrame({'major':PARENTS,'display':DISPLAY,'cells':np.bincount(pcode,minlength=8),'color_hex':PALETTE})
    pc.to_csv(tables/'Figure1c_major_counts.tsv',sep='\t',index=False)
    per_cancer=np.bincount(canc,minlength=len(cancers))
    actual=pd.Series(per_cancer,index=cancers).sort_index()
    expected=pd.read_csv(a.atlas_dir/'source_cancer_counts.tsv',sep='\t').groupby('cancer').cells.sum().sort_index()
    assert actual.index.equals(expected.index) and np.array_equal(actual.values,expected.values)

    stats=pd.read_csv(a.atlas_dir/'annotation_marker_statistics.tsv',sep='\t');stats['major']=stats.identity.map(mapping)
    s=stats[stats.gene.isin(GENES)].copy()
    s['weighted_mean']=s.mean_log1p_1e4*s.n_eligible_cells
    g=s.groupby(['major','gene']).agg(eligible=('n_eligible_cells','sum'),positive=('n_positive','sum'),weighted=('weighted_mean','sum'))
    data=[]
    for p in PARENTS:
        for gene in GENES:
            r=g.loc[(p,gene)];assert r.eligible>0
            data.append({'major':p,'gene':gene,'n_eligible_cells':int(r.eligible),'n_positive':int(r.positive),
                'fraction_positive':r.positive/r.eligible,'mean_log1p_1e4':r.weighted/r.eligible})
    dots=pd.DataFrame(data)
    maxima=dots.groupby('gene').mean_log1p_1e4.transform('max')
    dots['gene_scaled_mean']=dots.mean_log1p_1e4/maxima
    assert dots.fraction_positive.between(0,1).all() and dots.gene_scaled_mean.between(0,1).all()
    dots.to_csv(tables/'Figure1e_parent_marker_statistics.tsv',sep='\t',index=False)
    (tables/'resource_configuration.json').write_text(json.dumps(cfg,indent=2)+'\n')

    fig=plt.figure(figsize=(180/25.4,225/25.4))
    panel(fig,'a','Pan-cancer coverage',.025,.980)
    axa=fig.add_axes([.015,.716,.405,.247]);cov=coverage(axa,cancer)
    cov.to_csv(tables/'Figure1a_cancer_coverage.tsv',sep='\t',index=False)
    panel(fig,'b','Linked atlas resources',.445,.980)
    resources(fig.add_axes([.445,.736,.535,.217]),cfg)
    panel(fig,'c','Eight major cell compartments',.035,.690)
    fig.text(.064,.670,'Resolved TME reference · 1,974,309 cells',fontsize=6.3,color=MUTED,va='top')
    axc=fig.add_axes([.035,.319,.438,.334]);scatter(axc,xy,pcode,parent_colors,shuffle)
    c_labels=add_umap_labels(axc,xy,pcode,DISPLAY,fontsize=6.8,nudges={4:(-4,-3),3:(4,0)})
    panel(fig,'d','Forty-four resolved identities',.535,.690)
    fig.text(.564,.670,'Same cells and frozen coordinates',fontsize=6.3,color=MUTED,va='top')
    axd=fig.add_axes([.535,.319,.438,.334]);scatter(axd,xy,node,fine,shuffle)
    d_labels=add_umap_labels(axd,xy,node,[str(i) for i in range(1,45)],fontsize=6.1)
    assert len(c_labels)==8 and len(d_labels)==44
    panel(fig,'e','Native RNA support for the major compartments',.035,.294)
    axe=fig.add_axes([.126,.083,.735,.177]);xx,yy=np.meshgrid(np.arange(len(GENES)),np.arange(len(PARENTS)))
    frac=dots.fraction_positive.to_numpy().reshape(8,len(GENES))
    scale=dots.gene_scaled_mean.to_numpy().reshape(8,len(GENES))
    im=axe.scatter(xx,yy,s=frac*44,c=scale,cmap='Blues',vmin=0,vmax=1,lw=0)
    axe.set_xlim(-.7,len(GENES)-.3);axe.set_ylim(7.65,-.65)
    axe.set_xticks(range(len(GENES)),GENES,rotation=55,ha='right',fontsize=6.3)
    axe.set_yticks(range(8),DISPLAY,fontsize=6.8)
    axe.tick_params(axis='both',length=0,pad=4)
    for spine in axe.spines.values():spine.set_visible(False)
    for i in range(8):
        axe.plot([-.64,-.49],[i,i],color=PALETTE[i],lw=3,clip_on=False)
    cax=fig.add_axes([.910,.169,.013,.083]);cb=fig.colorbar(im,cax=cax,ticks=[0,.5,1]);cb.ax.tick_params(labelsize=5.7,length=2);cb.outline.set_linewidth(.4)
    fig.text(.895,.280,'Gene-scaled\nmean RNA',fontsize=5.8,va='top',color=MUTED)
    fig.text(.895,.150,'RNA-positive\ncells',fontsize=5.8,va='top',color=MUTED)
    for i,fraction in enumerate([.1,.5,1]):
        axe.scatter([],[],s=fraction*44,c='#7C8993',label=f'{int(fraction*100)}%')
    axe.legend(loc='upper left',bbox_to_anchor=(1.035,.27),frameon=False,labelspacing=.55,fontsize=5.8,handletextpad=.4,borderpad=0)
    fig.text(.035,.015,'c–d: one frozen global embedding.   e: qualified source–gene denominators; zeros retained.',fontsize=5.7,color=MUTED)
    # Keep text labels editable; every embedded point is represented in the scatter layers.
    for suffix in ['png','pdf','svg']:
        fig.savefig(out/('Figure1_atlas_redesigned.'+suffix),dpi=420)
        print('RENDERED',suffix,flush=True)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    outside=[]
    for text in list(fig.texts)+[t for ax in fig.axes for t in ax.texts]:
        if not text.get_text():continue
        box=text.get_window_extent(renderer)
        if box.xmin<0 or box.ymin<0 or box.xmax>fig.bbox.width or box.ymax>fig.bbox.height:
            outside.append(text.get_text())
    assert not outside,('Text outside figure',outside)
    label_audit=[]
    for pane,texts,ax in [('c',c_labels,axc),('d',d_labels,axd)]:
        boxes=[t.get_window_extent(renderer) for t in texts]
        collisions=[]
        for i in range(len(boxes)):
            for j in range(i+1,len(boxes)):
                if boxes[i].overlaps(boxes[j]):collisions.append([texts[i].get_text(),texts[j].get_text()])
        label_audit.append({'panel':pane,'labels':[t.get_text() for t in texts],'text_overlap_pairs':collisions})
    plt.close(fig)

    # Companion identity key has an independent canvas and frozen numbering.
    fig=plt.figure(figsize=(180/25.4,185/25.4))
    fig.text(.04,.955,'Companion key · 44 resolved TME identities',fontsize=12,fontweight='bold',color=INK)
    fig.text(.04,.922,'Numbers and colours match Figure 1d. Underlying identity labels are unchanged.',fontsize=7,color=MUTED)
    for row in key.itertuples():
        col=(row.number-1)//22;y=.868-((row.number-1)%22)*.036;x=.04+col*.50
        text=row.identity.replace('Endo_','EC ').replace('Str_','').replace('Mye_','').replace('_',' ')
        fig.text(x,y,'●',fontsize=9,color=row.color_hex,va='center')
        fig.text(x+.027,y,f'{row.number:02d}',fontsize=7.5,fontweight='bold',color=INK,va='center')
        fig.text(x+.069,y,text,fontsize=6.9,color=INK,va='center')
    for suffix in ['png','pdf','svg']:fig.savefig(out/('Figure1_companion_identity_key.'+suffix),dpi=330)
    plt.close(fig)
    receipt={'UTC':datetime.now(timezone.utc).isoformat(),'embedding_refit':False,'downsampling':False,
        'all_embedding_points':len(xy),'broad_compartments':8,'fine_identities':44,'reference_sources':70,
        'harmonised_cells':sum(cancer.values()),'harmonised_sources':len(registry),'harmonised_cancer_labels':len(cancer),
        'label_audit':label_audit,'input_sha256':{str(x.name):sha(x) for x in [a.registry,a.resources_json,a.node_map,embed,a.atlas_dir/'annotation_marker_statistics.tsv']},
        'source_code_sha256':sha(__file__),'outputs':{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/'render_receipt_PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('COMPLETE',json.dumps({'label_overlaps':sum(len(x['text_overlap_pairs']) for x in label_audit)}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--atlas-dir',type=Path,required=True);p.add_argument('--node-map',type=Path,required=True)
    p.add_argument('--registry',type=Path,required=True);p.add_argument('--resources-json',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args())
