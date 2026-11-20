#!/usr/bin/env python3
"""Render English R4 manuscript figures using frozen server-side evidence."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import FancyBboxPatch, Rectangle
import networkx as nx
import numpy as np
import pandas as pd

PARENTS=['B','Plasma','CD4T','CD8T','ILC','Myeloid','Endothelial','Stromal']
COLORS=dict(zip(PARENTS,['#0072B2','#56B4E9','#009E73','#117733','#CC79A7','#E69F00','#D55E00','#7C6DAD']))
CM=[f'CM{i:02d}' for i in range(1,16)]
TIP='Endo_tip_ESM1';PC='Str_Pericyte_RGS5';MAC='Mye_Mac_FOLR2'
CANCER_SHORT={'Hepatocellular_carcinoma':'HCC','Intrahepatic_cholangiocarcinoma':'iCCA',
 'NSCLC':'NSCLC','Colorectal_cancer':'Colorectal','Uveal_melanoma':'Uveal melanoma',
 'Basal_cell_carcinoma':'Basal cell','Osteosarcoma':'Osteosarcoma',
 'Pancreatic_ductal_adenocarcinoma':'PDAC','Renal_cell_carcinoma':'RCC','Breast_cancer':'Breast',
 'Nasopharyngeal_carcinoma':'Nasopharyngeal','Ovarian_cancer':'Ovarian','Gastric_cancer':'Gastric',
 'Diffuse_glioma':'Diffuse glioma','Head_and_neck_squamous_cell_carcinoma':'HNSCC',
 'Endometrial_cancer':'Endometrial','Prostate_cancer':'Prostate','Papillary_thyroid_carcinoma':'Papillary thyroid',
 'Acral_melanoma':'Acral melanoma','Cervical_cancer':'Cervical','Esophageal_cancer':'Esophageal',
 'Soft_tissue_sarcoma':'Soft-tissue sarcoma','Cutaneous_melanoma':'Cutaneous melanoma',
 'Bladder_cancer':'Bladder','Mesothelioma':'Mesothelioma','Pancreatic_neuroendocrine_tumor':'Pancreatic NET',
 'Renal_medullary_carcinoma':'Renal medullary'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
 'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,
 'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'savefig.facecolor':'white','figure.facecolor':'white','axes.titleweight':'bold'})


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def short(x):
    return str(x).replace('Endo_','EC ').replace('Str_','').replace('Mye_','').replace('_',' ')


def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=10)
    ax.text(-.08,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=13,va='top')


def clean(ax):
    ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values():s.set_visible(False)


def render(root, run, output):
    output.mkdir(parents=True,exist_ok=False)
    tables=output/'source_tables';tables.mkdir()
    source_sha={};figures=[]
    r4=root/'w2o_annotation_r4_covarnet_k15_20261115_v1';v=root/'CM_validation_20261118_v1'
    spat=root/'spatial_CoVarNet_R4_20261115_v1';atlas=run/'atlas'
    assert json.loads((atlas/'complete.json').read_text())['status'].startswith('PASS')
    def read(p,**kw):
        source_sha[str(p.relative_to(root))]=sha(p)
        return pd.read_csv(p,sep='\t',**kw)
    def export(name,data):data.to_csv(tables/(name+'.tsv'),sep='\t',index=False)
    def save(fig,name,description):
        for suffix in ['png','pdf','svg']:
            fig.savefig(output/(name+'.'+suffix),dpi=360,bbox_inches='tight',pad_inches=.13)
        plt.close(fig)
        figures.append({'id':name,'description':description,'files':{s:sha(output/(name+'.'+s)) for s in ['png','pdf','svg']}})
        (output/'progress.json').write_text(json.dumps({'finished':[x['id'] for x in figures],'UTC':datetime.now(timezone.utc).isoformat()}))
        print('RENDERED',name,flush=True)
    mapping=read(r4/'covarnet/node_mapping.tsv').set_index('node').majorCluster.to_dict()
    order=[n for p in PARENTS for n,c in mapping.items() if c==p]
    icount=read(atlas/'identity_counts.tsv').set_index('name').cells.reindex(order)
    source=read(atlas/'source_cancer_counts.tsv')
    audit=read(r4/'covarnet/unit_audit.tsv',dtype={'unit_id':str})
    review=read(r4/'covarnet/paper/cm_review.tsv')
    med=read(r4/'covarnet/paper/cm_missing_denominator_sensitivity.tsv')
    ranks=read(r4/'covarnet/paper/rank_survey.tsv')
    edges=read(r4/'covarnet/export/module_edges.tsv')
    member=read(r4/'covarnet/export/module_network_nodes.tsv')
    dot=read(atlas/'annotation_marker_statistics.tsv')
    ar=json.loads((atlas/'complete.json').read_text())

    # Figure 1: reproducible cohort flow and complete coverage in supplements.
    fig=plt.figure(figsize=(8.4,10));gs=fig.add_gridspec(2,2,height_ratios=[.85,1.3],hspace=.42,wspace=.54)
    ax=fig.add_subplot(gs[0,:]);ax.set_xlim(0,1);ax.set_ylim(0,1);clean(ax);panel(ax,'a','Study design and nested analysis populations')
    labels=[('Integrated TME','3,347,659 cells'),('R4 working labels',f"{ar['R4_working_cells']:,} cells"),
            ('Tier-qualified cells',f"{ar['R4_tier_main_cells']:,} cells"),('44 resolved identities',f"{ar['R4_44_identity_cells']:,} cells")]
    for j,(a,b) in enumerate(labels):
        x=.01+j*.25
        ax.add_patch(FancyBboxPatch((x,.59),.215,.27,boxstyle='round,pad=0.01',facecolor='#EDF3F7',edgecolor='#75889B'))
        ax.text(x+.1075,.77,a,ha='center',va='center',fontsize=8.1,fontweight='bold')
        ax.text(x+.1075,.66,b,ha='center',va='center',fontsize=8.3)
        if j<3:ax.annotate('',xy=(x+.245,.72),xytext=(x+.215,.72),arrowprops={'arrowstyle':'->','color':'#4B5965'})
    arms=[('Composition','1,196 units | 70 sources\n44 identities | 15 candidate CMs'),
          ('Spatial organization','4 Visium + 8 legacy ST sections\n6 breast imaging sections + 1 lung donor'),
          ('Molecular and clinical','CPDB and patient-level LR\nTCGA and three ICI cohorts')]
    for j,(a,b) in enumerate(arms):
        x=.02+j*.33
        ax.add_patch(FancyBboxPatch((x,.08),.29,.29,boxstyle='round,pad=0.01',facecolor='#F7F7F7',edgecolor='#A0A7AD'))
        ax.text(x+.145,.29,a,ha='center',fontsize=9,fontweight='bold')
        ax.text(x+.145,.18,b,ha='center',va='center',fontsize=6.5,linespacing=1.6)
        ax.annotate('',xy=(x+.145,.38),xytext=(.86,.56),arrowprops={'arrowstyle':'->','color':'#A0A7AD','lw':.8})
    a=source.groupby('source').cells.sum().sort_values(ascending=False);top=a.head(20)
    ax=fig.add_subplot(gs[1,0]);ax.barh(np.arange(len(top)),top.values/1000,color='#567F98');ax.invert_yaxis()
    ax.set_yticks(range(len(top)),top.index,fontsize=7);ax.set_xlabel('Resolved-identity cells (thousands)');panel(ax,'b','Twenty largest contributing sources')
    ax.text(0,-.19,f'All {len(a)} sources are retained; complete coverage in Figure S01.',transform=ax.transAxes,fontsize=7)
    can=source.groupby('cancer').cells.sum().sort_values(ascending=False)
    ax=fig.add_subplot(gs[1,1]);ax.barh(range(len(can)),can.values/1000,color='#8F849D');ax.invert_yaxis()
    ax.set_yticks(range(len(can)),[CANCER_SHORT.get(x,x) for x in can.index],fontsize=7);ax.set_xlabel('Resolved-identity cells (thousands)');panel(ax,'c','Current cancer-label coverage')
    save(fig,'Figure1_workflow_and_sources','Nested populations and source/cancer coverage; units are not universally verified patients.')
    export('Figure1_source_cancer_counts',source)

    # Figure 2: all points use the same frozen global embedding.
    with np.load(atlas/'atlas_embedding_PRIVATE.npz') as z:
        xy=z['xy'];nc=z['node'];nodes=z['nodes'].tolist();src=z['source'];canc=z['cancer'];sources=z['sources'].tolist();cancers=z['cancers'].tolist()
    assert nodes==order and len(xy)==1974309
    rng=np.random.default_rng(20261120);shuffle=rng.permutation(len(xy))
    pcode=np.array([PARENTS.index(mapping[n]) for n in nodes])[nc]
    parent_colors=np.array([colors.to_rgba(COLORS[p]) for p in PARENTS])
    fine_colors=[]
    for p in PARENTS:
        ns=[n for n in order if mapping[n]==p]
        base=np.array(colors.to_rgb(COLORS[p]))
        for k in range(len(ns)):fine_colors.append(tuple(base*(.6+.4*k/max(1,len(ns)-1))+np.ones(3)*(.4-.4*k/max(1,len(ns)-1))))
    fine_colors=np.clip(np.array(fine_colors),0,1)
    def umap(ax,codes,palette,numbers=False):
        ax.scatter(xy[shuffle,0],xy[shuffle,1],c=palette[codes[shuffle]],s=.12,linewidths=0,rasterized=True,alpha=.7)
        ax.set_xlabel('Frozen UMAP 1');ax.set_ylabel('Frozen UMAP 2');ax.set_xticks([]);ax.set_yticks([])
        if numbers:
            for j in range(len(nodes)):
                q=xy[nc==j];center=np.median(q,axis=0)
                ax.text(*center,str(j+1),ha='center',va='center',fontsize=5.4,fontweight='bold',bbox={'boxstyle':'round,pad=.12','fc':'white','ec':'none','alpha':.8})
    fig=plt.figure(figsize=(8.4,8.7));gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=.56,wspace=.28)
    ax=fig.add_subplot(gs[0,0]);umap(ax,pcode,parent_colors);panel(ax,'a','Eight parent compartments')
    handles=[plt.Line2D([],[],marker='o',ls='',color=COLORS[p],label=p,markersize=4) for p in PARENTS]
    ax.legend(handles=handles,fontsize=6.8,ncol=2,frameon=False,loc='lower center',bbox_to_anchor=(.5,-.31))
    ax=fig.add_subplot(gs[0,1]);umap(ax,nc,fine_colors,False);panel(ax,'b','Forty-four resolved identities')
    ax.text(.5,-.21,'Full identity color key: Figure S06',ha='center',transform=ax.transAxes,fontsize=7)
    selected=['CD3D','FOXP3','GZMK','FCGR3A','CD79A','MS4A1','MZB1','LYZ','FOLR2','VWF','CA4','ESM1','PI16','COL1A1','RGS5','MYH11']
    dd=dot.copy();dd['parent']=dd.identity.map(mapping)
    parent=dd.groupby(['parent','gene']).agg(n=('n_eligible_cells','sum'),positive=('n_positive','sum'))
    weighted=dd.assign(s=dd.mean_log1p_1e4*dd.n_eligible_cells).groupby(['parent','gene']).s.sum()
    means=np.array([[weighted.get((p,g),np.nan)/parent.loc[(p,g),'n'] if parent.loc[(p,g),'n']>0 else np.nan for g in selected] for p in PARENTS])
    frac=np.array([[parent.loc[(p,g),'positive']/parent.loc[(p,g),'n'] if parent.loc[(p,g),'n']>0 else np.nan for g in selected] for p in PARENTS])
    scaled=means/np.nanmax(means,axis=0)
    ax=fig.add_subplot(gs[1,:]);xx,yy=np.meshgrid(np.arange(len(selected)),np.arange(8));im=ax.scatter(xx,yy,s=frac*190,c=scaled,cmap='Blues',vmin=0,vmax=1,edgecolors='none')
    ax.set_xticks(range(len(selected)),selected,rotation=55,ha='right',fontsize=8);ax.set_yticks(range(8),PARENTS);ax.invert_yaxis();ax.set_xlim(-.8,len(selected)-.2)
    panel(ax,'c','Native RNA support across parent compartments')
    fig.colorbar(im,ax=ax,fraction=.018,pad=.025,label='Gene-scaled mean log RNA')
    for fraction,label in [(.1,'10%'),(.5,'50%'),(1,'100%')]:ax.scatter([],[],s=fraction*190,color='#7D8892',label=label)
    ax.legend(title='RNA-positive cells',loc='upper center',bbox_to_anchor=(.5,-.29),ncol=3,frameon=False,fontsize=7,title_fontsize=7)
    save(fig,'Figure2_unified_atlas','Frozen global embedding subset and full-native RNA dotplot; source analysis-feature-qualified denominators.')

    # Full fine-identity dotplot is retained at readable size.
    genes=list(dict.fromkeys([n.split('_')[-1] for n in nodes]+dot.gene.unique().tolist()));stats=dot.set_index(['identity','gene'])
    mean=np.array([[stats.loc[(n,g),'mean_log1p_1e4'] for g in genes] for n in nodes]);fraction=np.array([[stats.loc[(n,g),'fraction_positive'] for g in genes] for n in nodes])
    scaled=mean/np.maximum(np.nanmax(mean,axis=0),1e-12)
    fig,ax=plt.subplots(figsize=(11.8,10.4));fig.subplots_adjust(left=.27,bottom=.16,right=.94,top=.95)
    xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(44));im=ax.scatter(xx,yy,s=fraction*105,c=scaled,cmap='Blues',vmin=0,vmax=1,linewidths=0)
    ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=8);ax.set_yticks(range(44),[f'{j+1:02d}  {short(n)}' for j,n in enumerate(nodes)],fontsize=7.8);ax.invert_yaxis();ax.set_xlim(-1,len(genes))
    y=0
    for p in PARENTS:
        n=sum(mapping[x]==p for x in nodes);ax.axhline(y-.5,color='#D5DADF',lw=.6);ax.add_patch(Rectangle((-.9,y-.5),.25,n,color=COLORS[p],clip_on=False));y+=n
    bad=~np.isfinite(mean);ax.scatter(xx[bad],yy[bad],marker='x',color='#777777',s=13)
    fig.colorbar(im,ax=ax,fraction=.015,pad=.02,label='Gene-scaled mean log RNA')
    ax.set_title('R4 annotation marker support across all 44 identities',loc='left',fontsize=13)
    fig.text(.27,.035,'Dot area: RNA-positive fraction. Means include zeros within feature-qualified sources.\nGene-specific source masks define eligibility. Complete denominators accompany the figure.',fontsize=8)
    save(fig,'FigureS02_annotation_dotplot','All 44 identities and their frozen naming-marker genes, independent normalization and explicit per-gene coverage.')
    export('FigureS02_annotation_marker_statistics',dot)
    export('FigureS02_identity_numbers',pd.DataFrame({'number':range(1,45),'identity':nodes,'parent':[mapping[n] for n in nodes],'cells':icount.values}))

    # Figure 3: full module roster, rank and missingness sensitivity.
    fig=plt.figure(figsize=(8.4,9.5));gs=fig.add_gridspec(2,2,hspace=.45,wspace=.4,height_ratios=[1.2,1])
    roster=pd.crosstab(member.subCluster,member.cm).reindex(index=nodes,columns=CM,fill_value=0)
    ax=fig.add_subplot(gs[0,0]);ax.imshow(roster,aspect='auto',cmap=colors.ListedColormap(['#F0F2F4','#476C85']),vmin=0,vmax=1,interpolation='nearest')
    ax.set_yticks(range(44),[str(i+1) for i in range(44)],fontsize=5.6);ax.set_xticks(range(15),[str(i) for i in range(1,16)],fontsize=7);ax.set_xlabel('R4 CM');ax.set_ylabel('Identity number (Figure S02)');panel(ax,'a','Overlapping module membership')
    ax=fig.add_subplot(gs[0,1]);ax.plot(ranks['rank'],ranks.cophenetic,'o-',color='#567F98',ms=3);ax.axvline(15,color='#B75B40',ls='--');ax.set_xticks([2,5,10,15,20]);ax.set_xlabel('Factorization rank');ax.set_ylabel('Cophenetic coefficient');panel(ax,'b','Frozen rank survey')
    ax=fig.add_subplot(gs[1,0]);
    for i,cm in enumerate(CM):
        vals=med.loc[med.cm==cm,'cosine'].to_numpy();ax.scatter(np.full(len(vals),i)+np.linspace(-.13,.13,len(vals)),vals,s=16,color='#8B7E9B');ax.plot([i-.23,i+.23],[np.median(vals)]*2,color='#353A40',lw=1.2)
    ax.set_ylim(-.02,1.03);ax.set_xticks(range(15),[str(i) for i in range(1,16)]);ax.set_xlabel('R4 CM');ax.set_ylabel('Matched cosine');panel(ax,'c','Missing-denominator sensitivity')
    ax=fig.add_subplot(gs[1,1]);parents=[f'n_classes_ge20'];freq=audit.n_classes_ge20.value_counts().reindex(range(9),fill_value=0);ax.bar(freq.index,freq,color='#B68C63');ax.set_xlabel('Parent classes with at least 20 cells');ax.set_ylabel('Analysis units');ax.set_xticks(range(9));panel(ax,'d','Capture availability across eight parents')
    fig.text(.13,.015,'Sensitivity fits evaluate the frozen discovery population under explicit denominator assumptions.',fontsize=7.5)
    save(fig,'Figure3_modules_and_robustness','All candidate modules, frozen rank scan and denominator sensitivity; no claim of universal module stability.')
    export('Figure3_module_review',review);export('Figure3_missing_denominator_fits',med)
    fig,axes=plt.subplots(3,5,figsize=(15,10));
    for ax,cm in zip(axes.flat,CM):
        frame=edges[edges.cm==cm];g=nx.from_pandas_edgelist(frame,'subCluster1','subCluster2');g.add_nodes_from(member.loc[member.cm==cm,'subCluster'])
        pos=nx.spring_layout(g,seed=20261120,k=1.2)
        nx.draw_networkx_edges(g,pos,ax=ax,edge_color='#A7ADB4',width=1,alpha=.8)
        nx.draw_networkx_nodes(g,pos,ax=ax,node_size=105,node_color=[COLORS[mapping[n]] for n in g],edgecolors='white',linewidths=.6)
        nx.draw_networkx_labels(g,pos,ax=ax,labels={n:str(nodes.index(n)+1) for n in g},font_size=6,font_color='white',font_weight='bold')
        ax.set_title('R4 '+cm,fontsize=11);clean(ax)
    fig.subplots_adjust(wspace=.2,hspace=.25);fig.suptitle('All R4 module networks | identity numbers match Figure S02',fontsize=15)
    save(fig,'FigureS03_all_module_networks','Frozen topology, constant node sizes and parent colors; no model loadings exported as node sizes.')

    # Figure 4: original-coordinate Visium maps, all12 section effects and all N.
    technical=read(spat/'review/Visium_technical.tsv')
    density=read(spat/'review/Visium_density_sensitivity.tsv')
    ncore=read(v/'spatial_melanoma/N_sensitivity_all24_v1/core_edges_all24.tsv')
    fig=plt.figure(figsize=(8.4,10.1));gs=fig.add_gridspec(3,4,height_ratios=[1,1.25,1.05],hspace=.65,wspace=.58)
    for j,row in technical.sort_values('section').reset_index(drop=True).iterrows():
        xy0=read(spat/f'sections/{row.section}/coordinates_PRIVATE.tsv.gz',index_col=0)
        ab=read(spat/f'sections/{row.section}/means_abundance_PRIVATE.tsv.gz',index_col=0).loc[xy0.index]
        a=ab[TIP].to_numpy();p=ab[PC].to_numpy()
        a=np.clip(a/max(np.quantile(a,.99),1e-9),0,1);p=np.clip(p/max(np.quantile(p,.99),1e-9),0,1)
        rgb=np.column_stack([1-.6*p,1-.75*np.maximum(a,p),1-.6*a])
        ax=fig.add_subplot(gs[0,j]);ax.scatter(xy0.x,xy0.y,c=rgb,s=6,linewidths=0,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();clean(ax)
        ax.set_title(f'{row.dataset}\n{row.section}',fontsize=8)
        if j==0:ax.text(-.2,1.35,'a',transform=ax.transAxes,fontweight='bold',fontsize=13)
    fig.text(.13,.66,'Native capture coordinates | inferred tip EC (red), pericyte (blue), overlap (purple)\nWithin-section 99th-percentile display scaling; no H&E registration or contact inference.',fontsize=7)
    ax=fig.add_subplot(gs[1,:2]);core=density[(density.cm=='CM12')&(((density.a==TIP)&(density.b==PC))|((density.b==TIP)&(density.a==PC)))].copy()
    vals=core.groupby('section').partial_rank_correlation_control_total_and_epithelial_fraction.mean()
    ax.scatter(vals.values,range(len(vals)),color='#AD604D',s=35);ax.set_yticks(range(len(vals)),vals.index);ax.axvline(0,color='#B7BDC3',lw=.8);ax.set_xlim(-.1,.8);ax.set_xlabel('Composition-adjusted rank correlation');panel(ax,'b','Visium vascular core')
    ax=fig.add_subplot(gs[1,2:]);sel=ncore[(ncore.cm=='CM12')&(ncore.posterior=='means')&(((ncore.node1==TIP)&(ncore.node2==PC))|((ncore.node2==TIP)&(ncore.node1==PC)))]
    for j,n in enumerate([5,30,60]):
        q=sel[sel.assumed_N==n].sort_values('sample');x=j+np.linspace(-.18,.18,len(q));ax.scatter(x,q.partial_rank_total_melanoma_other_nuisance,color=['#93A7BB','#AD604D','#7D7896'][j],s=22)
    ax.axhline(0,color='#B7BDC3',lw=.8);ax.set_xticks(range(3),['N = 5','N = 30','N = 60']);ax.set_ylim(-.35,.8);ax.set_ylabel('Composition-adjusted rank correlation');panel(ax,'c','All eight legacy ST sections')
    ax=fig.add_subplot(gs[2,:2]);pairs=[(TIP,PC,'Tip EC - pericyte'),(TIP,MAC,'Tip EC - FOLR2 macrophage'),(MAC,PC,'FOLR2 macrophage - pericyte')]
    rows=[]
    for i,(a,b,label) in enumerate(pairs):
        q=ncore[(ncore.cm=='CM12')&(ncore.posterior=='means')&(ncore.assumed_N==30)&(((ncore.node1==a)&(ncore.node2==b))|((ncore.node2==a)&(ncore.node1==b)))]
        yy=i+np.linspace(-.17,.17,len(q));ax.scatter(q.partial_rank_total_melanoma_other_nuisance,yy,s=22,color=['#AD604D','#C39B5F','#7D7896'][i]);rows.extend(q.to_dict('records'))
    ax.axvline(0,color='#B7BDC3',lw=.8);ax.set_yticks(range(3),[x[2] for x in pairs],fontsize=7);ax.set_xlim(-.25,.8);ax.set_xlabel('Composition-adjusted rank correlation');panel(ax,'d','Core and macrophage-associated edges')
    ax=fig.add_subplot(gs[2,2:]);q=sel[sel.assumed_N==30];ax.scatter(q.partial_rank_total_melanoma_other_nuisance,q.parent_fraction_spearman,c='#6B8796',s=25);ax.axhline(0,color='#B7BDC3',lw=.8);ax.axvline(0,color='#B7BDC3',lw=.8);ax.set_xlabel('Composition-adjusted correlation');ax.set_ylabel('Within-parent fraction correlation');panel(ax,'e','Alternative abundance denominator')
    save(fig,'Figure4_spatial_organization','Four Visium maps and full melanoma parameter sensitivity; section-level descriptive effects, not independent-patient tests.')
    export('Figure4_Visium_density_effects',density);export('Figure4_legacy_ST_core_edges_all24',ncore)

    # Figure 5: actual imaging coordinates and all primary comparison directions.
    imaging=read(v/'imaging_lung_core_v1/joint_imaging_tests.tsv')
    fig=plt.figure(figsize=(8.4,9.2));gs=fig.add_gridspec(2,2,hspace=.4,wspace=.4,height_ratios=[1,1.3])
    for j,(rel,label) in enumerate([('imaging_breast_core_cells_v1/S01_core_PRIVATE.npz','Breast section S01: broad mural proxy'),('imaging_lung_core_v1/EC2_nucleated_cells_PRIVATE.npz','Lung: auxiliary-marker proxies')]):
        path=v/rel;source_sha[str(path.relative_to(root))]=sha(path)
        with np.load(path,allow_pickle=True) as a:
            xx=a['xy'];ec=a['ec'] if 'ec' in a.files else a['EC'];pc=a['mural'] if 'mural' in a.files else a['PC'];ep=ec&(a['ESM1']>=1);pp=pc&(a['RGS5']>=1)
        ax=fig.add_subplot(gs[0,j]);ax.scatter(xx[ec|pc,0],xx[ec|pc,1],s=.35,c='#CED4DA',rasterized=True,linewidths=0)
        ax.scatter(xx[ep,0],xx[ep,1],s=2,c='#B9564A',label='ESM1+ EC proxy',rasterized=True,linewidths=0);ax.scatter(xx[pp,0],xx[pp,1],s=2,c='#5879A8',label='RGS5+ mural proxy',rasterized=True,linewidths=0)
        ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('Native x (micrometres)');ax.set_ylabel('Native y (micrometres)');panel(ax,chr(97+j),label);ax.legend(fontsize=6.3,frameon=False,loc='upper right',markerscale=2)
    primary=imaging[(imaging.threshold==1)&(imaging.tile_um==200)&(imaging.radius_um==30)]
    for j,direction in enumerate(['ESM1_EC_to_RGS5_mural','RGS5_mural_to_ESM1_EC']):
        ax=fig.add_subplot(gs[1,j]);q=primary[(primary.platform_cohort=='breast')&(primary.current_state=='pc_low_MYH11')&(primary.direction==direction)].copy()
        lungdir='EC_to_PC' if j==0 else 'PC_to_EC'
        lung=primary[(primary.platform_cohort=='lung')&(primary.current_state=='EC2_nucleated')&(primary.direction==lungdir)].copy()
        q=pd.concat([q,lung],ignore_index=True);labels=[('Breast '+str(r.section)) if r.platform_cohort=='breast' else 'Lung (one donor)' for _,r in q.iterrows()]
        assert len(q)==7,(direction,len(q))
        y=np.arange(len(q));sig=q.q_joint_imaging<.05
        for i,(_,r) in enumerate(q.iterrows()):
            if not np.isfinite(r.q_joint_imaging):
                ax.text(.08,i,'Not testable',transform=ax.get_yaxis_transform(),ha='left',fontsize=7,color='#777777')
                continue
            ax.plot([0,r.effect_fraction],[i,i],color='#D5DADE',lw=1);ax.scatter(r.effect_fraction,i,s=40,facecolor='#567F98' if sig.iloc[i] else 'white',edgecolor='#567F98',zorder=3)
            ax.text(.97,i,f'q={r.q_joint_imaging:.2g}',transform=ax.get_yaxis_transform(),ha='right',fontsize=6.8)
        ax.axvline(0,color='#AAB3BB',ls='--',lw=.8);ax.set_yticks(y,labels,fontsize=7);ax.invert_yaxis();ax.set_xlabel('Observed - null expected positive fraction');panel(ax,chr(99+j),'EC -> PC proxy' if j==0 else 'PC -> EC proxy')
    fig.text(.12,.015,'Primary geometry: 30 um radius, 200 um spatial blocks, >=1 RNA. Filled: joint q < 0.05; open: joint q >= 0.05.\nBlank rows are ineligible. Direction denotes source-parent neighborhood geometry; effects use the frozen conditional null.',fontsize=7)
    save(fig,'Figure5_imaging_neighborhoods','Native imaging overview and all six breast plus one lung primary comparisons; proxy identities and unequal inference units explicit.')
    export('Figure5_complete_imaging_family',imaging)

    # Figure 6: communication specificity and qualified patient contrasts.
    overlap=read(v/'communication_joint15_v2/core_overlap_summary.tsv')
    markers=read(v/'highlow_joint15_v1/fixed_marker_results.tsv')
    pat=read(v/'differential_LR_three_source_joint_v1/source_module_summary.tsv')
    fig=plt.figure(figsize=(8.4,9.6));gs=fig.add_gridspec(2,2,hspace=.52,wspace=.48)
    ax=fig.add_subplot(gs[0,:]);q=overlap[overlap.core_module.isin(['CM12','CM01'])].copy();y=np.arange(len(q));ax.barh(y,q.focal_shared_with_other,color='#63859B',label='Also detected in another CMT');ax.barh(y,q.focal_only_significant,left=q.focal_shared_with_other,color='#C89360',label='Detected only in focal CMT')
    labels=[f"{r.core_module}  {short(r.source_identity)} -> {short(r.target_identity)}" for _,r in q.iterrows()];ax.set_yticks(y,labels,fontsize=7);ax.invert_yaxis();ax.set_xlabel('Joint-q-positive CPDB event combinations');ax.legend(frameon=False,fontsize=7,loc='lower right');panel(ax,'a','Most focal communication detections are shared')
    ax=fig.add_subplot(gs[1,0]);fixed=[(TIP,'ESM1'),(TIP,'KDR'),(TIP,'PDGFB'),(TIP,'VEGFA'),(PC,'RGS5'),(PC,'PDGFRB'),(PC,'CXCL12')];labels=[]
    for i,(node,gene) in enumerate(fixed):
        rr=markers[(markers.cm=='CM12')&(markers.node==node)&(markers.gene==gene)&(markers.family=='unit')];assert len(rr)==1
        a=rr.iloc[0];labels.append(gene);ax.errorbar(a.high_minus_low_mean_log_expression,i,xerr=1.96*a.HC3_se,fmt='o',color='#7C6DAD',markersize=4,capsize=2);ax.text(.97,i,f'q={a.q_joint15_model:.2g}',transform=ax.get_yaxis_transform(),ha='right',fontsize=6.7)
    ax.axvline(0,color='#AAB3BB',lw=.8);ax.set_yticks(range(len(fixed)),labels);ax.invert_yaxis();ax.set_xlabel('Adjusted mean log-RNA difference (95% interval)');panel(ax,'b','CM12 same-identity High/Low RNA')
    ax=fig.add_subplot(gs[1,1]);p=pat.groupby('source').agg(tests=('tested_variable_scores','sum'),q05=('q05','sum'));ax.barh(range(len(p)),p.tests,color='#8E9FA9');ax.set_yticks(range(len(p)),p.index,fontsize=7);ax.set_xlabel('Variable, patient-qualified LR hypotheses');panel(ax,'c','Three source-verified patient designs')
    for i,(_,rr) in enumerate(p.iterrows()):ax.text(rr.tests+p.tests.max()*.03,i,f"{int(rr.tests):,}",ha='left',va='center',fontsize=8,color='#333333')
    ax.set_xlim(0,p.tests.max()*1.23)
    fig.text(.13,.02,'CPDB tests cell-label specificity, not High/Low differential communication. Unit RNA intervals are HC3 Wald intervals.\nPatients are verified only in the three-source LR design; unavailable endpoints are not negative evidence.',fontsize=7)
    save(fig,'Figure6_communication_evidence','Focal/shared CPDB events, fixed marker unit models and all three patient-qualified differential LR source summaries.')
    export('Figure6_CPDB_shared_events',overlap);export('Figure6_fixed_markers',markers);export('Figure6_patient_LR_coverage',pat)

    # Figure 7: all15 immune associations plus frozen clinical candidates/ICI.
    immune=read(v/'bulk_R4_candidates/TCGA_leukocyte_purity_v3_nonredacted/primary_branch_purity_age_sex.tsv')
    surv=read(v/'bulk_R4_candidates/TCGA_survival_v1_nonredacted/piecewise_one_year_v1/piecewise_associations.tsv')
    ici=read(v/'bulk_R4_candidates/rose_response_v1/three_cohort_descriptive_comparison.tsv')
    fig=plt.figure(figsize=(8.4,9.3));gs=fig.add_gridspec(2,2,height_ratios=[1.2,1],hspace=.65,wspace=.5)
    q=immune[(immune.branch=='original_network_members')&(immune.model=='purity_age_sex')&(immune.policy=='LF_single')];cancers0=sorted(q.cancer.unique());mat=q.pivot(index='cm',columns='cancer',values='partial_rho').reindex(index=CM,columns=cancers0);qv=q.pivot(index='cm',columns='cancer',values='q_joint_with_unadjusted').reindex_like(mat)
    ax=fig.add_subplot(gs[0,:]);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#D8DDE1');im=ax.imshow(mat,aspect='auto',cmap=cmap,vmin=-1,vmax=1);ax.set_yticks(range(15),CM,fontsize=7);ax.set_xticks(range(len(cancers0)),cancers0,rotation=90,fontsize=6.4)
    for i,j in zip(*np.where(qv.to_numpy()<.05)):ax.text(j,i,'*',ha='center',va='center',fontsize=7)
    fig.colorbar(im,ax=ax,fraction=.018,pad=.02,label='Adjusted rank correlation');panel(ax,'a','Cancer-specific association with leukocyte fraction')
    ax=fig.add_subplot(gs[1,0]);s=surv[(surv.cm=='CM13')&(surv.cancer=='LIHC')&(surv.endpoint=='OS')&(surv.branch=='original_network_members')].set_index('model').reindex(['base','stage','composition','stage_composition']);assert s.status.eq('FIT_OK').all()
    for j,period,color in [(0,'early','#A9624D'),(1,'late','#547F9D')]:
        for i,(_,rr) in enumerate(s.iterrows()):
            y=i+(.12 if j else -.12);a=rr[f'HR_{period}'];lo=rr[f'CI_{period}_low'];hi=rr[f'CI_{period}_high'];ax.errorbar(a,y,xerr=[[a-lo],[hi-a]],fmt='o',color=color,ms=4,capsize=2,label=period if i==0 else None)
    ax.axvline(1,color='#AAB3BB',ls='--',lw=.8);ax.set_xscale('log');ax.set_yticks(range(4),['Age + sex','+ stage','+ purity / LF','+ stage + purity / LF'],fontsize=7);ax.invert_yaxis();ax.set_xlabel('Hazard ratio (95% interval)');ax.legend(frameon=False,fontsize=7);panel(ax,'b','Exploratory CM13 / LIHC overall survival')
    ax=fig.add_subplot(gs[1,1]);a=ici[ici.cm=='CM12'].copy();assert len(a)==3
    for i,(_,rr) in enumerate(a.iterrows()):ax.errorbar(rr.auc,i,xerr=[[rr.auc-rr.CI_low],[rr.CI_high-rr.auc]],fmt='o',color='#7C6DAD',capsize=2,ms=5);ax.text(.98,i,f'n={int(rr.patients)}',transform=ax.get_yaxis_transform(),ha='right',fontsize=7)
    ax.axvline(.5,color='#AAB3BB',ls='--',lw=.8);ax.set_xlim(0,1);ax.set_yticks(range(3),a.cohort);ax.invert_yaxis();ax.set_xlabel('Fixed-direction response AUC (95% descriptive interval)');panel(ax,'c','CM12: three ICI cohorts')
    fig.text(.12,.015,'A: original-network signature; purity/age/sex adjustment; stars use the 5,400-test joint family. Gray: unavailable.\nB: fixed one-year split; four models share a cohort. C: no joint-significant ICI response association; Rose timing unknown.',fontsize=7)
    save(fig,'Figure7_cancer_and_clinical_context','All15 cancer-specific immune associations and exploratory survival/response contexts; no validated predictor claim.')
    export('Figure7_all15_immune_associations',immune);export('Figure7_all_piecewise_models',surv);export('Figure7_ICI_descriptive_intervals',ici)

    # Supplements: complete coverage, source/cancer UMAP, model and edge sensitivities.
    source_matrix=source.pivot_table(index='source',columns='cancer',values='cells',fill_value=0).reindex(index=source.groupby('source').cells.sum().sort_values(ascending=False).index)
    fig,ax=plt.subplots(figsize=(11,14));im=ax.imshow(np.log10(1+source_matrix),aspect='auto',cmap='Blues');ax.set_yticks(range(len(source_matrix)),source_matrix.index,fontsize=7);ax.set_xticks(range(len(source_matrix.columns)),source_matrix.columns,rotation=90,fontsize=7);fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label='log10(1 + cells)');ax.set_title('Complete source-by-cancer coverage in the 44-identity view',loc='left');save(fig,'FigureS01_complete_source_coverage','Every source and current cancer label; cells, not patient counts.')
    fig,axes=plt.subplots(1,2,figsize=(12,5.2))
    for ax,codes,labels,title in [(axes[0],src,sources,'Source'),(axes[1],canc,cancers,'Cancer label')]:
        most=np.argsort(np.bincount(codes,minlength=len(labels)))[::-1][:10]
        remap=np.full(len(labels),10);remap[most]=np.arange(10);pal=np.vstack([plt.get_cmap('tab20')(np.arange(0,20,2)),colors.to_rgba('#D9DDE1')])
        umap(ax,remap[codes],pal[:,:3]);ax.set_title('Ten largest '+title.lower()+' categories',loc='left')
        key=[labels[k] if title=='Source' else CANCER_SHORT.get(labels[k],labels[k]) for k in most]+['Other categories']
        handles=[plt.Line2D([],[],marker='o',ls='',color=pal[j],label=t,markersize=4) for j,t in enumerate(key)]
        ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.1),ncol=2,frameon=False,fontsize=7)
    save(fig,'FigureS04_source_and_cancer_embedding','Same cells and coordinates; ten largest categories by cell count and all other cells in gray, not a batch-removal test.')
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,metric,title in zip(axes,['raw_spearman','partial_rank_total_melanoma_other_nuisance','parent_fraction_spearman'],['Raw correlation','Composition-adjusted','Within-parent fraction']):
        for i,(a,b,label) in enumerate(pairs):
            for j,n in enumerate([5,30,60]):
                q=ncore[(ncore.posterior=='means')&(ncore.assumed_N==n)&(ncore.cm=='CM12')&(((ncore.node1==a)&(ncore.node2==b))|((ncore.node2==a)&(ncore.node1==b)))]
                ax.scatter(i+np.linspace(-.25,.25,len(q)),q[metric],s=13,alpha=.65,color=['#93A7BB','#AD604D','#7D7896'][j],marker=['o','s','^'][j],label=f'N = {n}' if i==0 else None)
        ax.axhline(0,color='#BBBBBB',lw=.8);ax.set_xticks(range(3),['EC-PC','EC-Mac','Mac-PC']);ax.set_ylim(-1,1);ax.set_title(title,fontsize=10)
    axes[0].legend(frameon=False,fontsize=8,loc='lower left')
    save(fig,'FigureS05_spatial_denominator_sensitivity','All melanoma core-edge N settings; repeated settings not independent sections or patients.')
    lr=read(v/'communication_joint15_v2/core_edges_all_CMT.tsv.gz')
    fixed=['VEGFA_KDR','VEGFA_FLT1','PDGFB_PDGFRB','CXCL12_CXCR4','ANGPT1_TEK','COL1A1_CD44','LGALS3_ENG','APP_CD74','SPP1_CD44','GRN_TNFRSF1B','LGALS1_PTPRC','CCL5_CCR1']
    focal=lr[lr.core_module=='CM12'];pairs0=focal[['source_identity','target_identity']].drop_duplicates().apply(tuple,axis=1).tolist();arr=np.full((len(pairs0),len(fixed)),np.nan)
    for i,(a,b) in enumerate(pairs0):
        for j,name in enumerate(fixed):
            q=focal[(focal.source_identity==a)&(focal.target_identity==b)&(focal.interacting_pair==name)]
            if len(q):arr[i,j]=q.loc[q.q_joint15<.05,'CMT'].nunique()
    fig,ax=plt.subplots(figsize=(11,4.2));cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#D8DDE1');im=ax.imshow(arr,aspect='auto',vmin=0,vmax=15,cmap=cmap);ax.set_xticks(range(len(fixed)),[s.replace('_',' - ') for s in fixed],rotation=55,ha='right',fontsize=8);ax.set_yticks(range(len(pairs0)),[short(a)+' -> '+short(b) for a,b in pairs0],fontsize=7);fig.colorbar(im,ax=ax,fraction=.02,pad=.02,label='CMTs with joint-q-positive detection');ax.set_title('Prespecified LR records across CM12 core endpoint pairs',loc='left');fig.text(.13,-.22,'Gray: not present in eligible exported records; zero: evaluated but no joint-positive CMT. Detection is not differential signaling.',fontsize=8);save(fig,'FigureS06_LR_record_coverage','Fixed candidate record coverage without effect-directed selection; empty exports retained as unavailable.')
    export('FigureS06_fixed_LR_record_coverage',pd.DataFrame(arr,index=[a+' -> '+b for a,b in pairs0],columns=fixed).rename_axis('ordered_pair').reset_index())
    fig,axes=plt.subplots(2,2,figsize=(10,8));base=imaging[(imaging.threshold==1)&(imaging.tile_um==200)&(imaging.radius_um==30)]
    for ax,state in zip(axes.flat,['original_broad_mural','pc_low_MYH11','MYH11_contractile','EC2_nucleated']):
        q=base[base.current_state==state].reset_index(drop=True)
        valid=q.q_joint_imaging.notna();ax.scatter(q.loc[valid,'effect_fraction'],np.flatnonzero(valid),s=22,c=np.where(q.loc[valid,'q_joint_imaging']<.05,'#567F98','#CCCCCC'));
        for row in np.flatnonzero(~valid):ax.text(.03,row,'Not testable',transform=ax.get_yaxis_transform(),fontsize=6,color='#777777')
        ax.axvline(0,color='#BBBBBB',lw=.8);ax.set_yticks(range(len(q)),[str(r.section)+' '+r.direction for _,r in q.iterrows()],fontsize=6);ax.set_title(state,fontsize=10);ax.set_xlabel('Observed - null expected fraction')
    fig.subplots_adjust(wspace=.9,hspace=.4);fig.text(.12,-.015,'Blue: joint q < 0.05; gray: tested, not significant; blank rows: not testable.',fontsize=8);save(fig,'FigureS07_imaging_identity_sensitivity','All primary-geometry records under each identity rule; joint family values preserved.')
    const=read(v/'bulk_R4_candidates/TCGA_survival_v1_nonredacted/piecewise_one_year_v1/prior_associations_joint_q.tsv')
    fig,axes=plt.subplots(1,2,figsize=(13,6))
    for ax,endpoint in zip(axes,['OS','PFI']):
        q=const[(const.branch=='original_network_members')&(const.model=='base')&(const.endpoint==endpoint)]
        mat=q.pivot(index='cm',columns='cancer',values='HR').reindex(index=CM);qq=q.pivot(index='cm',columns='cancer',values='q_joint_piecewise').reindex_like(mat)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#D8DDE1')
        im=ax.imshow(np.log2(mat),cmap=cmap,vmin=-1.5,vmax=1.5,aspect='auto');ax.set_xticks(range(len(mat.columns)),mat.columns,rotation=90,fontsize=6);ax.set_yticks(range(15),CM,fontsize=7)
        for i,j in zip(*np.where(qq.to_numpy()<.05)):ax.text(j,i,'*',ha='center',va='center',fontsize=7)
        ax.set_title(endpoint+' | age/sex model',fontsize=11)
    fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.015,pad=.015,label='log2 hazard ratio',extend='both');save(fig,'FigureS08_all15_survival_context','All15 base-model survival estimates with updated extended-family q values; unavailable models remain missing.')
    fig,ax=plt.subplots(figsize=(7.5,4));status=imaging.groupby(['platform_cohort','status']).size().unstack(fill_value=0);status.plot.bar(stacked=True,ax=ax,color=['#567F98','#BFC5CA','#C89365']);ax.set_xlabel('Imaging cohort');ax.set_ylabel('Planned settings');ax.set_title('Imaging test eligibility, including insufficient settings',loc='left');ax.legend(frameon=False,fontsize=7);save(fig,'FigureS09_imaging_coverage','480 planned settings, preserving all insufficient rows and the 466-test inference family.')
    fig=plt.figure(figsize=(12.5,9));ax=fig.add_axes([.06,.12,.58,.78]);umap(ax,nc,fine_colors,False);ax.set_title('All 44 identities on the frozen global UMAP',loc='left',fontsize=13)
    handles=[plt.Line2D([],[],marker='o',ls='',color=fine_colors[j],label=f'{j+1:02d}  {short(n)}',markersize=5) for j,n in enumerate(nodes)]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.65,.5),frameon=False,fontsize=8.3,labelspacing=.65)
    save(fig,'FigureS10_full_identity_UMAP','Full fine-identity color key with unchanged coordinates and complete resolved-identity population.')

    manifest={'status':'FIGURES_RENDERED_PENDING_VISUAL_AND_INDEPENDENT_AUDIT','UTC':datetime.now(timezone.utc).isoformat(),
              'R4_version':'2026-11-15_R4_K15','figures':figures,'input_sha256':source_sha,
              'render_code_sha256':sha(Path(__file__)),'new_inferential_tests':False,
              'notes':['All scientific source tables stay on server.','English figures shared by English and Chinese manuscripts.',
                       'No private cell arrays or patients exported as figure-source spreadsheets.','UMAP is frozen; dotplot masks have explicit per-point denominators.']}
    (output/'figure_manifest_PRIVATE.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('ALL_FIGURES_RENDERED',len(figures),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();render(a.root,a.run,a.output)
