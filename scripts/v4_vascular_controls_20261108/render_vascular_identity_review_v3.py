"""Native marker review figures: every marker is paired with expected identity."""
from pathlib import Path
import sys,json,hashlib,os
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
r=Path(sys.argv[1]);t=pd.read_csv(r/'author_marker_review.tsv',sep='\t');t=t[t.level.eq('minor')&t.n_cells.ge(20)];out=r/os.environ.get('VASCULAR_FIGURE_SUBDIR','figures');out.mkdir(exist_ok=False)
for name,groups,pattern in [('immune',['B','T','NK_anchor','NK_effector'],'^B cells|^T cells|^NK|^NKT|Cycling T|Plasmablast'),('vascular_stromal',['EC','Lymphatic','Mural','Fibroblast','Venular_context','Arterial_context','Capillary_context'],'Endothelial|PVL|CAFs')]:
 z=t[t.author_label.str.contains(pattern,case=False,regex=True)];rows=sorted(z.author_label.unique());cols=[];expected={}
 for group in groups:
  for g in t.loc[t.expected_identity_or_axis.eq(group),'gene'].unique():
   if g not in cols:cols.append(g);expected[g]=group
 if not rows or not cols:continue
 fig,ax=plt.subplots(figsize=(max(10,len(cols)*.44),max(4,len(rows)*.3+2.4)));
 if 'SYNTHETIC' in r.name:fig.suptitle('SYNTHETIC SOFTWARE TEST ONLY',color='red',fontsize=14)
 maxval=max(1,float(z.mean_log1p_CP10K.max()));norm=Normalize(0,maxval)
 for i,label in enumerate(rows):
  for j,g in enumerate(cols):
   v=z[(z.author_label==label)&(z.gene==g)].iloc[0]
   if not v.matrix_gene_available:ax.scatter(j,i,marker='x',s=12,c='#999',lw=.7);continue
   if np.isfinite(v.fraction_positive):ax.scatter(j,i,s=180*v.fraction_positive,c=[v.mean_log1p_CP10K],cmap='viridis',norm=norm,edgecolors='none')
 ax.set_xticks(range(len(cols)));ax.set_xticklabels([g+'\n'+expected[g].replace('_',' ') for g in cols],rotation=65,ha='right',fontsize=7)
 ax.set_yticks(range(len(rows)));ax.set_yticklabels([f'{a} (n={int(z.loc[z.author_label.eq(a),"n_cells"].iloc[0]):,})' for a in rows],fontsize=8);ax.invert_yaxis();ax.set_xlim(-.7,len(cols)-.3);ax.set_title('Author-transferred labels: native marker review ('+name.replace('_',' ')+')',pad=12)
 fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),ax=ax,fraction=.025,pad=.025,label='Mean log1p CP10K (available genes)')
 handles=[ax.scatter([],[],s=180*f,c='#777',label=f'{int(f*100)}%') for f in [.1,.5,.9]];ax.legend(handles=handles,title='Cells with >=1 UMI',loc='upper left',bbox_to_anchor=(1.12,1),frameon=False,fontsize=8)
 ax.grid(alpha=.12);ax.set_axisbelow(True);fig.text(.015,.015,'X: input column unavailable; assay status is recorded separately; not zero expression. Rows are author mapping labels, not newly frozen subtypes.\nMarkers used for identity gating are diagnostic support, not independent validation. IGHD/XCL1 are absent from this panel.',fontsize=7.5)
 fig.subplots_adjust(left=.25,right=.82,bottom=.34,top=.88)
 for ext in ['png','svg']:
  p=out/(name+'_marker_identity_review.'+ext);fig.savefig(p,dpi=180)
  if ext=='svg':p.write_text('\n'.join(a.rstrip() for a in p.read_text().splitlines())+'\n')
 plt.close(fig)
# Spatial overview is a fixed random display subset; calculations use every QC cell.
c=pd.read_csv(r/'cell_identity_view_PRIVATE.tsv.gz',sep='\t');c=c[c.QC_keep].reset_index(drop=True);rng=np.random.default_rng(2026110802);sel=np.sort(rng.choice(len(c),min(len(c),100000),replace=False));d=c.iloc[sel]
colors={'B':'#3a78bc','T':'#8c4a99','NK':'#e1a336','Blood_EC':'#bb343b','Lymphatic_EC':'#ca7786','Mural':'#25866d','Fibroblast':'#8c7653','Other':'#d9d9d9'};marks={'B':'MS4A1/CD79A','T':'CD3E/CD3G','NK':'NCR1/KLRD1','Blood_EC':'PECAM1/VWF/CDH5','Lymphatic_EC':'PROX1/PDPN/LYVE1/FLT4','Mural':'RGS5/PDGFRB/MYH11','Fibroblast':'DCN/PDGFRA','Other':'other QC cells'}
fig,ax=plt.subplots(figsize=(10,7))
for typ in colors:
 v=d[d.cell_type.eq(typ)]
 if len(v):ax.scatter(v.x_um,v.y_um,s=.7 if typ=='Other' else 2,c=colors[typ],rasterized=True,label=typ+' | '+marks[typ],alpha=.6 if typ=='Other' else .85)
ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('Native x (micrometres)');ax.set_ylabel('Native y (micrometres)');ax.set_title(('SYNTHETIC SOFTWARE TEST ONLY | ' if 'SYNTHETIC' in r.name else '')+'Native QC identity view; centroid positions, not cell contact')
ax.legend(markerscale=5,bbox_to_anchor=(1.02,1),loc='upper left',frameon=False,fontsize=8);fig.subplots_adjust(right=.66,bottom=.12);fig.savefig(out/'native_identity_overview.png',dpi=180);plt.close(fig)
# Display all three predefined axes; all LR-positive QC cells are shown, no effect-based ROI choice.
raw=pd.read_csv(r/'cell_identity_view_PRIVATE.tsv.gz',sep='\t');keep=raw.QC_keep.to_numpy(bool);dat=np.load(r/'selected_native_counts_PRIVATE.npz');gx=list(dat['genes']);xx=dat['counts'][keep];assert len(xx)==len(c)
for ligand,receptor in [('CCL19','CCR7'),('CCL21','CCR7'),('CXCL12','CXCR4')]:
 fig,axes=plt.subplots(1,4,figsize=(17,4.8),sharex=True,sharey=True)
 for a in axes:a.scatter(d.x_um,d.y_um,s=.4,color='#ddd',rasterized=True);a.set_aspect('equal');a.set_xlabel('Native x (micrometres)')
 ec=c.cell_type.eq('Blood_EC').to_numpy()
 if 'ACKR1' in gx:
  ack=xx[:,gx.index('ACKR1')]>=1
  for mask,color,label in [(ec&ack,'#ba3b3b','ACKR1 detected EC'),(ec&~ack,'#555','ACKR1 undetected EC')]:axes[0].scatter(c.loc[mask,'x_um'],c.loc[mask,'y_um'],s=2,color=color,label=label,rasterized=True)
 else:axes[0].scatter(c.loc[ec,'x_um'],c.loc[ec,'y_um'],s=2,c='#777',label='Blood EC; ACKR1 input unavailable')
 axes[0].set_title('Vascular context | PECAM1/VWF/CDH5',fontsize=9);axes[0].legend(fontsize=6,markerscale=3,loc='upper right')
 for ax,gene in zip(axes[1:3],[ligand,receptor]):
  if gene not in gx:ax.set_title(gene+' | input column unavailable',fontsize=9);continue
  values=xx[:,gx.index(gene)];mask=values>=1;im=ax.scatter(c.loc[mask,'x_um'],c.loc[mask,'y_um'],s=2,c=np.log1p(values[mask]),cmap='magma',rasterized=True);fig.colorbar(im,ax=ax,fraction=.04,pad=.02,label='log1p native UMI');ax.set_title(gene+f'+ all identities (n={int(mask.sum()):,})',fontsize=9)
 if ligand in gx and receptor in gx:
  send=c.cell_type.isin(['Mural','Fibroblast','Blood_EC']).to_numpy()&(xx[:,gx.index(ligand)]>=1);receive=c.cell_type.isin(['B','T','NK']).to_numpy()&(xx[:,gx.index(receptor)]>=1)
  axes[3].scatter(c.loc[send,'x_um'],c.loc[send,'y_um'],s=4,c='#24846b',marker='^',label=ligand+'+ Mural/Fibroblast/EC',rasterized=True);axes[3].scatter(c.loc[receive,'x_um'],c.loc[receive,'y_um'],s=3,c='#3971b7',label=receptor+'+ B/T/NK',rasterized=True);axes[3].legend(fontsize=6,markerscale=3,loc='upper right')
 axes[3].set_title('Candidate senders and receivers | identity shown',fontsize=9);axes[0].set_ylabel('Native y (micrometres)');axes[0].invert_yaxis();fig.suptitle(('SYNTHETIC SOFTWARE TEST ONLY | ' if 'SYNTHETIC' in r.name else '')+ligand+' - '+receptor+' RNA localization; not proof of contact or recruitment',fontsize=12);fig.subplots_adjust(top=.82,bottom=.15,wspace=.3);fig.savefig(out/(ligand+'_'+receptor+'_native_RNA_localization.png'),dpi=160);plt.close(fig)
(r/('figure_manifest.json' if out.name=='figures' else 'figure_manifest_'+out.name+'.json')).write_text(json.dumps({'status':'GENERATED_PENDING_VISUAL_REVIEW','displayed_cells':len(d),'total_QC_cells':len(c),'all_cells_used_in_statistics':True,'files':[{'name':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(out.iterdir())]},indent=2)+'\n')
