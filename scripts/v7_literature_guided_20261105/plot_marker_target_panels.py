"""Scientific dot plots with explicit marker-to-type/state panels and row labels."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
P=argparse.ArgumentParser();P.add_argument('--marker-table',type=Path,required=True);P.add_argument('--row-labels',type=Path,required=True);P.add_argument('--panels',type=Path,required=True);P.add_argument('--output-dir',type=Path,required=True);args=P.parse_args()
args.output_dir.mkdir(parents=True,exist_ok=False)
e=pd.read_csv(args.marker_table,sep='\t');rows=pd.read_csv(args.row_labels,sep='\t');spec=json.loads(args.panels.read_text());labels=rows.set_index('cluster').display_annotation.to_dict();order=sorted(e.cluster.unique())
assert set(order)==set(labels) and not e.duplicated(['cluster','gene']).any()
assert ((e.n_measured>=0)&(e.n_measured<=e.n_cells)).all();assert e.groupby('cluster').n_cells.nunique().max()==1
ncell=e.groupby('cluster').n_cells.first().to_dict();plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','font.size':10})
short={'mast_identity':'Mast\nidentity','MCTC':'MC_TC-like\n(candidate)','MCT':'MC_T-like\n(candidate)','IL1RL1':'IL1RL1/lipid\nmediator axis','cycling':'Cycling Mast\n(state)','IFN':'IFN response\n(state)','HLAII':'HLA-II antigen\npresentation','T_competitor':'T cells\n(competitor)','myeloid_competitor':'Myeloid cells\n(competitor)','stromal_competitor':'Fibroblasts\n(competitor)','epithelial_competitor':'Epithelial cells\n(competitor)'}
palette={'parent':'#bdd7e7','candidate_subtype':'#cbc9e2','state_axis':'#fdd49e','state':'#fdd49e','state_program':'#fdd49e','competing_parent':'#fcbba1','published_reference_state':'#c7e9c0'}
def plot(groups,name,title):
 positions=[];bands=[];x=0
 for p in groups:
  lo=x
  for gene in p['genes']:positions.append((x,gene));x+=1
  bands.append((lo,x-1,p));x+=.8
 assert set(g for _,g in positions)<=set(e.gene)
 fig=plt.figure(figsize=(max(19,len(positions)*.44+8),13));ax=fig.add_axes([.30,.15,.57,.66]);stats={}
 for g in set(g for _,g in positions):
  v=np.log1p(e.loc[e.gene.eq(g),'mean_cp10k_measured'].to_numpy());v=v[np.isfinite(v)];stats[g]=(float(v.mean()),float(v.std())) if len(v) else (0,1)
 xs=[];ys=[];size=[];color=[];missing=[]
 for xx,g in positions:
  mu,sd=stats[g]
  for row in e[e.gene.eq(g)].to_dict('records'):
   y=order.index(row['cluster'])
   if row['coverage_fraction']<.8 or not np.isfinite(row['pct_expressing_measured']):missing.append((xx,y));continue
   xs.append(xx);ys.append(y);size.append(125*row['pct_expressing_measured']);color.append(np.clip((np.log1p(row['mean_cp10k_measured'])-mu)/(sd or 1),-2,2))
 sc=ax.scatter(xs,ys,s=size,c=color,cmap='RdBu_r',vmin=-2,vmax=2,edgecolors='none')
 if missing:xx,yy=zip(*missing);ax.scatter(xx,yy,s=16,c='#969696',marker='x')
 ax.set_xticks([v[0] for v in positions],[v[1] for v in positions],rotation=65,ha='right',fontsize=10)
 ax.set_yticks(range(len(order)),[f'C{k:02d} | {labels[k]}  (n={ncell[k]:,})' for k in order],fontsize=10)
 ax.set_ylim(len(order)-.7,-.7);ax.set_xlim(-.65,x-.9);ax.grid(axis='y',alpha=.14);ax.set_axisbelow(True)
 for text,k in zip(ax.get_yticklabels(),order):
  if '[review]' in labels[k]:text.set_color('#a63630')
 top=ax.secondary_xaxis('top');top.set_xticks([(lo+hi)/2 for lo,hi,p in bands],[short.get(p['id'],p['target']) for lo,hi,p in bands],fontsize=9);top.tick_params(length=0,pad=19)
 for lo,hi,p in bands:
  ax.add_patch(Rectangle((lo-.42,1.015),hi-lo+.84,.017,transform=ax.get_xaxis_transform(),facecolor=palette[p.get('level','published_reference_state')],clip_on=False,edgecolor='none'))
  ax.plot([lo-.42,hi+.42],[1.007,1.007],transform=ax.get_xaxis_transform(),color='#777777',lw=.7,clip_on=False)
 cbax=fig.add_axes([.90,.26,.012,.38]);cb=fig.colorbar(sc,cax=cbax);cb.set_label('Gene-scaled log1p(mean CP10K)',fontsize=10)
 handles=[ax.scatter([],[],s=125*f,color='#777777',label=f'{int(f*100)}% detected') for f in [.1,.5,1]];ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.025,1.01),frameon=False,fontsize=10)
 fig.text(.30,.965,title,fontsize=19,weight='bold');fig.text(.30,.925,'Top: corresponding identity / candidate subtype / state. Left: working annotation and cell count.',fontsize=11)
 fig.text(.30,.035,'Size: detection among source-measured cells; x: source coverage <80%. Color: relative expression within each gene.\nRows are exploratory partitions, not accepted fine subtypes. MC_T/MC_TC are RNA-based candidates.\nReference C0-C9 indicate published pan-MC panels, not the V7 row IDs. Original V7 annotation is retained.',fontsize=10)
 for suffix in ['png','svg']:fig.savefig(args.output_dir/(name+'.'+suffix),dpi=180)
 plt.close(fig)
main=[p for p in spec['panels'] if p['level']!='competing_parent'];qc=[p for p in spec['panels'] if p['level'] in ['parent','competing_parent']]
plot(main,'Mast_identity_subtype_state_v3','V7 Mast | identity, candidate subtype and state markers')
plot(qc,'Mast_competing_lineage_QC_v3','V7 Mast | competing lineage markers and QC review')
plot([spec['panels'][0]]+spec['panMC_reference_panels']+[p for p in spec['panels'] if p['id']=='T_competitor'],'Mast_panMC_reference_states_v3','V7 Mast | published pan-MC state marker panels')
manifest={'n_clusters':len(order),'n_cells':int(sum(ncell.values())),'n_marker_records':len(e),'row_labels_explicit':True,'marker_target_labels_explicit':True,'candidate_status_explicit':True,'input_marker_sha256':hashlib.sha256(args.marker_table.read_bytes()).hexdigest(),'panel_spec_sha256':hashlib.sha256(args.panels.read_bytes()).hexdigest(),'outputs':[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in args.output_dir.iterdir()]}
(args.output_dir/'plot_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({'plots':3,'cells':manifest['n_cells'],'output_dir':str(args.output_dir)}))
