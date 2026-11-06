from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[1];o=r/'results/figures';o.mkdir(exist_ok=True);s=r/'results/summary'
# Empty panels are explicitly marked unestimable, never plotted as zero.
if (s/'global_CM_CCI.tsv').exists():
 a=pd.read_csv(s/'global_CM_CCI.tsv',sep='\t').pivot(index='source_CM',columns='target_CM',values='mean_CCI');fig,ax=plt.subplots(figsize=(7,6));im=ax.imshow(a, cmap='viridis');ax.set_xticks(range(9),a.columns,rotation=90);ax.set_yticks(range(9),a.index);ax.set_title('20260825 V5 K9: global mean candidate CCI\nInclusive CM membership; non-additive');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(o/'global_CM_CCI.png',dpi=180);plt.close(fig)
if (s/'all_CM_CCI_pairs.tsv').exists():
 a=pd.read_csv(s/'all_CM_CCI_pairs.tsv',sep='\t');a=a[(a.branch=='network')&(~a.self_pair)]
 for comp in ['pooled','median_CMT']:
  fig,axs=plt.subplots(3,3,figsize=(11,10))
  for i,ax in enumerate(axs.flat,1):
   cm=f'CM{i:02d}';g=a[a.cm==cm][['High','Low_'+comp]].dropna()
   if len(g):
    for _,v in g.iterrows():ax.plot([0,1],[v.iloc[1],v.iloc[0]],color='grey',alpha=.25,lw=.5)
    ax.scatter(np.zeros(len(g)),g.iloc[:,1],s=5);ax.scatter(np.ones(len(g)),g.High,s=5);ax.plot([0,1],[g.iloc[:,1].mean(),g.High.mean()],c='black',lw=2)
   else:ax.text(.5,.5,'Pending / not evaluable',ha='center',transform=ax.transAxes)
   ax.set_title(cm+f' (pairs={len(g)})');ax.set_xticks([0,1],['Low','High']);ax.set_ylabel('Candidate CCI count')
  fig.suptitle('20260825 V5 K9 / '+comp+' / ordered non-self pairs (not patients)');fig.tight_layout();fig.savefig(o/f'HighLow_CCI_{comp}.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,4));ax.axis('off');ax.text(.5,.6,'Patient-level effects: not estimable under frozen eligibility',ha='center',fontsize=13);ax.text(.5,.4,'No source × cancer × material × known-timepoint stratum\nwith at least 10 documented patients in BOTH High and Low.\nUnknown patients/timepoints are not treated as independent untreated cases.',ha='center');fig.savefig(o/'patient_effect_eligibility.png',dpi=180,bbox_inches='tight');plt.close(fig)
c=r/'results/cytokine'
if (c/'response_full_family.tsv').exists():
 a=pd.read_csv(c/'response_full_family_reviewed.tsv',sep='\t');a=a[a.branch=='project_directional_measured_universe'];a['label']=a.cm+' / '+a.identity+' / '+a.direction
 # All tested cytokines shown, not significant-only selection.
 v=a.pivot_table(index='label',columns='cytokine',values='q_full_family',aggfunc='min');fig,ax=plt.subplots(figsize=(18,max(8,len(v)*.16)));im=ax.imshow(-np.log10(v.clip(lower=1e-20)),aspect='auto',cmap='magma',vmin=0,vmax=5);ax.set_yticks(range(len(v)),v.index,fontsize=5);ax.set_xticks(range(len(v.columns)),v.columns,rotation=90,fontsize=6);ax.set_title('Response inference: full-family q; cell-level, not patient or causal evidence');fig.colorbar(im,ax=ax,label='-log10(q) capped at 5');fig.tight_layout();fig.savefig(o/'cytokine_response_all_reviewed.png',dpi=180);plt.close(fig)
 e=pd.read_csv(c/'production_response_network_subunit_review.tsv',sep='\t');e=e[e.status=='all_required_subunits_expressed'];fig,axs=plt.subplots(3,3,figsize=(15,12))
 for i,ax in enumerate(axs.flat,1):
  cm=f'CM{i:02d}';g=e[e.cm==cm];prods=sorted(g.producer.unique());rec=sorted(g.receiver.unique());ax.set_title(cm);ax.axis('off')
  if g.empty:ax.text(.5,.5,'No qualified producer-response edge',ha='center',transform=ax.transAxes);continue
  for j,n in enumerate(prods):ax.text(0,(j+.5)/len(prods),n,fontsize=5,ha='right')
  for j,n in enumerate(rec):ax.text(1,(j+.5)/len(rec),n,fontsize=5,ha='left')
  cyts=sorted(g.cytokine.unique())
  for j,n in enumerate(cyts):ax.text(.5,(j+.5)/len(cyts),n,fontsize=5,ha='center',bbox=dict(facecolor='white',edgecolor='none',pad=0))
  for row in g.drop_duplicates(['producer','cytokine','receiver']).itertuples():
   yy=(cyts.index(row.cytokine)+.5)/len(cyts);ax.plot([0,.5],[(prods.index(row.producer)+.5)/len(prods),yy],lw=.4,alpha=.2,color='purple');ax.plot([.5,1],[yy,(rec.index(row.receiver)+.5)/len(rec)],lw=.4,alpha=.2,color='teal')
 fig.suptitle('Production–response inference (all qualified edges; not causal LR communication)');fig.tight_layout();fig.savefig(o/'production_response_network_qualified.png',dpi=180,bbox_inches='tight');plt.close(fig)
print('figures generated')
