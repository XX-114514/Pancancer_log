"""Author-style unsigned response and separately represented cytokine subunits.
This is a reconstruction with documented mapping, DE and enrichment implementation deviations.
"""
from pathlib import Path
import pandas as pd,numpy as np,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[1];o=r/'results/cytokine';figdir=r/'results/figures';a=pd.read_csv(o/'response_full_family_reviewed.tsv',sep='\t');a=a[a.branch=='author_unsigned_annotation_universe'];z=np.load(r/'private/CMT_identity_stats.npz');ix=pd.read_csv(r/'private/CMT_identity_stats_index.tsv',sep='\t');members=pd.read_csv(r/'audit/network_members.tsv',sep='\t');alias={'IFNg':['IFNG'],'TNFa':['TNF'],'IL12':['IL12A','IL12B'],'IL23':['IL23A','IL12B'],'IL27':['IL27','EBI3'],'TGFb1':['TGFB1']};edges=[]
for row in a[a.q_within_identity<.05].itertuples():
 genes=alias.get(row.cytokine_symbol_mouse,[row.cytokine_symbol_mouse.upper()]);responders=set(a[(a.cm==row.cm)&(a.cytokine==row.cytokine)&(a.q_within_identity<.05)].identity)
 for prod in members[members.cm==row.cm].subCluster:
  if prod in responders:continue
  mask=(ix.CMT==row.cm)&(ix.identity==prod);n=z['n'][mask].sum()
  if n<2:continue
  for gene in genes:
   if gene not in z['genes']:continue
   mean=float(z['sums'][mask][:,list(z['genes']).index(gene)].sum()/n)
   if mean>.1:edges.append(dict(cm=row.cm,producer=prod,cytokine=row.cytokine,subunit=gene,receiver=row.identity,mean_expression=mean,response_q_within_identity=row.q_within_identity,response_q_full_family=row.q_full_family))
e=pd.DataFrame(edges);e.to_csv(o/'author_unsigned_subunit_network.tsv',sep='\t',index=False)
v=a.pivot_table(index='cytokine',columns='cm',values='q_within_identity',aggfunc='min');fig,ax=plt.subplots(figsize=(6,10));im=ax.imshow(-np.log10(v.clip(lower=1e-20)),aspect='auto',cmap='magma',vmin=0,vmax=5);ax.set_xticks(range(len(v.columns)),v.columns,rotation=45);ax.set_yticks(range(len(v)),v.index,fontsize=8);ax.set_title('Author-style unsigned response reconstruction\nMinimum subset q is descriptive, not CM-adjusted q',fontsize=10);fig.colorbar(im,ax=ax,label='-log10(minimum subset q), capped at 5');fig.tight_layout();fig.savefig(figdir/'author_unsigned_response_minimum_q.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(3,3,figsize=(15,12))
for i,ax in enumerate(axs.flat,1):
 cm=f'CM{i:02d}';g=e[e.cm==cm].copy();ax.set_title(cm);ax.axis('off')
 if g.empty:ax.text(.5,.5,'No qualified subunit-response edge',ha='center',transform=ax.transAxes);continue
 g['cytokine_subunit']=g.cytokine+' ['+g.subunit+']';prods=sorted(g.producer.unique());rec=sorted(g.receiver.unique());cyts=sorted(g.cytokine_subunit.unique())
 for j,n in enumerate(prods):ax.text(0,(j+.5)/len(prods),n,fontsize=5,ha='right')
 for j,n in enumerate(rec):ax.text(1,(j+.5)/len(rec),n,fontsize=5,ha='left')
 for j,n in enumerate(cyts):ax.text(.5,(j+.5)/len(cyts),n,fontsize=5,ha='center',bbox=dict(facecolor='white',edgecolor='none',pad=0))
 for row in g.drop_duplicates(['producer','cytokine_subunit','receiver']).itertuples():
  y=(cyts.index(row.cytokine_subunit)+.5)/len(cyts);ax.plot([0,.5],[(prods.index(row.producer)+.5)/len(prods),y],lw=.4,alpha=.25,c='purple');ax.plot([.5,1],[y,(rec.index(row.receiver)+.5)/len(rec)],lw=.4,alpha=.25,c='teal')
fig.suptitle('Author-style unsigned response / individual cytokine subunits / not causal');fig.tight_layout();fig.savefig(figdir/'author_unsigned_subunit_network.png',dpi=160,bbox_inches='tight');plt.close(fig)
receipt={'branch':'author_unsigned_annotation_universe','q_rule':'within_CM_identity_q_lt_0.05; implementation deviations documented','n_response_associations_within_q':int((a.q_within_identity<.05).sum()),'n_response_associations_full_family_q':int((a.q_full_family<.05).sum()),'n_subunit_network_edges':len(e),'separate_subunits_not_complete_cytokine_proof':True};(r/'audit/author_response_network_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
