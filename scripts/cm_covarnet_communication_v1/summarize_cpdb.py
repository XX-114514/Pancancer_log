from pathlib import Path
import pandas as pd,numpy as np,json,itertools
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
r=Path(__file__).resolve().parents[1];o=r/'results/summary';o.mkdir(exist_ok=True);mem=pd.read_csv(r/'audit/member_comparison.tsv',sep='\t');cache={}
def read(tag):
 if tag in cache:return cache[tag]
 d=r/'results'/tag
 if not (d/'complete.json').exists():return None
 a=pd.read_csv(next(d.glob('statistical_analysis_significant_means_*.txt')),sep='\t');cols=[c for c in a if '|' in c];counts=a[cols].notna().sum();cache[tag]=(counts,a);return cache[tag]
globalres=read('global')
if globalres:
 counts,a=globalres;rows=[]
 for c,n in counts.items():
  x,y=c.split('|');shared=[cm for cm,g in mem[mem.network_member].groupby('cm') if x in set(g.subCluster) and y in set(g.subCluster)];rows.append(dict(source=x,target=y,candidate_CCI=int(n),shared_CM=';'.join(shared),within_any=bool(shared)))
 g=pd.DataFrame(rows);g.to_csv(o/'global_cellpair_CCI.tsv',sep='\t',index=False)
 cmrows=[]
 for c1,m1 in mem[mem.network_member].groupby('cm'):
  for c2,m2 in mem[mem.network_member].groupby('cm'):
   x=g[g.source.isin(m1.subCluster)&g.target.isin(m2.subCluster)];cmrows.append(dict(source_CM=c1,target_CM=c2,n_pairs=len(x),mean_CCI=x.candidate_CCI.mean(),sum_CCI=x.candidate_CCI.sum(),counting='inclusive_nonadditive'))
 pd.DataFrame(cmrows).to_csv(o/'global_CM_CCI.tsv',sep='\t',index=False)
allrows=[]
for branch in ['network','old_top10']:
 for cm,t in mem.groupby('cm'):
  nodes=t.loc[t.network_member if branch=='network' else t.top10_rank.notna(),'subCluster'].tolist();hi=read(cm+'_High');lo=read(cm+'_Low');others=[read(x+'_High') for x in sorted(mem.cm.unique()) if x!=cm]
  for s,tg in itertools.product(nodes,repeat=2):
   pair=s+'|'+tg;h=hi[0].get(pair,np.nan) if hi else np.nan;l=lo[0].get(pair,np.nan) if lo else np.nan;vals=[x[0].get(pair,np.nan) if x else np.nan for x in others];med=float(np.median(vals)) if np.isfinite(vals).all() else np.nan
   allrows.append(dict(branch=branch,cm=cm,source=s,target=tg,self_pair=s==tg,High=h,Low_pooled=l,Low_median_CMT=med,delta_pooled=h-l,delta_median_CMT=h-med,status='evaluated' if np.isfinite(h+l) else 'not_evaluated'))
  # complete LR comparison for network pairs; no rank gate, preserve missing cells and database IDs.
  if branch=='network' and hi and lo:
   ah=hi[1].set_index('id_cp_interaction');al=lo[1].set_index('id_cp_interaction');lr=[]
   mh=pd.read_csv(next((r/'results'/(cm+'_High')).glob('statistical_analysis_means_*.txt')),sep='\t').set_index('id_cp_interaction');ml=pd.read_csv(next((r/'results'/(cm+'_Low')).glob('statistical_analysis_means_*.txt')),sep='\t').set_index('id_cp_interaction')
   for s,tg in itertools.product(nodes,repeat=2):
    col=s+'|'+tg;ids=ah.index.union(al.index);hv=ah[col].reindex(ids) if col in ah else pd.Series(np.nan,index=ids);lv=al[col].reindex(ids) if col in al else pd.Series(np.nan,index=ids)
    status=np.select([hv.notna()&lv.notna(),hv.notna(),lv.notna()],['both','High_only','Low_only'],default='neither');status=status.astype(object);status[:]='not_evaluable_missing_cellpair' if col not in ah or col not in al else status;status[~ids.isin(ah.index.intersection(al.index))]='not_evaluable_LR_not_retained_both_groups';lr.append(pd.DataFrame({'cm':cm,'source':s,'target':tg,'id_cp_interaction':ids,'High_significant_mean':hv.values,'Low_significant_mean':lv.values,'High_mean':mh[col].reindex(ids).values if col in mh else np.nan,'Low_mean':ml[col].reindex(ids).values if col in ml else np.nan,'detection':status,'cellpair_present_both':col in ah and col in al}))
   pd.concat(lr).to_csv(o/f'{cm}_LR_detection.tsv.gz',sep='\t',index=False)
pairs=pd.DataFrame(allrows);pairs.to_csv(o/'all_CM_CCI_pairs.tsv',sep='\t',index=False);stats=[]
for (branch,cm),g in pairs[~pairs.self_pair].groupby(['branch','cm']):
 for alt in ['pooled','median_CMT']:
  x=g[['High','Low_'+alt]].dropna();d=x.High-x['Low_'+alt];p=float(wilcoxon(x.High,x['Low_'+alt],alternative='two-sided',zero_method='wilcox',method='approx',correction=True).pvalue) if len(x) and (d!=0).any() else (1. if len(x) else np.nan)
  stats.append(dict(branch=branch,cm=cm,comparison=alt,n_cellpairs=len(x),n_patients='not_applicable',median_delta=d.median(),mean_delta=d.mean(),p_cellpair_Wilcoxon=p,status='evaluated_descriptive_pairs_dependent' if len(x) else 'pending_or_unmeasurable'))
st=pd.DataFrame(stats);st['q_extension_all_CM']=np.nan
for _,g in st.groupby(['branch','comparison']):st.loc[g.index,'q_extension_all_CM']=multipletests(g.p_cellpair_Wilcoxon.fillna(1),method='fdr_bh')[1]
st.loc[st.p_cellpair_Wilcoxon.isna(),'q_extension_all_CM']=np.nan
st.to_csv(o/'all_CM_CCI_tests.tsv',sep='\t',index=False);print(st.to_string(index=False))

# Decompose only after the complete prespecified High/Low family is available.
if all((r/'results'/f'CM{i:02d}_{g}'/'complete.json').exists() for i in range(1,10) for g in ['High','Low']):
 import subprocess,sys
 subprocess.run([sys.executable,str(r/'scripts/summarize_environment_drivers.py')],check=True)
