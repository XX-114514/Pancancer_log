"""Prespecified all-CM decomposition of candidate-count changes, not differential-LR tests."""
from pathlib import Path
import pandas as pd,numpy as np
r=Path(__file__).resolve().parents[1];o=r/'results/summary';pairs=pd.read_csv(o/'all_CM_CCI_pairs.tsv',sep='\t');out=[];details=[]
for cm in sorted(pairs.cm.unique()):
 p=pairs[(pairs.cm==cm)&(pairs.branch=='network')&(~pairs.self_pair)].dropna(subset=['High','Low_pooled']);delta=p.High-p.Low_pooled;denom=float(abs(delta).sum());lr=pd.read_csv(o/(cm+'_LR_detection.tsv.gz'),sep='\t');lr=lr[(lr.source!=lr.target)&lr.cellpair_present_both];lr=lr[lr.detection.isin(['both','High_only','Low_only','neither'])];lr['candidate_delta']=lr.High_significant_mean.notna().astype(int)-lr.Low_significant_mean.notna().astype(int);s=lr.groupby('id_cp_interaction').candidate_delta.sum();ldenom=float(abs(s).sum());out.append({'cm':cm,'nonself_pairs':len(p),'net_candidate_change':float(delta.sum()),'sum_absolute_pair_changes':denom,'top3_pair_fraction_absolute_change':float(abs(delta).nlargest(3).sum()/denom) if denom else np.nan,'top10_LR_fraction_absolute_net_change':float(abs(s).nlargest(10).sum()/ldenom) if ldenom else np.nan,'interpretation':'descriptive_count_decomposition; not_patient_or_LR_difference_p'})
 for key,v in s.items():details.append(dict(cm=cm,id_cp_interaction=key,net_candidate_event_change=int(v)))
pd.DataFrame(out).to_csv(o/'HighLow_concentration.tsv',sep='\t',index=False);pd.DataFrame(details).to_csv(o/'HighLow_LR_contributions.tsv',sep='\t',index=False)
