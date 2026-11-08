import pathlib,sys,json
import pandas as pd,numpy as np
h,r=map(pathlib.Path,sys.argv[1:3]);records=[]
for pool in ['T_NK','B_Plasma','Stromal','Myeloid','Endothelial']:
 s=pd.read_csv(h/f'p0/consensus/{pool}.seeds.tsv.gz',sep='\t',dtype=str).set_index('cell_id')
 d=pd.read_csv(r/f'consensus_v1/{pool}.cells.tsv.gz',sep='\t').set_index('cell_id');assert d.index.equals(s.index)
 names=sorted(d.review_label.unique());trans=pd.DataFrame(0.,index=names,columns=names);pred=[]
 for col in [f'seed_{n}' for n in range(2026110800,2026110810)]:
  ct=pd.crosstab(s[col],d.review_label).reindex(columns=names,fill_value=0);winner=ct.idxmax(axis=1);p=s[col].map(winner);pred.append(p.values)
  trans+=pd.crosstab(d.review_label,p).reindex(index=names,columns=names,fill_value=0)
 rates=trans.div(d.review_label.value_counts()*10,axis=0)
 pairs=pd.read_csv(r/f'consensus_v1/{pool}.pairs.tsv',sep='\t')
 for row in pairs.itertuples():
  a=rates.loc[row.label_a,row.label_b];b=rates.loc[row.label_b,row.label_a]
  assert np.isclose(a,row.a_to_b,atol=1e-14) and np.isclose(b,row.b_to_a,atol=1e-14)
  assert ((a+b)/2>=.25)==row.merge
 components=json.loads((r/f'consensus_v1/{pool}.components.json').read_text());lookup={x:k for k,v in components.items() for x in v}
 P=np.stack([pd.Series(p).map(lookup).values for p in pred],axis=1)
 # Check reported winner directly against independently mapped votes, and maximum count.
 reported=(P==d.consensus_label.values[:,None]).sum(1);maxcount=np.zeros(len(d),dtype=int)
 for comp in components:maxcount=np.maximum(maxcount,(P==comp).sum(1))
 assert np.array_equal(reported,maxcount);assert np.allclose(reported/10,d.stability)
 assert not (d.include_main_subtype&((d.stability<.7)|d.low_depth|~d.retained_pool|d.tied_mode)).any()
 records.append({'pool':pool,'rows':len(d),'independent_crosstab_replay':'passed'});print(json.dumps(records[-1]),flush=True)
(r/'consensus_replay_acceptance.json').write_text(json.dumps({'status':'passed','checks':records},indent=2))
