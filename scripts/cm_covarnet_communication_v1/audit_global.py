from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
r=Path(__file__).resolve().parents[1];d=r/'results/global';receipt=json.loads((d/'complete.json').read_text());assert receipt['iterations']==1000
for name,sha in receipt['outputs'].items():assert hashlib.sha256((d/name).read_bytes()).hexdigest()==sha
p=pd.read_csv(next(d.glob('*_pvalues_*.txt')),sep='\t').set_index('id_cp_interaction');m=pd.read_csv(next(d.glob('*_means_*.txt')),sep='\t').set_index('id_cp_interaction').reindex(p.index);s=pd.read_csv(next(d.glob('*_significant_means_*.txt')),sep='\t').set_index('id_cp_interaction').reindex(p.index);cols=[c for c in p if '|' in c];expected=(p[cols]<=.05);assert (s[cols].notna()==expected).all().all();n=s[cols].notna().sum();pairs=pd.DataFrame([{'source':k.split('|')[0],'target':k.split('|')[1],'CCI':int(v)} for k,v in n.items()]);members=pd.read_csv(r/'audit/network_members.tsv',sep='\t');rows=[]
for cm,g in members.groupby('cm'):
 a=pairs.source.isin(g.subCluster);b=pairs.target.isin(g.subCluster)
 for name,mask in [('inside',a&b),('sender_inside_receiver_outside',a&~b),('sender_outside_receiver_inside',~a&b),('both_outside',~a&~b)]:
  x=pairs[mask&(pairs.source!=pairs.target)];rows.append(dict(cm=cm,background=name,n_nonself_pairs=len(x),mean_CCI=x.CCI.mean(),median_CCI=x.CCI.median(),sum_CCI=x.CCI.sum(),zero_CCI_pairs=int((x.CCI==0).sum())))
pd.DataFrame(rows).to_csv(r/'audit/global_inside_vs_outside.tsv',sep='\t',index=False)
res={'status':'technical_and_descriptive_readback_pass','database_measurable_rows':len(p),'ordered_pairs':len(cols),'candidate_events':int(n.sum()),'nonzero_pairs':int((n>0).sum()),'max_CCI':int(n.max()),'permutation_p_not_environment_p':True};(r/'audit/global_acceptance.json').write_text(json.dumps(res,indent=2));print(res)
