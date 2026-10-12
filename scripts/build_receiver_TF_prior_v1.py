"""External only: signed shortest receiver-component -> TF paths <=3.
No expression, effects, p-values, or biological validation in this stage.
"""
from pathlib import Path
from collections import defaultdict,deque
from datetime import datetime,timezone
import json,os,hashlib,re
import pandas as pd

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def TFname(g):
 u=g.upper()
 return 'AP1' if u.startswith(('JUN','FOS')) else 'NFKB' if u.startswith(('REL','NFKB')) else g

def main():
 planp=Path(os.environ['PANCANCER_SIGNALING_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'];rawp=Path(os.environ['PANCANCER_SIGNALING_RAW']);tfp=Path(os.environ['PANCANCER_SIGNED_TF_NET']);lrp=Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']);O=Path(os.environ['PANCANCER_SIGNALING_PRIOR_OUT']);assert not O.exists();assert sha(rawp)==plan['raw_sha256'] and sha(tfp)==plan['TF_sha256'] and sha(lrp)==plan['LR_sha256'];start=datetime.now(timezone.utc).isoformat();raw=pd.read_csv(rawp,sep='\t');tf=pd.read_csv(tfp,sep='\t');lr=pd.read_csv(lrp,sep='\t');TFs=set(tf.source);gene=re.compile(r'^[A-Z0-9][A-Za-z0-9.-]*$');signs=defaultdict(set);rows=0
 for x in raw.itertuples(index=False):
  if not (x.type=='post_translational' and x.consensus_direction and (x.consensus_stimulation != x.consensus_inhibition) and isinstance(x.references,str) and x.references):continue
  if not (isinstance(x.source_genesymbol,str) and isinstance(x.target_genesymbol,str) and gene.fullmatch(x.source_genesymbol) and gene.fullmatch(x.target_genesymbol)):continue
  if x.source_genesymbol==x.target_genesymbol:continue
  rows+=1;signs[(x.source_genesymbol,x.target_genesymbol)].add(1 if x.consensus_stimulation else -1)
 graph=defaultdict(list);edges=[]
 for (a,b),v in sorted(signs.items()):
  if len(v)==1:
   sign=next(iter(v));graph[a].append((b,sign));edges.append(dict(source=a,target=b,sign=sign))
 receptors=sorted(lr.receptor.unique());roots=sorted({g for r in receptors for g in r.split('_')});reach={}
 for root in roots:
  dist={root:0};counts={(root,1):1};examples={(root,1):[root]};queue=deque([root])
  while queue:
   a=queue.popleft();d=dist[a]
   if d>=plan['max_edges']:continue
   for b,sg in graph[a]:
    if b not in dist:dist[b]=d+1;queue.append(b)
    if dist[b]!=d+1:continue
    for sign in [-1,1]:
     count=counts.get((a,sign),0)
     if count:
      key=(b,sign*sg);counts[key]=counts.get(key,0)+count;examples.setdefault(key,examples[(a,sign)]+[b])
  terms=defaultdict(list)
  for b,d in dist.items():
   if d>0 and TFname(b) in TFs:terms[TFname(b)].append((b,d))
  out={}
  for name,ends in terms.items():
   depth=min(d for _,d in ends);sel=[b for b,d in ends if d==depth];pos=sum(counts.get((b,1),0) for b in sel);neg=sum(counts.get((b,-1),0) for b in sel);example=next(examples[(b,sg)] for b in sel for sg in [1,-1] if (b,sg) in examples);out[name]=dict(depth=depth,positive_shortest_paths=pos,negative_shortest_paths=neg,example_path='>'.join(example))
  reach[root]=out
 result=[]
 for receptor in receptors:
  components=receptor.split('_');names={name for root in components for name in reach[root]}
  for name in sorted(names):
   hits=[(root,reach[root][name]) for root in components if name in reach[root]];depth=min(v['depth'] for _,v in hits);hits=[(root,v) for root,v in hits if v['depth']==depth];pos=sum(v['positive_shortest_paths'] for _,v in hits);neg=sum(v['negative_shortest_paths'] for _,v in hits);result.append(dict(receptor=receptor,TF=name,min_edges=depth,signal_components=';'.join(root for root,_ in hits),positive_shortest_paths=pos,negative_shortest_paths=neg,sign_ambiguous=bool(pos and neg),expected_sign=0 if pos and neg else 1 if pos else -1,example_path=hits[0][1]['example_path'],candidate_compatibility_only=True))
 O.mkdir();pd.DataFrame(edges).to_csv(O/'private_unique_signed_edges.tsv',sep='\t',index=False);d=pd.DataFrame(result);d.to_csv(O/'private_all_receptor_TF_paths.tsv',sep='\t',index=False);summary=d.groupby('receptor').agg(TFs_with_short_path=('TF','size'),TFs_unambiguous_sign=('sign_ambiguous',lambda x:int((~x).sum())),min_edges=('min_edges','min')).reindex(receptors).fillna(0).reset_index();summary.to_csv(O/'receptor_prior_coverage_summary.tsv',sep='\t',index=False)
 j=dict(status='COMPLETE_EXTERNAL_PRIOR_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),eligible_raw_rows=rows,unique_signed_edges=len(edges),ambiguous_gene_edges_excluded=sum(len(v)>1 for v in signs.values()),resource_receptors=len(receptors),receptors_with_short_TF_paths=int((summary.TFs_with_short_path>0).sum()),receptor_TF_records=len(d),unambiguous_shortest_sign_records=int((~d.sign_ambiguous).sum()),resource_LR_pairs_with_receptor_paths=int(lr.receptor.isin(d.receptor).sum()),no_expression_or_p=True,no_binding_or_activation_validation=True,outputs={f.name:sha(f) for f in O.iterdir() if f.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
if __name__=='__main__':main()
