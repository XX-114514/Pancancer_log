"""Alternative explicit bounded simple-path enumeration; no producer import."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import os,json,hashlib,re
import pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 O=Path(os.environ['PANCANCER_SIGNALING_PRIOR_OUT']);plan=json.loads(Path(os.environ['PANCANCER_SIGNALING_PLAN']).read_text());assert not (O/'independent_audit.json').exists() and sha(__file__)==plan['auditor_sha256'];rawp=Path(os.environ['PANCANCER_SIGNALING_RAW']);tfp=Path(os.environ['PANCANCER_SIGNED_TF_NET']);lrp=Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']);assert sha(rawp)==plan['raw_sha256'] and sha(tfp)==plan['TF_sha256'] and sha(lrp)==plan['LR_sha256'];raw=pd.read_csv(rawp,sep='\t');TFs=set(pd.read_csv(tfp,sep='\t').source);lr=pd.read_csv(lrp,sep='\t');edges=defaultdict(set)
 for r in raw.to_dict('records'):
  a=r['source_genesymbol'];b=r['target_genesymbol'];valid=r['type']=='post_translational' and r['consensus_direction'] and r['consensus_stimulation']!=r['consensus_inhibition'] and isinstance(r['references'],str) and bool(r['references']);valid=valid and isinstance(a,str) and isinstance(b,str) and bool(re.fullmatch('[A-Z0-9][A-Za-z0-9.-]*',a)) and bool(re.fullmatch('[A-Z0-9][A-Za-z0-9.-]*',b)) and a!=b
  if valid:edges[(a,b)].add(1 if r['consensus_stimulation'] else -1)
 graph=defaultdict(list);single={k:next(iter(v)) for k,v in edges.items() if len(v)==1}
 for (a,b),s in single.items():graph[a].append((b,s))
 observed=pd.read_csv(O/'private_unique_signed_edges.tsv',sep='\t');assert {(a,b):s for a,b,s in observed.itertuples(index=False,name=None)}==single
 mapping={}
 for g in set(b for a,b in single):
  alias='AP1' if g.upper().startswith(('JUN','FOS')) else 'NFKB' if g.upper().startswith(('REL','NFKB')) else g
  if alias in TFs:mapping[g]=alias
 receptors=sorted(lr.receptor.unique());roots={g for r in receptors for g in r.split('_')};rootresults={}
 for root in sorted(roots):
  out={};stack=[([root],1)]
  while stack:
   path,sign=stack.pop()
   if len(path)>1 and path[-1] in mapping:
    name=mapping[path[-1]];d=len(path)-1
    if name not in out or d<out[name][0]:out[name]=[d,0,0]
    if d==out[name][0]:out[name][1 if sign>0 else 2]+=1
   if len(path)-1<plan['max_edges']:
    for b,s in graph[path[-1]]:
     if b not in path:stack.append((path+[b],sign*s))
  rootresults[root]=out
 expected={}
 for receptor in receptors:
  rr=[rootresults[g] for g in receptor.split('_')];names=set().union(*(set(x) for x in rr))
  for name in names:
   hits=[v[name] for v in rr if name in v];d=min(v[0] for v in hits);positive=sum(v[1] for v in hits if v[0]==d);negative=sum(v[2] for v in hits if v[0]==d);expected[(receptor,name)]=(d,positive,negative)
 actual=pd.read_csv(O/'private_all_receptor_TF_paths.tsv',sep='\t');assert not actual.duplicated(['receptor','TF']).any();assert {(x.receptor,x.TF):(x.min_edges,x.positive_shortest_paths,x.negative_shortest_paths) for x in actual.itertuples(index=False)}==expected
 for x in actual.itertuples(index=False):
  path=x.example_path.split('>');assert path[0] in x.receptor.split('_') and mapping[path[-1]]==x.TF and len(path)-1==x.min_edges and len(path)==len(set(path));assert all((a,b) in single for a,b in zip(path,path[1:]));assert bool(x.sign_ambiguous)==bool(x.positive_shortest_paths and x.negative_shortest_paths)
 inv=json.loads((O/'inventory.json').read_text());assert len(actual)==inv['receptor_TF_records'];j=dict(status='PASS_ALL_SIGNED_RECEPTOR_TF_SHORTEST_PATH_COUNTS_AND_EXAMPLES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),inventory_sha256=sha(O/'inventory.json'),signed_edges_checked=len(single),roots_checked=len(roots),receptor_TF_records_checked=len(expected),no_expression_or_p_or_mechanism_validation=True);(O/'independent_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
if __name__=='__main__':main()
