import os
from pathlib import Path
import csv,json,shutil
r=Path(__file__).resolve().parents[1]; src=(Path(os.environ['CODE_ROOT']) / 'releases/annotation_v7_20260905'); old=(Path(os.environ['PC_ROOT']) / 'Project_v5/frozen/ontology')
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(name,rows):
 with (r/'tables'/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={}
for name in ['major_v5_to_v7_transition.tsv','malignancy_v3_to_v7_transition.tsv']:
 rows=read(src/'tables'/name); first=list(rows[0])[0]; cols=list(rows[0])[1:]; total=sum(int(x[c]) for x in rows for c in cols); assert total==4676787
 long=[]; marg=[]
 for x in rows:
  n=sum(int(x[c]) for c in cols)
  for c in cols:
   if int(x[c]):long.append(dict(v5=x[first],v7=c,n_cells=int(x[c]),percent_of_v5=100*int(x[c])/n))
 labels=sorted(set(cols)|{x[first] for x in rows})
 for label in labels:
  before=sum(int(x[c]) for x in rows if x[first]==label for c in cols);after=sum(int(x.get(label,0)) for x in rows)
  retained=sum(int(x.get(label,0)) for x in rows if x[first]==label)
  marg.append(dict(label=label,V5=before,V7=after,delta=after-before,retained=retained,inflow=after-retained,outflow=before-retained))
 write(name.replace('.tsv','_long.tsv'),sorted(long,key=lambda x:-x['n_cells']));write(name.replace('_transition','_marginals'),marg)
 same=sum(int(x.get(x[first],0)) for x in rows);summary[name]={'total':total,'exact_label_retained':same,'exact_label_changed':total-same,'changed_pct':100*(total-same)/total}
 shutil.copyfile(src/'tables'/name,r/'tables'/name)
for a,b,name in [(old/'identity_annotation_final_counts.tsv',src/'tables/final_annotation_v7_counts.tsv','fine_identity_count_comparison.tsv'),(old/'malignancy_tier_v3_counts.tsv',src/'tables/malignancy_tier_v7_counts.tsv','malignancy_tier_count_comparison.tsv')]:
 def counts(p):
  rows=read(p); k=list(rows[0])[0];return {x[k]:int(x['n_cells']) for x in rows}
 x,y=counts(a),counts(b);assert sum(x.values())==sum(y.values())==4676787
 write(name,[dict(label=k,V5=x.get(k,0),V7=y.get(k,0),delta=y.get(k,0)-x.get(k,0)) for k in sorted(x.keys()|y.keys())])
(r/'tables/summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
for row in read(r/'tables/major_v5_to_v7_marginals.tsv'):
 if row['delta']!='0':print(row['label'],row['V5'],row['V7'],row['delta'])
