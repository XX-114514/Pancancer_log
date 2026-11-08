"""Literal documented 10-seed majority mapping and bidirectional 0.25 union rule.
No seed re-fit. All reviewed labels participate in diagnostics; QC/state/low-depth
components are explicitly ineligible for main identity proportions.
"""
import sys,pathlib,json
import numpy as np,pandas as pd

def consensus(original,seeds,merge_threshold=.25):
 names=np.array(sorted(set(original)),dtype=str); lookup={x:i for i,x in enumerate(names)}
 y=np.array([lookup[x] for x in original],dtype=np.int32); k=len(names); votes=np.empty(seeds.shape,dtype=np.int32)
 mapping=[]
 for s in range(seeds.shape[1]):
  cl,codes=np.unique(seeds[:,s],return_inverse=True); tab=np.bincount(codes*k+y,minlength=len(cl)*k).reshape(len(cl),k)
  win=tab.argmax(1);votes[:,s]=win[codes]
  for j,c in enumerate(cl):mapping.append(dict(seed=s,cluster=c,label=names[win[j]],majority_fraction=float(tab[j,win[j]]/tab[j].sum()),tied_majority=bool((tab[j]==tab[j].max()).sum()>1)))
 trans=np.zeros((k,k),dtype=np.int64)
 for s in range(seeds.shape[1]):trans+=np.bincount(y*k+votes[:,s],minlength=k*k).reshape(k,k)
 den=np.bincount(y,minlength=k)*seeds.shape[1];rate=trans/den[:,None]
 parent=list(range(k))
 def find(a):
  while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
  return a
 pairs=[]
 for a in range(k):
  for b in range(a+1,k):
   ab=float((rate[a,b]+rate[b,a])/2);merge=ab>=merge_threshold
   pairs.append(dict(label_a=names[a],label_b=names[b],a_to_b=float(rate[a,b]),b_to_a=float(rate[b,a]),bidirectional_mean=ab,merge=merge))
   if merge:parent[find(b)]=find(a)
 components={}
 for i in range(k):components.setdefault(find(i),[]).append(i)
 remap=np.zeros(k,dtype=np.int32);compnames=[];members={}
 for j,ids in enumerate(components.values()):
  label='|'.join(names[ids]);compnames.append(label);members[label]=list(names[ids]);remap[ids]=j
 mv=remap[votes];freq=np.zeros((len(y),len(compnames)),dtype=np.int16)
 for s in range(mv.shape[1]):np.add.at(freq,(np.arange(len(y)),mv[:,s]),1)
 mode=freq.argmax(1);mx=freq.max(1);ties=(freq==mx[:,None]).sum(1)>1
 return np.array(compnames)[mode],mx/seeds.shape[1],ties,pairs,mapping,members

def main(h,r):
 assert json.loads((r/'input_audit_v1/acceptance.json').read_text())['status']=='passed'
 out=r/'consensus_v1';out.mkdir(exist_ok=False)
 L=pd.read_csv(h/'p0/cell_labels_r3.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('cell_id')
 R2=pd.read_csv(h/'p1/provenance/repo/inventories/annotation_v2/lineage_r2_reviewed.tsv',sep='\t',dtype=str,keep_default_na=False)
 R3=pd.read_csv(r/'review_v1/lineage_r3_reviewed.tsv',sep='\t',dtype=str,keep_default_na=False)
 sums=[]
 for pool in ['T_NK','B_Plasma','Stromal','Myeloid','Endothelial']:
  rnd='r3' if pool in ['Myeloid','Endothelial'] else 'r2';labcol=rnd+'_label'
  rv=(R3 if rnd=='r3' else R2);rv=rv[rv.pool==pool].set_index('cluster')
  s=pd.read_csv(h/f'p0/consensus/{pool}.seeds.tsv.gz',sep='\t',dtype=str).set_index('cell_id');assert s.index.is_unique
  original=s.run_r1.map(rv[labcol]);actions=s.run_r1.map(rv.action);assert original.notna().all()
  sc=[f'seed_{n}' for n in range(2026110800,2026110810)];assert not s[sc].isna().any().any()
  label,stability,ties,pairs,maps,members=consensus(original.values,s[sc].values)
  ll=L.loc[s.index];low=(ll.r1_action=='keep_flag_low_depth')|(actions=='keep_flag_low_depth')
  retained=(ll.major_label=='r3_'+pool)&actions.str.startswith('keep')
  # A merged component containing any nonidentity action is not accepted as an identity.
  allowed=set(rv.loc[rv.action.isin(['keep','keep_tissue']),labcol]);bad=set(rv.loc[~rv.action.isin(['keep','keep_tissue']),labcol]);allowed-=bad
  safe={name:all(x in allowed for x in vals) for name,vals in members.items()}
  eligible=retained&~low&actions.isin(['keep','keep_tissue'])&(stability>=.7)&~ties&pd.Series(label,index=s.index).map(safe)
  d=pd.DataFrame({'cell_id':s.index,'pool':pool,'review_round':rnd,'cluster':s.run_r1.values,'review_label':original.values,'review_action':actions.values,'consensus_label':label,'stability':stability,'tied_mode':ties,'low_depth':low.values,'retained_pool':retained.values,'include_main_subtype':eligible.values})
  d.to_csv(out/f'{pool}.cells.tsv.gz',sep='\t',index=False,compression='gzip')
  pd.DataFrame(pairs).to_csv(out/f'{pool}.pairs.tsv',sep='\t',index=False);pd.DataFrame(maps).to_csv(out/f'{pool}.seed_mapping.tsv',sep='\t',index=False)
  (out/f'{pool}.components.json').write_text(json.dumps(members,indent=2))
  d.groupby(['consensus_label','review_action','include_main_subtype']).size().rename('cells').reset_index().to_csv(out/f'{pool}.summary.tsv',sep='\t',index=False)
  sums.append(dict(pool=pool,cells=len(d),retained=int(retained.sum()),stable=int((stability>=.7).sum()),main_subtype_cells=int(eligible.sum()),components=len(members),merged_pairs=sum(x['merge'] for x in pairs),tied_modes=int(ties.sum())))
  print(json.dumps(sums[-1]),flush=True)
 pd.DataFrame(sums).to_csv(out/'summary.tsv',sep='\t',index=False)
 (out/'DONE.json').write_text(json.dumps({'status':'computed_not_final_annotation','n_seeds':10,'merge_threshold':.25,'stable_threshold':.7,'majority_ties':'lexicographic_tiebreak_logged','modal_ties':'exclude_main','exchange_denominator':'original_label_cell_count_times_10','scope':'all reviewed labels; nonidentity components excluded from main'},indent=2))
if __name__=='__main__':main(pathlib.Path(sys.argv[1]),pathlib.Path(sys.argv[2]))
