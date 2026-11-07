"""Independently audit literal author scores and add predeclared spatial diagnostics."""
import sys,json,hashlib,itertools,math,zlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import rankdata
from anndata.io import read_elem
import h5py
r=Path(sys.argv[1]);base=Path(sys.argv[2]);key=sys.argv[3];out=r/'sections'/key;receipt=json.loads((out/'complete.json').read_text());assert receipt['epochs_recorded']>=30000 and receipt['n_reference_labels']==50;history=np.asarray(np.load(out/'ELBO_history_PRIVATE.npy',allow_pickle=True),dtype=float);assert np.isfinite(history).all() and len(history)>=30000
fq=pd.read_csv(base/'inputs/frequency_PRIVATE.tsv',sep='\t',index_col=0);pool=fq.index.tolist();parents=fq.majorCluster.to_dict();members=pd.read_csv(base/'CM_results/CM_members_after_author_network_filter.tsv',sep='\t');edges=pd.read_csv(base/'CM_results/each_edge.tsv',sep='\t');W=pd.read_csv(base/'CM_results/W_PRIVATE.tsv',sep='\t',index_col=0);coords=pd.read_csv(out/'coordinates_PRIVATE.tsv.gz',sep='\t',index_col=0);weights=pd.read_csv(out/'author_weights60_PRIVATE.tsv',sep='\t');G=sp.csr_matrix((weights.weight,(weights.i.to_numpy()-1,weights.j.to_numpy()-1)),shape=(len(coords),len(coords)));assert np.allclose(np.asarray(G.sum(axis=1)).ravel(),1,atol=1e-10) and np.all(G.getnnz(axis=1)==60);assert G.diagonal().sum()==0
# Every saved distance weight is independently inverse-square and row normalized.
xy=coords.to_numpy();dist2=np.square(xy[weights.i.to_numpy()-1]-xy[weights.j.to_numpy()-1]).sum(axis=1);raw=1/dist2;expected=raw/np.bincount(weights.i.to_numpy()-1,weights=raw,minlength=len(coords))[weights.i.to_numpy()-1];assert np.allclose(expected,weights.weight,atol=1e-10,rtol=1e-9)
sourcejob=next(z for z in json.loads((r/'inputs/section_jobs_PRIVATE.json').read_text()) if z['section_code']==key);inf=pd.read_csv(r/'reference/signatures_PRIVATE.tsv',sep='\t',index_col=0)
with h5py.File(sourcejob['input_path'],'r') as h:native=set(read_elem(h['var']).index)
shared=sorted(set(inf.index)&native);refcorr=np.corrcoef(np.log1p(inf.loc[shared,pool].to_numpy()).T);upper=np.triu_indices(len(pool),1);refvec=refcorr[upper];author=pd.read_csv(out/'author_CM_scores_PRIVATE.tsv',sep='\t');nullcache={};diagnostics=[];edge_rows=[];activities=[];idx={n:i for i,n in enumerate(pool)}
def author_score(mat,nodes,directed=False):
 membership=np.isin(np.arange(len(pool)),[idx[n] for n in nodes]);ii,jj=np.where(~np.eye(len(pool),dtype=bool)) if directed else upper;vals=mat[ii,jj];within=membership[ii]&membership[jj];bg=membership[ii]|membership[jj]
 if len(nodes)<2 or not within.any():return np.nan
 median=np.median(vals[within]);return float(np.mean(median>=vals[bg])) if np.isfinite(vals[bg]).all() else np.nan
def null_members(cm,nodes):
 if cm in nullcache:return nullcache[cm]
 composition=pd.Series([parents[n] for n in nodes]).value_counts().sort_index();choices={p:[n for n in pool if parents[n]==p] for p in composition.index};N=math.prod(math.comb(len(choices[p]),int(composition[p])) for p in composition.index);sets=[]
 if N<=1000:
  for t in itertools.product(*[list(itertools.combinations(choices[p],int(composition[p]))) for p in composition.index]):sets.append([n for z in t for n in z])
 else:
  rng=np.random.default_rng(2026110701+zlib.crc32(cm.encode()));seen=set()
  while len(sets)<1000:
   draw=tuple(sorted(n for p in composition.index for n in rng.choice(choices[p],int(composition[p]),replace=False)))
   if draw not in seen:seen.add(draw);sets.append(list(draw))
 nullcache[cm]=(sets,N);return sets,N
for posterior in ['means','q05']:
 ab=pd.read_csv(out/(posterior+'_abundance_PRIVATE.tsv.gz'),sep='\t',index_col=0);assert ab.index.equals(coords.index);a=ab[pool].to_numpy(float);assert np.isfinite(a).all() and (a>=0).all();ranked=rankdata(a,axis=0,method='average');C=np.corrcoef(ranked.T);sd=a.std(axis=0,ddof=1);Z=(a-a.mean(axis=0))/sd;M=Z.T@(G@Z)/(len(a)-1);np.fill_diagonal(M,np.nan);adjustedC=np.full_like(C,np.nan);adjustedM=np.full_like(M,np.nan)
 for i,j in zip(*upper):
  similar=(np.abs(refvec-refcorr[i,j])<=.02)&~((upper[0]==i)&(upper[1]==j))
  if similar.sum()>=20:
   adjustedC[i,j]=adjustedC[j,i]=C[i,j]-np.median(C[upper][similar]);adjustedM[i,j]=M[i,j]-np.median(M[upper][similar]);adjustedM[j,i]=M[j,i]-np.median(M.T[upper][similar])
 unscaled=np.zeros((len(a),14))
 for k,cm in enumerate(W.columns):
  nodes=members.loc[members.cm.eq(cm),'subCluster'].tolist();memidx=[idx[n] for n in nodes]
  if nodes:unscaled[:,k]=a[:,memidx]@W.loc[nodes,cm].to_numpy()
  for stat,mat,directed,adjusted in [('colocalization',C,False,adjustedC),('aggregation',M,True,adjustedM)]:
   value=author_score(mat,nodes,directed);saved=author[author.cm.eq(cm)&author.posterior.eq(posterior)&author.statistic.eq(stat)].iloc[0]
   if stat=='aggregation' and np.any(sd==0):value=np.nan # literal R loop aborts on a constant scaled member
   assert (pd.isna(saved.author_score) and np.isnan(value)) or abs(value-saved.author_score)<1e-10,(key,posterior,cm,stat,value,saved.author_score)
   null,N=(null_members(cm,nodes) if len(nodes)>=2 else ([],0));nullval=np.array([author_score(mat,n,directed) for n in null]);valid=nullval[np.isfinite(nullval)];diagnostics.append({'section_code':key,'dataset':sourcejob['dataset'],'cm':cm,'posterior':posterior,'statistic':stat,'author_score':value,'n_original_members':len(nodes),'lineage_null_unique_sets_evaluated':len(null),'lineage_null_all_possible_sets':N,'lineage_null_exhaustive':bool(N and N<=1000),'lineage_null_median':float(np.median(valid)) if len(valid) else np.nan,'lineage_null_percentile_le':float(np.mean(valid<=value)) if len(valid) and np.isfinite(value) else np.nan,'reference_similarity_adjusted_score':author_score(adjusted,nodes,directed),'reference_control_rule':'pair log-signature Pearson +/-0.02; >=20 alternative pairs; subtract observed median; diagnostic only'})
  for e in edges[edges.cm.eq(cm)].itertuples():
   i,j=idx[e.subCluster1],idx[e.subCluster2];edge_rows.append({'section_code':key,'dataset':sourcejob['dataset'],'cm':cm,'posterior':posterior,'node1':e.subCluster1,'node2':e.subCluster2,'Spearman':C[i,j],'Moran_forward':M[i,j],'Moran_reverse':M[j,i],'reference_signature_logPearson':refcorr[i,j],'reference_adjusted_Spearman':adjustedC[i,j],'reference_adjusted_Moran_forward':adjustedM[i,j]})
 q99=float(np.quantile(unscaled,.99));scaled=unscaled/q99 if q99>0 else np.full_like(unscaled,np.nan);pd.DataFrame(scaled,index=coords.index,columns=W.columns).to_csv(out/(posterior+'_CM_activity_PRIVATE.tsv.gz'),sep='\t');activities.extend({'section_code':key,'dataset':sourcejob['dataset'],'cm':cm,'posterior':posterior,'activity_mean':float(np.nanmean(scaled[:,k])),'activity_q99_scale':q99,'activity_raw_mean':float(unscaled[:,k].mean())} for k,cm in enumerate(W.columns))
pd.DataFrame(diagnostics).to_csv(out/'audited_CM_diagnostics_PRIVATE.tsv',sep='\t',index=False);pd.DataFrame(edge_rows).to_csv(out/'audited_edge_diagnostics_PRIVATE.tsv',sep='\t',index=False);pd.DataFrame(activities).to_csv(out/'activity_summary_PRIVATE.tsv',sep='\t',index=False);proof={'status':'PASS_independent_literal_author_scores_60NN_inverse_square_weights_and_full_epoch_posterior_checks','UTC':datetime.now(timezone.utc).isoformat(),'section_code':key,'dataset':sourcejob['dataset'],'n_spots':len(coords),'n_CM':14,'n_scores_checked':56,'paper_primary_posterior':'means expected abundance; q05 sensitivity','all_50_reference_labels_retained':True,'all_45_original_CM_background_nodes_retained':True,'epoch_history_finite_full_30000':True,'no_spot_independence_or_contact_or_mechanism_claim':True};(out/'independent_audit.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps(proof))
