import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
"""Supplement legacy-unscored recorded QC units using full-gene raw counts."""
import os,json,hashlib,warnings,time
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np,pandas as pd,scipy.sparse as sp,zarr,h5py
from anndata.io import read_elem
import scrublet as scr

def fit_unit(job):
 key,x,positions=job;results=[];scores=[];calls=[];sims=[]
 for seed in [0,37,20261106]:
  try:
   with warnings.catch_warnings(record=True) as ww:
    model=scr.Scrublet(x,expected_doublet_rate=.04,random_state=seed);sc,pr=model.scrub_doublets(min_counts=2,min_cells=3,min_gene_variability_pctl=85,n_prin_comps=min(30,max(5,x.shape[0]//5)),use_approx_neighbors=False,verbose=False)
   assert pr is not None and np.isfinite(sc).all() and np.isfinite(model.threshold_)
   scores.append(np.asarray(sc));calls.append(np.asarray(pr,bool));sims.append(np.asarray(model.doublet_scores_sim_));results.append({'seed':seed,'status':'completed','threshold':float(model.threshold_),'detected_fraction':float(np.mean(pr)),'detectable_fraction':float(model.detectable_doublet_fraction_),'n_PC':int(model.manifold_obs_.shape[1]),'n_features':int(len(model._gene_filter)),'warnings':[str(w.message)[:160] for w in ww]})
  except Exception as e:results.append({'seed':seed,'status':'failed','error':str(e)[:300]})
 valid=len(scores)==3
 if valid:valid=all(z['detectable_fraction']>=.1 and z['detected_fraction']<=.15 for z in results) and np.ptp([z['detected_fraction'] for z in results])<=.05
 vote=np.stack(calls).sum(axis=0) if len(calls)==3 else np.zeros(x.shape[0],int);median=np.median(np.stack(scores),axis=0) if len(scores)==3 else np.full(x.shape[0],np.nan)
 return key,positions,results,valid,vote,median,scores,sims

def main():
 root=PROJECT_ROOT;run=root/'runs/20261106_v7_QC_clean_reannotation_v1';o=run/'stage_b_unscored_supplement';o.mkdir(exist_ok=False);started=datetime.now(timezone.utc).isoformat();g=zarr.open_group(str(root/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r')
 meta=pd.DataFrame({c:read_elem(g['obs/'+c]) for c in ['gse_id','sampleID','doublet_method','doublet_status']});index=pd.Index(read_elem(g['obs/_index']));take=meta.doublet_status.astype(str).eq('not_run').to_numpy();positions=np.flatnonzero(take);q=meta.iloc[positions].copy();q.index=index[positions];q['global_position']=positions;q['supplement_status']='pending';q['supplement_votes']=0;q['supplement_score_median']=np.nan;q['supplement_operational_exclude']=False
 jobs=[];records=[];small_units=[]
 per=root/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'
 assert len(q)==28999
 for gs,z in q.groupby('gse_id',observed=True):
  with h5py.File(per/(str(gs)+'_symbol_counts.h5ad'),'r') as f:
   obs=f['obs'];ids=pd.Index(read_elem(obs[obs.attrs['_index']]));mapping=ids.get_indexer(z.index);assert (mapping>=0).all();m=f['X'];ptr=m['indptr'][:]
   for key,sub in z.groupby(['sampleID','doublet_method'],observed=True,dropna=False):
    unit=hashlib.sha256((str(gs)+'|'+str(key)).encode()).hexdigest();rec={'unit_hash':unit,'gse_id':str(gs),'n_cells':len(sub)}
    if len(sub)<100:
     q.loc[sub.index,'supplement_status']='not_evaluable_lt100';rec.update(status='not_evaluable_lt100',n_excluded=0);records.append(rec);small_units.append(unit);continue
    rr=ids.get_indexer(sub.index);outptr=np.r_[0,np.cumsum(ptr[rr+1]-ptr[rr])].astype(np.int64);dd=[];ii=[]
    for j in rr:
     lo,hi=int(ptr[j]),int(ptr[j+1]);dd.append(m['data'][lo:hi]);ii.append(m['indices'][lo:hi])
    x=sp.csr_matrix((np.concatenate(dd),np.concatenate(ii),outptr),shape=(len(rr),int(m.attrs['shape'][1])));assert (x.data>=0).all();assert np.all(x.data==np.round(x.data));jobs.append((unit,x,q.index.get_indexer(sub.index)));rec['status']='ready';records.append(rec)
 (o/'started.json').write_text(json.dumps({'UTC_start':started,'n_cells':len(q),'n_units':len(records),'n_eligible_units':len(jobs),'n_lt100_units':len(small_units),'seeds':[0,37,20261106],'threads_per_worker':1,'workers':4,'model_scope':'Full recorded QC unit, not lineage subset.','min_cells':100,'new_mask_rule':'Three-seed consensus in numerically valid, stable units; operational QC decision subject to score-histogram review.'},indent=2))
 bykey={z['unit_hash']:z for z in records};hist=[];done=0
 with ProcessPoolExecutor(max_workers=4) as pool:
  futures=[pool.submit(fit_unit,j) for j in jobs]
  for fut in as_completed(futures):
   key,pos,res,valid,votes,med,scs,sim=fut.result();ids=q.index[pos];q.loc[ids,'supplement_votes']=votes;q.loc[ids,'supplement_score_median']=med;q.loc[ids,'supplement_status']='scored_stable_unit' if valid else 'model_or_threshold_review';q.loc[ids,'supplement_operational_exclude']=valid&(votes==3);rec=bykey[key];rec.update(status='scored_stable_unit' if valid else 'model_or_threshold_review',n_excluded=int((valid&(votes==3)).sum()),seed_results=res)
   for seed,sc,si in zip([z['seed'] for z in res if z['status']=='completed'],scs,sim):
    bins=np.linspace(0,1,51);hist.append({'unit_hash':key,'seed':seed,'observed_counts':np.histogram(sc,bins=bins)[0].tolist(),'simulated_counts':np.histogram(si,bins=bins)[0].tolist(),'bin_edges':bins.tolist()})
   done+=1;(o/'progress.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'units_done':done,'units_total':len(jobs),'n_provisional_consensus_exclusions':int(q.supplement_operational_exclude.sum())}))
 q.to_pickle(o/'unscored_supplement_PRIVATE.pkl.gz');(o/'unit_fit_records_PRIVATE.json').write_text(json.dumps(records,indent=2));(o/'score_histograms_PRIVATE.json').write_text(json.dumps(hist));agg=q.groupby(['gse_id','supplement_status','supplement_operational_exclude'],observed=True,dropna=False).size().reset_index(name='n_cells');agg.to_csv(o/'supplement_by_source.tsv',sep='\t',index=False);status=q.supplement_status.value_counts().to_dict();result={'status':'completed_legacy_unscored_supplement_candidate_review','UTC_start':started,'UTC_end':datetime.now(timezone.utc).isoformat(),'n_cells':len(q),'n_units':len(records),'n_fitted_units':len(jobs),'status_counts':status,'n_consensus_exclusion_candidates':int(q.supplement_operational_exclude.sum()),'histogram_adjudication_pending':True,'old_scores_and_labels_modified':False};(o/'complete.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
if __name__=='__main__':main()
