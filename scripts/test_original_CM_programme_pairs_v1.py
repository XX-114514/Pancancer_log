"""Frozen original CM pairs: block-level expression-programme associations, not LR proof."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd
from scipy.stats import rankdata,pearsonr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def residual(x,c):
 q=np.linalg.qr(c,mode='reduced')[0];return x-q@(q.T@x)
def corrmat(x,y):
 den=np.linalg.norm(x,axis=0)[:,None]*np.linalg.norm(y,axis=0)[None,:]
 return np.divide(x.T@y,den,out=np.full(den.shape,np.nan),where=den>1e-18)
def scorecorr(x,y,c):
 xx=rankdata(x,axis=0);yy=rankdata(y,axis=0);xx-=xx.mean(axis=0);yy-=yy.mean(axis=0)
 nx=np.linalg.norm(xx,axis=0);ny=np.linalg.norm(yy,axis=0)
 if c is not None:xx=residual(xx,c);yy=residual(yy,c)
 xx[:,np.linalg.norm(xx,axis=0)<=1e-12*np.maximum(nx,1)]=0;yy[:,np.linalg.norm(yy,axis=0)<=1e-12*np.maximum(ny,1)]=0
 # Exact-constant columns must not become spurious floating point residuals.
 xx[:,np.ptp(x,axis=0)==0]=0;yy[:,np.ptp(y,axis=0)==0]=0
 return corrmat(xx,yy),xx,yy
def bh(p):
 out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p));o=ix[np.argsort(p[ix],kind='stable')];q=p[o]*len(o)/np.arange(1,len(o)+1);out[o]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);return out
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path('/tmp/original_CM_programme_pair_plan_v1.json');plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)==plan['dictionary_sha256'];rec={r['programme_id']:r for r in json.loads(D.read_text())['programmes']};A=V/'all_source_native_UCell_audit_aggregation_20261014_v2';S=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';ai=json.loads((A/'inventory.json').read_text());si=json.loads((S/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL4532933') and si['audit_inventory_sha256']==sha(A/'inventory.json');assert sha(S/'inventory.json')==plan['sensitivity_inventory_sha256'];F=R/'all_CM_state_native_input_preflight_20261014_v1';fi=json.loads((F/'inventory.json').read_text());U=F/'private_all_frozen_sample_node_units.tsv.gz';assert sha(U)==fi['private_unit_table_sha256'];u=pd.read_csv(U,sep='\t',keep_default_na=False);assert len(u)==35567
 E=AR/'20261010_expression_annotation_downstream_v1/covarnet/tables/covarnet_each_edges.tsv';assert sha(E)==plan['edges_sha256'];edges=pd.read_csv(E,sep='\t');pairs=sorted({tuple(sorted([r.subCluster1,r.subCluster2])) for r in edges.itertuples()});assert len(pairs)==60
 # All identities retained; no programme selection by association statistics.
 arrays={};major_for={};position={};totals=u.groupby(['sample_uid','major']).n_cells.sum();variants=['scores_maxrank1000_mean','scores_maxrank3000_mean','TPM_signature_scores_maxrank1500_mean']
 for af in ai['private_aggregate_files']:
  p=A/af['relative_file'];sf=next(f for f in si['private_outputs'] if f['major']==af['major']);sp=S/sf['relative_file'];assert sha(p)==af['sha256'] and sha(sp)==sf['sha256'];a=dict(np.load(p,allow_pickle=False));s=dict(np.load(sp,allow_pickle=False));assert np.array_equal(a['unit_index'],s['unit_index']) and np.array_equal(a['programme_ids'],s['programme_ids']);ids=np.array(a['programme_ids']);keep=np.array([rec[str(x)]['algorithm_stable_candidate'] for x in ids],bool);idx=a['unit_index'];md=u.iloc[idx];assert md.major.eq(af['major']).all();assert np.array_equal(a['cells'],md.n_cells);assert np.allclose(a['programme_mean'],s['scores_maxrank1500_mean'],atol=1e-10,equal_nan=True)
  arrays[af['major']]=dict(ids=ids[keep],score=a['programme_mean'][:,keep],valid=a['valid_score_cells'][:,keep],support=a['informative_gene_support_cells'][:,keep],RNA=s['mean_full_RNA'],detected=s['mean_full_detected_genes'],variants={v:s[v][:,keep] for v in variants})
  for j,k in enumerate(idx):major_for[int(k)]=af['major'];position[int(k)]=j
 O=V/'original_CM_programme_pair_associations_20261015_v1';assert not O.exists();O.mkdir(mode=0o700);results=[];coverage=[];audits=[];rng=np.random.default_rng(plan['seed']);context_index=0
 for pair_index,(na,nb) in enumerate(pairs):
  left=u[u.subCluster.eq(na)].copy();right=u[u.subCluster.eq(nb)].copy();assert left.sample_uid.is_unique and right.sample_uid.is_unique;joined=left.merge(right,on=['sample_uid','gse_id','stratum','effective_resampling_block'],suffixes=('_A','_B'),validate='one_to_one');mods=sorted(set(edges.loc[((edges.subCluster1==na)&(edges.subCluster2==nb))|((edges.subCluster1==nb)&(edges.subCluster2==na)),'cm']))
  for (sid,stratum),jall in joined.groupby(['gse_id','stratum'],sort=True):
   j=jall[(jall.n_cells_A>=20)&(jall.n_cells_B>=20)&jall.analysis_eligible_A&jall.analysis_eligible_B].copy();nlibs=len(j);blocks=j.effective_resampling_block.to_numpy(str);bkeys=sorted(set(blocks));n=len(bkeys);context_index+=1;cov=dict(context=context_index,pair_index=pair_index,endpoint_A=na,endpoint_B=nb,CMs=';'.join(mods),gse_id=sid,stratum=stratum,paired_libraries_before_cell_gate=len(jall),paired_libraries=nlibs,analysis_blocks=n,reviewed_patient_blocks=int(j.loc[j.unit_kind_A.eq('reviewed_patient'),'effective_resampling_block'].nunique()),status='PENDING');coverage.append(cov)
   if n<plan['min_blocks']:cov['status']='INSUFFICIENT_BLOCKS';continue
   ma=str(j.major_A.iloc[0]);mb=str(j.major_B.iloc[0]);aa=arrays[ma];bb=arrays[mb];ia=np.array([position[int(k)] for k in j.unit_index_A]);ib=np.array([position[int(k)] for k in j.unit_index_B]);xx=aa['score'][ia];yy=bb['score'][ib];ka=(aa['valid'][ia]>=20).all(axis=0)&np.isfinite(xx).all(axis=0);kb=(bb['valid'][ib]>=20).all(axis=0)&np.isfinite(yy).all(axis=0);cov.update(stable_programmes_A=len(ka),stable_programmes_B=len(kb),measured_programmes_A=int(ka.sum()),measured_programmes_B=int(kb.sum()))
   if not ka.any() or not kb.any():cov['status']='GENE_OR_SCORE_COVERAGE_INSUFFICIENT';continue
   xa=xx[:,ka];yb=yy[:,kb];ida=aa['ids'][ka];idb=bb['ids'][kb];groups=[np.flatnonzero(blocks==k) for k in bkeys];collapse=lambda z:np.stack([np.mean(z[t],axis=0) for t in groups]);x=collapse(xa);y=collapse(yb)
   # Equal library means within each reviewed-patient/sample-fallback block.
   pa=np.array([nc/totals.loc[(sample,ma)] for nc,sample in zip(j.n_cells_A,j.sample_uid)]);pb=np.array([nc/totals.loc[(sample,mb)] for nc,sample in zip(j.n_cells_B,j.sample_uid)]);qc=np.column_stack([np.log1p(j.n_cells_A.to_numpy(float)),np.log1p(j.n_cells_B.to_numpy(float)),pa,pb,np.log1p(aa['RNA'][ia]),np.log1p(bb['RNA'][ib]),np.log1p(aa['detected'][ia]),np.log1p(bb['detected'][ib])]);z=collapse(qc);sd=z.std(axis=0)
   if (sd==0).any():cov['status']='CONSTANT_PRESPECIFIED_COVARIATE';continue
   c=np.column_stack([np.ones(n),(z-z.mean(axis=0))/sd]);rank=int(np.linalg.matrix_rank(c));cov['design_rank']=rank;cov['design_columns']=c.shape[1]
   if rank<c.shape[1]:cov['status']='COLLINEAR_PRESPECIFIED_DESIGN';continue
   if n-rank-1<plan['min_residual_df']:cov['status']='INSUFFICIENT_RESIDUAL_DF';continue
   obs,ex,ey=scorecorr(x,y,c);raw,_,_=scorecorr(x,y,None);valid=np.isfinite(obs);cov['finite_programme_pairs']=int(valid.sum());cov['status']='ESTIMABLE' if valid.any() else 'CONSTANT_PROGRAMME_SCORES'
   if not valid.any():continue
   picks=np.unique(np.linspace(0,obs.size-1,min(3,obs.size),dtype=int));rx=rankdata(x,axis=0);ry=rankdata(y,axis=0)
   for pick in picks:
    a,b=divmod(int(pick),y.shape[1]);vx=rx[:,a]-c@np.linalg.lstsq(c,rx[:,a],rcond=None)[0];vy=ry[:,b]-c@np.linalg.lstsq(c,ry[:,b],rcond=None)[0];v=float(pearsonr(vx,vy).statistic) if np.isfinite(obs[a,b]) else np.nan;assert np.isclose(v,obs[a,b],atol=1e-8,equal_nan=True);audits.append(dict(context=context_index,programme_A=str(ida[a]),programme_B=str(idb[b]),OLS_partial_rank=v,QR_partial_rank=float(obs[a,b])))
   exceed=np.zeros(obs.shape,np.int64);fitted=rx-ex
   for k in range(plan['permutations']):
    perm=rng.permutation(n);pp=residual(fitted+ex[perm],c);pp[:,np.ptp(x,axis=0)==0]=0;v=corrmat(pp,ey);exceed+=np.isfinite(v)&(np.abs(v)>=np.abs(obs)-1e-12)
   pv=(exceed+1)/(plan['permutations']+1);pv[~valid]=np.nan
   boot=np.full((plan['bootstraps'],*obs.shape),np.nan,np.float32)
   for k in range(plan['bootstraps']):
    take=rng.integers(0,n,n);cc=c[take]
    if np.linalg.matrix_rank(cc)<rank:continue
    boot[k]=scorecorr(x[take],y[take],cc)[0]
   bn=np.isfinite(boot).sum(axis=0);ci=np.full((2,*obs.shape),np.nan)
   for a,b in zip(*np.where(bn>=plan['min_valid_bootstraps'])):ci[:,a,b]=np.quantile(boot[:,a,b][np.isfinite(boot[:,a,b])],[.025,.975])
   sens={}
   for vk in variants:
    xv=collapse(aa['variants'][vk][ia][:,ka]);yv=collapse(bb['variants'][vk][ib][:,kb]);ga=np.isfinite(xv).all(axis=0);gb=np.isfinite(yv).all(axis=0);mat=np.full(obs.shape,np.nan)
    if ga.any() and gb.any():mat[np.ix_(ga,gb)]=scorecorr(xv[:,ga],yv[:,gb],c)[0]
    sens[vk]=mat
   sa=np.median(collapse(aa['support'][ia][:,ka]/aa['valid'][ia][:,ka]),axis=0);sb=np.median(collapse(bb['support'][ib][:,kb]/bb['valid'][ib][:,kb]),axis=0)
   for a,b in zip(*np.where(valid)):
    ra=rec[str(ida[a])];rb=rec[str(idb[b])];r=dict(context=context_index,pair_index=pair_index,endpoint_A=na,endpoint_B=nb,CMs=';'.join(mods),gse_id=sid,stratum=stratum,programme_A=str(ida[a]),programme_B=str(idb[b]),analysis_blocks=n,paired_libraries=nlibs,reviewed_patient_blocks=cov['reviewed_patient_blocks'],raw_spearman=float(raw[a,b]),adjusted_partial_rank=float(obs[a,b]),permutation_p=float(pv[a,b]),CI_low=float(ci[0,a,b]),CI_high=float(ci[1,a,b]),valid_bootstraps=int(bn[a,b]),median_support_A=float(sa[a]),median_support_B=float(sb[b]),foreign_RNA_review_A=ra.get('foreign_RNA_sensitive',False),foreign_RNA_review_B=rb.get('foreign_RNA_sensitive',False),basic_expression_review_A=ra.get('basic_expression_review_flag',False),basic_expression_review_B=rb.get('basic_expression_review_flag',False),formal_state=False,LR_mechanism_tested=False)
    for vk in variants:r[vk+'_adjusted_rho']=float(sens[vk][a,b])
    results.append(r)
   (O/'progress.json').write_text(json.dumps(dict(contexts_seen=context_index,estimable_contexts=sum(r['status']=='ESTIMABLE' for r in coverage),tests=len(results),observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');del boot;gc.collect()
 out=pd.DataFrame(results)
 if len(out):out['BH_q']=bh(out.permutation_p.to_numpy(float))
 out.to_csv(O/'all_programme_pair_associations.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'all_context_coverage_reasons.tsv',sep='\t',index=False);pd.DataFrame(audits).to_csv(O/'independent_scalar_OLS_checks.tsv',sep='\t',index=False);inv=dict(status='COMPLETE_FIXED_ORIGINAL60_STATE_PROGRAMME_PAIRS_NOT_LR_MECHANISM',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),dictionary_sha256=sha(D),pairs=60,contexts=len(coverage),estimable_contexts=sum(r['status']=='ESTIMABLE' for r in coverage),tests=len(out),main_BH_records=int((out.BH_q<.05).sum()) if len(out) else 0,permutations=plan['permutations'],bootstraps=plan['bootstraps'],formal_states=False,LR_mechanism_tested=False,output_sha256={p.name:sha(p) for p in O.glob('*.tsv')});(O/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n');print(json.dumps(inv))
if __name__=='__main__':main()
