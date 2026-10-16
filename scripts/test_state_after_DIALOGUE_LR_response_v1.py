"""All LR x supported receiver MCP; exploratory full-block fit, internal holdout diagnostic."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict,Counter
import os,json,hashlib,time,gc
import numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def residual(x,c):
 q=np.linalg.qr(c,mode='reduced')[0];return x-q@(q.T@x)
def scorr(x,y,c=None):
 if x.ndim==1:x=x[:,None]
 if y.ndim==1:y=y[:,None]
 rx=rankdata(x,axis=0);ry=rankdata(y,axis=0);rx-=rx.mean(0);ry-=ry.mean(0);sx=np.linalg.norm(rx,axis=0);sy=np.linalg.norm(ry,axis=0)
 if c is not None:rx=residual(rx,c);ry=residual(ry,c)
 rx[:,(np.linalg.norm(rx,axis=0)<=1e-12*np.maximum(sx,1))|(np.ptp(x,axis=0)==0)]=0;ry[:,(np.linalg.norm(ry,axis=0)<=1e-12*np.maximum(sy,1))|(np.ptp(y,axis=0)==0)]=0
 den=np.linalg.norm(rx,axis=0)[:,None]*np.linalg.norm(ry,axis=0)[None,:];out=np.divide(rx.T@ry,den,out=np.full(den.shape,np.nan),where=den>1e-18);return out[:,0],rx,ry
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p));o=ix[np.argsort(p[ix],kind='stable')];out[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1]);return out

def main():
 V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);plan_path=Path('/tmp/state_after_DIALOGUE_LR_response_plan_v1.json');plan=json.loads(plan_path.read_text());assert sha(__file__)==plan['producer_sha256'];G=V/'state_after_DIALOGUE_gene_components_20261016_v1';start=time.monotonic()
 while not (G/'inventory.json').exists():
  e=json.loads((V/'state_after_DIALOGUE_gene_components_20261016_v1.execution.json').read_text());p=Path('/proc',str(e['child_pid']),'cmdline');assert p.exists() and e['producer'].encode() in p.read_bytes(),e;assert time.monotonic()-start<36*3600;time.sleep(20)
 gi=json.loads((G/'inventory.json').read_text())
 for n,h in gi['files'].items():assert sha(G/n)==h,n
 O=V/'state_after_DIALOGUE_LR_response_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
 L=V/'state_after_DIALOGUE_all_LR_preflight_20261016_v1';la=json.loads((L/'independent_audit.json').read_text());assert la['status']=='PASS_ALL129472_CACHE_PARTITION_NATIVE_GENE_AND_QR_DESIGN_RECORDS';assert la['producer_inventory_sha256']==sha(L/'inventory.json');E=R/'immune_same_name_family_source_review_v1/CM_edge_LR_extraction_v2';ei=json.loads((E/'inventory.json').read_text());assert sha(E/'inventory.json')==plan['LR_inventory_sha256'];A=pd.read_csv(L/'all28_cache_alignment_checks.tsv',sep='\t');Q=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';assert sha(Q/'inventory.json')==plan['QC_inventory_sha256'];qi=json.loads((Q/'inventory.json').read_text());qc={}
 for info in qi['private_outputs']:
  p=Q/info['relative_file'];assert sha(p)==info['sha256'];a=np.load(p,allow_pickle=False)
  for k,r,d in zip(a['unit_index'],a['mean_full_RNA'],a['mean_full_detected_genes']):qc[int(k)]=(float(r),float(d))
 components=pd.read_csv(G/'private_all_block_gene_components.tsv.gz',sep='\t');scope=pd.read_csv(V/'state_after_DIALOGUE_frozen_projection_20261016_v1/public_slot_scope.tsv',sep='\t');coverage=pd.read_csv(V/'original_CM_programme_pair_associations_20261015_v2/all_context_coverage_reasons.tsv',sep='\t').set_index('context');results=[];counts=Counter();checks=[];seen=0;rng=np.random.default_rng(plan['seed'])
 for ar in A.itertuples():
  ctx=int(ar.context);root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{ctx:04d}';ci=json.loads((root/'private_input_inventory.json').read_text());parts=pd.read_csv(root/'private_block_library_partition.tsv',sep='\t');reverse=ar.source==ci['endpoints'][1];side='A' if reverse else 'B';typ='type'+side;p=E/'private_group_arrays'/(ar.old_group+'.npz');assert sha(p)==ei['outputs'][str(p.relative_to(E))];z=np.load(p,allow_pickle=False);ri=pd.Index(z['samples'].astype(str)).get_indexer(parts.sample_uid.astype(str));assert (ri>=0).all() and list(z['patient_blocks'][ri])==list(parts.effective_resampling_block);unitS=parts['unit_index_B' if reverse else 'unit_index_A'];unitT=parts['unit_index_A' if reverse else 'unit_index_B'];qs=np.array([qc[int(k)] for k in unitS]);qt=np.array([qc[int(k)] for k in unitT]);cov=np.column_stack([np.log1p(z['source_n'][ri]),np.log1p(z['target_n'][ri]),z['source_within_major'][ri],z['target_within_major'][ri],np.log1p(qs[:,0]),np.log1p(qt[:,0]),np.log1p(qs[:,1]),np.log1p(qt[:,1])]);scores=z['scores'][ri];assert scores.shape==(len(parts),4624,4);observed=np.isfinite(scores[:,:,0]);assert np.array_equal(np.isfinite(scores),np.repeat(observed[:,:,None],4,axis=2));train=parts.partition.eq('discovery').to_numpy();test=parts.partition.eq('holdout').to_numpy();CMs=coverage.loc[ctx,'CMs']
  for mcp in ['MCP1','MCP2']:
   key0=(ctx,ar.source,ar.target,mcp);ss=scope[(scope.context==ctx)&scope.type.eq(typ)&scope.program.eq(mcp)]
   if len(ss)!=1 or ss.iloc[0].status!='FROZEN_RECONSTRUCTION_PASS':counts[key0+('NO_GENE_SUPPORTED_RECEIVER_SCORE',)]+=4624;seen+=4624;continue
   cc=components[(components.context==ctx)&components.type.eq(typ)&components.program.eq(mcp)];assert len(cc)>0;wt=cc.groupby('gene').weight.agg(['first','nunique']);assert wt['nunique'].eq(1).all();B=cc.pivot(index='sample',columns='gene',values='contribution').reindex(parts.effective_resampling_block);assert B.notna().all().all();weights=wt['first'].reindex(B.columns).to_numpy();mat=B.to_numpy();full=mat@weights;genes=np.array(B.columns,str);L1=float(np.abs(weights).sum());groups=defaultdict(list)
   for j,(lig,rec) in enumerate(zip(z['ligands'],z['receptors'])):
    seen+=1;mask=observed[:,j];n=int(mask.sum());remove=tuple(sorted(set(str(lig).split('_')+str(rec).split('_'))&set(genes)));keep=~np.isin(genes,remove);retained=float(np.abs(weights[keep]).sum()/L1)
    if not z['native_measured'][j]:reason='NATIVE_COMPLEX_SUBUNIT_UNMEASURED'
    elif n<20:reason='FEWER_THAN_20_OBSERVED_BLOCKS'
    elif n/len(parts)<.8:reason='OUTPUT_COVERAGE_BELOW_80PCT'
    elif keep.sum()<plan['minimum_retained_genes'] or retained<plan['minimum_retained_L1_fraction']:reason='RECEIVER_PROGRAMME_CHANGED_TOO_MUCH_BY_LR_GENE_MASK'
    elif (str(lig),str(rec)) in [('C1QA','CD93'),('INHBA','ACTR2')]:reason='DESCRIPTIVE_BINDING_CAUTION_NO_MAIN_INFERENCE'
    else:reason='GROUP_PENDING'
    if reason!='GROUP_PENDING':counts[key0+(reason,)]+=1;continue
    groups[(mask.tobytes(),remove)].append((j,retained,int(keep.sum())))
   for (maskbytes,remove),items in groups.items():
    mask=np.frombuffer(maskbytes,dtype=bool);ix=np.flatnonzero(mask);X=np.log1p(scores[ix][:,[v[0] for v in items],0]);keep=~np.isin(genes,remove);Y=(mat[:,keep]@weights[keep])[ix];cv=cov[ix];sd=cv.std(0);n=len(ix)
    if (sd<=1e-12).any():counts[key0+('CONSTANT_PRESPECIFIED_COVARIATE',)]+=len(items);continue
    D=np.column_stack([np.ones(n),(cv-cv.mean(0))/sd]);rank=int(np.linalg.matrix_rank(D))
    if rank!=9 or n-rank-1<10:counts[key0+('FULL_DESIGN_COLLINEAR_OR_RESIDUAL_DF_BELOW_10',)]+=len(items);continue
    obs,ex,ey=scorr(X,Y,D);raw,_,_=scorr(X,Y);unmasked,_,_=scorr(X,full[ix],D);valid=np.isfinite(obs)
    if not valid.any():counts[key0+('CONSTANT_ADJUSTED_LR_OR_RECEIVER_SCORE',)]+=len(items);continue
    for jj in np.unique(np.linspace(0,len(items)-1,min(3,len(items)),dtype=int)):
     rx=rankdata(X[:,jj]);ry=rankdata(Y);vx=rx-D@np.linalg.lstsq(D,rx,rcond=None)[0];vy=ry-D@np.linalg.lstsq(D,ry,rcond=None)[0];den=np.linalg.norm(vx)*np.linalg.norm(vy);v=float(vx@vy/den) if np.isfinite(obs[jj]) else np.nan;assert np.isclose(v,obs[jj],atol=1e-8,equal_nan=True);checks.append(dict(context=ctx,source=ar.source,target=ar.target,MCP=mcp,ligand=str(z['ligands'][items[jj][0]]),receptor=str(z['receptors'][items[jj][0]]),OLS_rho=v,QR_rho=float(obs[jj])))
    fitted=rankdata(X,axis=0)-ex;exceed=np.zeros(len(items),np.int64)
    for _ in range(plan['permutations']):
     pp=residual(fitted+ex[rng.permutation(n)],D);den=np.linalg.norm(pp,axis=0)*np.linalg.norm(ey[:,0]);null=np.divide(pp.T@ey[:,0],den,out=np.full(len(items),np.nan),where=den>1e-18);exceed+=np.isfinite(null)&(np.abs(null)>=np.abs(obs)-1e-12)
    pv=(1+exceed)/(plan['permutations']+1);pv[~valid]=np.nan;boot=np.full((plan['bootstraps'],len(items)),np.nan,np.float32)
    for b in range(plan['bootstraps']):
     take=rng.integers(0,n,n)
     if np.linalg.matrix_rank(D[take])==rank:boot[b]=scorr(X[take],Y[take],D[take])[0]
    bn=np.isfinite(boot).sum(0);ci0=np.full((2,len(items)),np.nan)
    for j in np.flatnonzero(bn>=plan['minimum_valid_bootstraps']):ci0[:,j]=np.quantile(boot[:,j][np.isfinite(boot[:,j])],[.025,.975])
    sens=[]
    for k in [1,2,3]:sens.append(scorr(np.log1p(scores[ix][:,[v[0] for v in items],k]) if k==1 else -scores[ix][:,[v[0] for v in items],k],Y,D)[0])
    # Internal-heldout diagnostic: nuisance fitted on discovery only; no heldout p/FDR claim.
    tr=train[ix];te=test[ix];held=np.full(len(items),np.nan);heldreason='FEWER_THAN_20_DISCOVERY_OR_10_HELDOUT_OBSERVED_BLOCKS'
    if tr.sum()>=20 and te.sum()>=10 and np.linalg.matrix_rank(D[tr])==9:
     by=np.linalg.lstsq(D[tr],Y[tr],rcond=None)[0];bx=np.linalg.lstsq(D[tr],X[tr],rcond=None)[0];held=scorr(X[te]-D[te]@bx,Y[te]-D[te]@by)[0];heldreason='DESCRIPTIVE_FROZEN_NUISANCE_INTERNAL_HOLDOUT_NOT_EXTERNAL_VALIDATION'
    for j,(lr,retained,ng) in enumerate(items):
     if not valid[j]:counts[key0+('CONSTANT_ADJUSTED_LR_OR_RECEIVER_SCORE',)]+=1;continue
     counts[key0+('ESTIMABLE_EXPLORATORY_FULL_BLOCK_ASSOCIATION',)]+=1;results.append(dict(context=ctx,gse_id=ar.gse_id,CMs=CMs,source=ar.source,target=ar.target,MCP=mcp,ligand=str(z['ligands'][lr]),receptor=str(z['receptors'][lr]),blocks=n,reviewed_patient_blocks=int(parts.loc[mask,'unit_kind_A'].eq('reviewed_patient').sum()) if 'unit_kind_A' in parts else int(z['reviewed_patient'][ri][mask].sum()),discovery_observed_blocks=int(tr.sum()),internal_holdout_observed_blocks=int(te.sum()),masked_genes=';'.join(remove),retained_programme_genes=ng,retained_L1_fraction=retained,raw_masked_rho=float(raw[j]),adjusted_masked_rho=float(obs[j]),adjusted_unmasked_rho=float(unmasked[j]),permutation_p=float(pv[j]),CI_low=float(ci0[0,j]),CI_high=float(ci0[1,j]),valid_bootstraps=int(bn[j]),lr_means_adjusted_rho=float(sens[0][j]),negative_specificity_rank_adjusted_rho=float(sens[1][j]),negative_magnitude_rank_adjusted_rho=float(sens[2][j]),internal_holdout_rho=float(held[j]),internal_holdout_status=heldreason,exploratory_learned_programme_same_cohort=True,independent_validation=False,causal_or_binding_validation=False))
    (O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),hypotheses_accounted=seen,estimable_results=len(results),current_context=ctx,source=ar.source,target=ar.target,MCP=mcp))+'\n');del boot;gc.collect()
 assert seen==258944 and sum(counts.values())==seen,(seen,sum(counts.values()));out=pd.DataFrame(results)
 if len(out):
  out['BH_q']=bh(out.permutation_p);finite=out.permutation_p.notna();assert np.allclose(out.loc[finite,'BH_q'],multipletests(out.loc[finite,'permutation_p'],method='fdr_bh')[1]);assert not out.duplicated(['context','source','target','MCP','ligand','receptor']).any()
 out.to_csv(O/'all_estimable_exploratory_LR_receiver_associations.tsv',sep='\t',index=False);pd.DataFrame([dict(context=k[0],source=k[1],target=k[2],MCP=k[3],status=k[4],records=v) for k,v in counts.items()]).to_csv(O/'all258944_hypothesis_scope_counts.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'independent_scalar_OLS_checks.tsv',sep='\t',index=False);inventory=dict(status='ALL258944_LR_RECEIVER_SLOTS_ACCOUNTED_EXPLORATORY_NOT_EXTERNAL_OR_CAUSAL_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(plan_path),component_inventory_sha256=sha(G/'inventory.json'),LR_preflight_audit_sha256=sha(L/'independent_audit.json'),hypotheses_accounted=seen,estimable_tests=len(out),BH_below_05=int((out.BH_q<.05).sum()) if len(out) else 0,files={p.name:sha(p) for p in O.iterdir() if p.is_file()},clinical_outcomes_used=False,formal_states=False);(O/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n');print(json.dumps(inventory))
if __name__=='__main__':main()
