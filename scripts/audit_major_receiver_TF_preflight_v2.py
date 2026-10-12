"""Separate raw-network reconstruction and all-score OLS audit; no producer import."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,ast
import numpy as np,pandas as pd
from scipy.sparse import load_npz,csr_matrix
from scipy.stats import t

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';O=P/'major_receiver_TF_preflight_v2';D=Path(os.environ['PANCANCER_TF_REFERENCE']);plan=json.loads(Path(os.environ['PANCANCER_TF_PLAN']).read_text());assert sha(__file__)==plan['auditor_sha256'] and not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text());native=json.loads((N/'inventory.json').read_text());assert sha(N/'inventory.json')==plan['major_inventory_sha256'] and sha(N/'independent_audit.json')==plan['major_audit_sha256'];coverage=pd.read_csv(O/'all_source_TF_coverage.tsv',sep='\t');assert not coverage.duplicated(['gse_id','variant','TF']).any() and len(coverage)==11*2*1163
 for n,h in plan['reference_files'].items():assert sha(D/n)==h
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 # Independently reconstruct duplicate inhibition-priority weights from raw rows.
 raw=pd.read_csv(D/'human_collectri_raw.tsv',sep='\t');weights={}
 for a,b,inh in raw[['source_genesymbol','target_genesymbol','is_inhibition']].itertuples(index=False,name=None):
  a='AP1' if a.upper().startswith(('JUN','FOS')) else 'NFKB' if a.upper().startswith(('REL','NFKB')) else a
  key=(a,b);weights[key]=min(weights.get(key,1),-1 if inh else 1)
 signed=pd.read_csv(D/'human_collectri_signed_complexes.tsv',sep='\t');assert dict(((a,b),w) for a,b,w in signed.itertuples(index=False,name=None))==weights
 TFs=sorted(set(a for a,b in weights));tc={a:j for j,a in enumerate(TFs)};lr=pd.read_csv(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']),sep='\t');lrgenes=set()
 for vals in lr.itertuples(index=False,name=None):
  for v in vals:lrgenes.update(v.split('_'))
 # Evaluate the actual official function bodies in a restricted numeric namespace.
 tree=ast.parse((D/'official_decoupler192_ULM_source.txt').read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['mat_cov','mat_cor','t_val','ulm']];assert len(nodes)==4;space={'np':np,'csr_matrix':csr_matrix,'t':t,'tqdm':lambda v,**kw:v};exec(compile(ast.Module(body=nodes,type_ignores=[]),'<official_decoupler192_numeric_functions>','exec'),space)
 reports=[];source_data={};all_values=0;official_values=0
 for sid in sorted(x['gse_id'] for x in inv['source_summary']):
  a=np.load(O/'private_source_TF_arrays'/(sid+'.npz'));genes=np.load(N/'private_source_major_arrays'/(sid+'_axes.npz'))['genes'];X=load_npz(N/'private_source_major_arrays'/(sid+'_meanlog.npz')).toarray().astype(float);PCT=load_npz(N/'private_source_major_arrays'/(sid+'_pct_positive.npz')).toarray();units=pd.read_csv(N/'private_source_major_arrays'/(sid+'_units.tsv'),sep='\t',keep_default_na=False);assert a['TFs'].tolist()==TFs and np.isfinite(X).all()
  nonempty=np.any(X!=0,axis=0);source_data[sid]=(units,a);checked=0;official_checked=0
  for variant,prefix in [('nonLR_primary','primary'),('full_gene_sensitivity','sensitivity')]:
   keep=nonempty & (np.array([g not in lrgenes for g in genes]) if prefix=='primary' else True);ids=np.flatnonzero(keep);assert np.array_equal(ids,a[prefix+'_gene_indices']);gg=genes[ids];cc={g:i for i,g in enumerate(gg)};edges=[(cc[b],tc[c],w) for (c,b),w in weights.items() if b in cc];W=csr_matrix(([w for _,_,w in edges],([i for i,_,_ in edges],[j for _,j,_ in edges])),shape=(len(ids),len(TFs)));cnt=np.asarray((W!=0).sum(0)).ravel();assert np.array_equal(cnt,a[prefix+'_target_counts']);eligible=cnt>=5;Y=X[:,ids];n=len(ids);yc=Y-Y.mean(1,keepdims=True);wc=W[:,eligible];ws=np.asarray(wc.sum(0)).ravel();wss=np.asarray(wc.multiply(wc).sum(0)).ravel()-ws*ws/n
   # All values from OLS coefficient/SSR rather than correlation/t-transform.
   cross=np.asarray(yc@wc);beta=cross/wss;SSR=np.sum(yc*yc,axis=1)[:,None]-beta*beta*wss;ols=beta/np.sqrt(SSR/(n-2)/wss);actual=a[prefix+'_scores'];assert np.isnan(actual[:,~eligible]).all();assert np.allclose(actual[:,eligible],ols,rtol=2e-8,atol=2e-8);checked+=int(ols.size)
   det=np.asarray((PCT[:,ids]>=.01).astype(np.int32)@(W!=0).astype(np.int32));assert np.array_equal(det,a[prefix+'_detected_target_counts']);cov=coverage[(coverage.gse_id==sid)&(coverage.variant==variant)].set_index('TF').loc[TFs];assert np.array_equal(cov.measured_nonempty_targets.to_numpy(int),cnt);assert np.array_equal(cov.min5_scored.to_numpy(bool),eligible);assert np.array_equal(cov.min20_target_gate.to_numpy(bool),cnt>=20);assert np.allclose(cov.median_detected_targets_ge1pct.to_numpy(float),np.median(det,axis=0));assert cov.gene_universe.eq(n).all() and cov.selected_major_units.eq(len(Y)).all()
   # Representative rows, all TFs: actual fixed official implementation parity.
   ix=np.unique([0,len(Y)//2,len(Y)-1]);official,_=space['ulm'](Y[ix],wc.toarray(),verbose=False);assert np.allclose(official,actual[ix][:,eligible],rtol=2e-8,atol=2e-8);official_checked+=int(official.size)
  all_values+=checked;official_values+=official_checked;reports.append(dict(gse_id=sid,all_score_values_checked=checked,actual_official_function_values_checked=official_checked,status='PASS_ALL_OLS_AND_OFFICIAL_ULM_MATH'));print(json.dumps(reports[-1]),flush=True)
 contexts=pd.read_csv(O/'all_corrected_context_readiness.tsv',sep='\t',keep_default_na=False);tfrows=pd.read_csv(O/'private_all_context_TF_readiness.tsv.gz',sep='\t',keep_default_na=False);assert not contexts.context_id.duplicated().any() and not tfrows.duplicated(['context_id','TF']).any();groups=pd.read_csv(N/'all276_major_direction_contexts.tsv',sep='\t',keep_default_na=False).set_index('group_id');count=0;qualified=0
 for cid,sub in tfrows.groupby('context_id',sort=False):
  row=contexts.set_index('context_id').loc[cid];g=groups.loc[row.group_id];a=np.load(N/'private_group_arrays'/(row.group_id+'.npz'));out=np.load(O/'private_context_arrays'/(cid+'.npz'));mask=a['corrected_material_for_next_inference']==row.material;assert mask.any();assert out['samples'].tolist()==a['samples'][mask].tolist();bb=a['corrected_blocks_for_next_inference'][mask];rv=a['corrected_identity_reviewed_for_next_inference'][mask];assert np.array_equal(bb,out['blocks']) and np.array_equal(rv,out['identity_reviewed']);n=int(mask.sum());b=len(set(bb));passed=n>=20 and b>=10;assert int(row.n_libraries)==n and int(row.n_analysis_blocks)==b and bool(row.coverage_ge20libs_ge10blocks)==passed
  units,src=source_data[g.gse_id];lookup={(x.sampleID,x.majorCluster):i for i,x in enumerate(units.itertuples(index=False))};ti=[lookup[(s,g.target)] for s in out['samples']]
  for prefix in ['primary','sensitivity']:assert np.array_equal(out[prefix+'_TF_scores'],src[prefix+'_scores'][ti],equal_nan=True)
  support=(src['primary_detected_target_counts'][ti]>=10).sum(0);vary=np.array([np.ptp(v[np.isfinite(v)])>1e-10 if np.isfinite(v).sum()>1 else False for v in out['primary_TF_scores'].T]);ready=(src['primary_target_counts']>=20)&(support>=.8*n)&vary&passed;assert np.array_equal(ready,out['TF_ready']);s=sub.set_index('TF').loc[TFs];assert np.array_equal(s.measured_targets.to_numpy(int),src['primary_target_counts']);assert np.array_equal(s.libraries_detected_at_least10targets.to_numpy(int),support);assert np.array_equal(s.TF_ready.to_numpy(bool),ready);assert int(row.TFs_primary_ready)==int(ready.sum())
  l=a['ligand_meanlog'][mask];lp=a['ligand_min_positive_fraction'][mask];rp=a['receptor_min_positive_fraction'][mask];lrs=a['native_measured']&(np.ptp(np.nan_to_num(l),axis=0)>1e-12)&((lp>=.01).sum(0)>=.5*n)&((rp>=.01).sum(0)>=.8*n)&passed;assert np.array_equal(lrs,out['LR_RNA_support_ready']) and int(row.LR_RNA_support_ready)==int(lrs.sum());count+=1;qualified+=int(passed)
 assert count==inv['corrected_contexts'] and len(tfrows)==count*1163;assert set(contexts.group_id)==set(groups.index);assert qualified==inv['qualified_contexts']
 result=dict(status='PASS_ALL_TF_SCORES_OLS_OFFICIAL_NUMERIC_FUNCTIONS_AND_CORRECTED_CONTEXTS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),all_score_values_checked=all_values,actual_official_function_values_checked=official_values,corrected_contexts_checked=count,context_TF_records_checked=len(tfrows),qualified_contexts_checked=qualified,source_reports=reports,no_response_p_or_mechanism_validation=True);(O/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
