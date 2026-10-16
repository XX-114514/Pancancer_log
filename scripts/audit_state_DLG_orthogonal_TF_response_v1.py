"""Independent official-ULM readback, scalar SVD effects and BH for frozen TF family."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,ast
import numpy as np,pandas as pd
from scipy.stats import rankdata,t
from scipy.linalg import lstsq
from scipy.sparse import csr_matrix
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);D=V/'state_DLG_orthogonal_TF_response_20261016_v1'
for tick in range(480):
 if (D/'inventory.json').exists():break
 e=V/'state_DLG_orthogonal_TF_response_20261016_v1.execution.json'
 if e.exists() and json.loads(e.read_text()).get('exit_code') is not None:raise RuntimeError('Producer terminated without result inventory')
 time.sleep(30)
else:raise TimeoutError('Bounded four-hour dependency wait expired')
iv=json.loads((D/'inventory.json').read_text())
for f,h in iv['files'].items():assert sha(D/f)==h
F=V/'state_DLG_orthogonal_TF_preflight_20261016_v1';fi=json.loads((F/'inventory.json').read_text());S=pd.read_csv(F/'private_fixed_TF_slots.tsv',sep='\t');keys=['context','source','target','MCP','ligand','receptor','TF'];fixed=set(S[keys].itertuples(index=False,name=None));assert len(fixed)==5448
T=pd.read_csv(D/'private_all_TF_association_results.tsv',sep='\t');assert not T.duplicated(keys).any();saved={tuple(r[k] for k in keys):r for r in T.to_dict('records')};scope=pd.read_csv(D/'public_all5448_response_scope.tsv',sep='\t');assert scope.records.sum()==5448
N=len(T);ix=np.argsort(T.permutation_p.to_numpy(),kind='stable');bh=np.empty(N);bh[ix]=np.minimum(1,np.minimum.accumulate((T.permutation_p.to_numpy()[ix]*N/np.arange(1,N+1))[::-1])[::-1]);assert np.allclose(bh,T.BH_q,atol=1e-12);assert T.loc[T.valid_bootstraps<1900,['CI_low','CI_high']].isna().all().all();assert T.loc[T.valid_bootstraps>=1900,['CI_low','CI_high']].notna().all().all()
C=R/'fullgene_receiver_target_code_v1/collectri_reference_v1';netp=C/'human_collectri_signed_complexes.tsv';off=C/'official_decoupler192_ULM_source.txt';assert sha(netp)==fi['prior_TF_network_sha256'] and sha(off)==fi['official_ULM_source_sha256'];net=pd.read_csv(netp,sep='\t');nets={k:dict(zip(g.target,g.weight)) for k,g in net.groupby('source')};fs=[n for n in ast.parse(off.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in ['mat_cov','mat_cor','t_val','ulm']];assert len(fs)==4;ns=dict(np=np,csr_matrix=csr_matrix,t=t);exec(compile(ast.Module(body=fs,type_ignores=[]),str(off),'exec'),ns);gweight=pd.read_csv(V/'state_after_DIALOGUE_frozen_projection_20261016_v1/public_fixed_gene_weights.tsv',sep='\t');A=V/'state_allgene_block_components_20261016_v1';ai=json.loads((A/'inventory.json').read_text());O=V/'state_DLG_orthogonal_TF_independent_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);seen=set();checked=set();official_errors=[];point_errors=[];held_errors=[];public=[]
def rho(x,y,D=None):
 x=rankdata(x).astype(float);y=rankdata(y).astype(float);x-=x.mean();y-=y.mean();nx=np.linalg.norm(x);ny=np.linalg.norm(y)
 if D is not None:x-=D@lstsq(D,x,lapack_driver='gelsd')[0];y-=D@lstsq(D,y,lapack_driver='gelsd')[0]
 if np.linalg.norm(x)<=1e-12*max(nx,1) or np.linalg.norm(y)<=1e-12*max(ny,1):return np.nan
 return float(x@y/(np.linalg.norm(x)*np.linalg.norm(y)))
for p in sorted((D/'private_case_arrays').glob('group_*.npz')):
 z=np.load(p,allow_pickle=False);cases=[json.loads(v) for v in z['case_json']];r0=cases[0];ctx=r0['context'];ci=json.loads((V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{ctx:04d}'/'private_input_inventory.json').read_text());typ='type'+['A','B'][ci['endpoints'].index(r0['target'])];gpath=A/(str(ctx)+'_'+typ)/'private_genewise_block_components.tsv.gz';assert sha(gpath)==ai['files'][str(gpath.relative_to(A))];B=pd.read_csv(gpath,sep='\t',index_col=0).T.reindex(z['sample_ids']);mcp=set(gweight.loc[(gweight.context==ctx)&gweight.type.eq(typ)&gweight.program.eq(r0['MCP']),'gene']);maxerr=0;ec=0
 for col,case in enumerate(cases):
  key=tuple(case[k] for k in keys);assert key in fixed and key not in seen;seen.add(key);exclude=mcp|set(case['ligand'].split('_')+case['receptor'].split('_'));bg=[g for g in B.columns if g not in exclude];w=np.array([nets[case['TF']].get(g,0) for g in bg]);assert len(bg)==case['background_genes']>=5000 and np.count_nonzero(w)==case['outside_targets']>=10
  actual,_=ns['ulm'](B[bg].to_numpy(float),w[:,None]);activity=z['TF_activity'][:,col];assert np.allclose(actual[:,0],activity,atol=1e-7,equal_nan=True);err=float(np.nanmax(abs(actual[:,0]-activity))) if np.isfinite(activity).any() else 0;official_errors.append(err);maxerr=max(maxerr,err);ec+=1
  if key not in saved:continue
  row=saved[key];checked.add(key);x=z['log_expr_prod'][:,col];mask=np.isfinite(x);n=int(mask.sum());assert n==row['blocks'] and n>=20 and n/len(x)>=.8;cov=z['covariates'][mask];sd=cov.std(0);assert (sd>1e-12).all();design=np.column_stack([np.ones(n),(cov-cov.mean(0))/sd]);assert np.linalg.matrix_rank(design)==9 and n-10>=10;y=activity[mask];x=x[mask];mc=z['LR_masked_MCP'][mask,col];calc=[rho(x,y),rho(x,y,design),rho(mc,y,design)];ref=[row['raw_LR_TF_rho'],row['adjusted_LR_TF_rho'],row['MCP_TF_adjusted_rho']];assert np.allclose(calc,ref,atol=1e-8,equal_nan=True);point_errors.extend(abs(a-b) for a,b in zip(calc,ref) if np.isfinite(a) and np.isfinite(b));tr=z['partitions'][mask]=='discovery';te=z['partitions'][mask]=='holdout'
  if tr.sum()>=20 and te.sum()>=10 and np.linalg.matrix_rank(design[tr])==9:
   rx=x[te]-design[te]@lstsq(design[tr],x[tr],lapack_driver='gelsd')[0];ry=y[te]-design[te]@lstsq(design[tr],y[tr],lapack_driver='gelsd')[0];h=rho(rx,ry);assert np.isclose(h,row['internal_holdout_rho'],atol=1e-8,equal_nan=True)
   if np.isfinite(h):held_errors.append(abs(h-row['internal_holdout_rho']))
  else:assert pd.isna(row['internal_holdout_rho'])
 public.append(dict(context=ctx,source=r0['source'],target=r0['target'],MCP=r0['MCP'],all_fixed_activity_checks=ec,max_official_activity_error=maxerr));(O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),activities_checked=len(seen),point_tests_checked=len(checked)))+'\n')
assert seen==fixed and checked==set(saved);pd.DataFrame(public).to_csv(O/'public_allgroup_activity_audit.tsv',sep='\t',index=False);result=dict(status='PASS_ALL5448_OFFICIAL_ULM_ACTIVITIES_ALL_ESTIMABLE_SVD_POINT_HELDOUT_BH_AND_CI_GATES',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_inventory_sha256=sha(D/'inventory.json'),auditor_sha256=sha(__file__),activities_checked=len(seen),estimable_point_tests_checked=len(checked),max_official_ULM_error=max(official_errors),max_point_error=max(point_errors),max_heldout_error=max(held_errors) if held_errors else None,BH_records=int((T.BH_q<.05).sum()),all_primary_permutations_or_bootstrap_endpoints_independently_replayed=False,formal_states_or_causal_validation=False,files={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
