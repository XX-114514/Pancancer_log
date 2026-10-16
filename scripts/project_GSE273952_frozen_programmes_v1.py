"""Frozen UCell programmes and fixed signed target RNA readouts on mixed spatial spots."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,shutil,gc,importlib.metadata
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
import pyucell as uc

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);planp=Path('/tmp/GSE273952_frozen_programme_projection_plan_v1.json');plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'];assert importlib.metadata.version('pyucell')=='0.8.0';A=V/'existing_GSE273952_Visium_input_audit_20261016_v1'
for tick in range(120):
 if (A/'inventory.json').exists():break
 e=V/'existing_GSE273952_Visium_input_audit_20261016_v1.execution.json'
 if e.exists() and json.loads(e.read_text()).get('exit_code') is not None:raise RuntimeError('Visium input audit terminal without acceptance')
 time.sleep(30)
else:raise TimeoutError('One-hour input-audit wait')
ai=json.loads((A/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL71398');I=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';iv=json.loads((I/'inventory.json').read_text());assert ai['producer_inventory_sha256']==sha(I/'inventory.json');schema=json.loads((I/'private_native_gene_and_library_schema.json').read_text());lookup=schema['feature_lookup'];dictionary=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(dictionary)==plan['dictionary_sha256'];records=json.loads(dictionary.read_text())['programmes'];ids=[r['programme_id'] for r in records];assert len(ids)==814;Q=np.load(A/'private_all_spot_derived_QC_and_grid.npz',allow_pickle=False);assert sha(A/'private_all_spot_derived_QC_and_grid.npz')==ai['private_QC_sha256'];qc=Q['RNA_total'];n=len(qc);assert n==71398
old=AR.parent.parent/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';path=Path(json.loads((old/'outputs/spatial/GSE273952_spatial_contract_audit.json').read_text())['input_h5ad']);assert sha(path)==iv['input_file_sha256'];T=V/'joint_NAMPT_DDIT3_spatial_measurement_20261016_v1/private_frozen_target_and_background_definition.json';assert sha(T)==plan['frozen_DDIT3_definition_sha256'];tf=json.loads(T.read_text());weights=tf['outside_target_weights'];assert len(weights)==74 and set(weights.values())<=set([-1,1]);assert all(g in lookup for g in weights);tfpos=sorted(g for g,w in weights.items() if w==1);tfneg=sorted(g for g,w in weights.items() if w==-1);extra=['RNA_DDIT3_TARGETS_POS','RNA_DDIT3_TARGETS_NEG','RNA_DDIT3_TARGETS_SIGNED'];allids=ids+extra;signatures={};indices=[];scorable=[];coverage=[]
for r in records:
 present=[g for g in r['primary_genes'] if g in lookup];jj=[lookup[g] for g in present];unique=list(dict.fromkeys(jj));gate=len(unique)>=10 and len(unique)/len(r['primary_genes'])>=.8;indices.append(np.array(unique,int));scorable.append(gate)
 if gate:signatures[r['programme_id']]=['f'+str(j) for j in unique]
 coverage.append(dict(programme_id=r['programme_id'],role='frozen_single_cell_candidate_programme_on_mixed_spots',signature_genes=len(r['primary_genes']),measured_physical_genes=len(unique),measurement_gate=gate,allowed_single_cell_majors=';'.join(r['allowed_majors']),spot_cell_identity_assigned=False,formal_state=False))
for key,genes in zip(extra,[tfpos,tfneg,sorted(weights)]):
 jj=np.array([lookup[g] for g in genes],int);indices.append(jj);scorable.append(True)
 if key==extra[-1]:signatures[key]=['f'+str(lookup[g])+('+' if w==1 else '-') for g,w in weights.items()]
 elif len(jj):signatures[key]=['f'+str(j) for j in jj]
 coverage.append(dict(programme_id=key,role='fixed_prior_target_RNA_proxy_not_original_ULM_or_protein',signature_genes=len(genes),measured_physical_genes=len(genes),measurement_gate=True,allowed_single_cell_majors='not_cell_intrinsic_in_mixed_spots',spot_cell_identity_assigned=False,formal_state=False))
# Fixed signed-UCell fixture: explicitly validate clipping/zero/ties before any tissue score.
fixture=np.array([[5,5,1,0,0,0,0,0],[0,0,8,8,0,0,0,0],[0,0,0,0,0,0,0,0],[5,3,7,2,1,0,0,0]],float);test=ad.AnnData(csr_matrix(fixture),var=pd.DataFrame(index=['f'+str(i) for i in range(8)]));uc.compute_ucell_scores(test,signatures={'P':['f0','f1'],'N':['f2','f3'],'S':['f0+','f1+','f2-','f3-']},max_rank=6,ties_method='average',missing_genes='skip',w_neg=1,n_jobs=1,device='cpu');expected=[]
for row in fixture:
 ranks=np.full(8,6,float);nz=np.flatnonzero(row>0);ranks[nz]=np.minimum(6,rankdata(-row[nz],method='average'));p=1-(ranks[:2].sum()-3)/9;g=1-(ranks[2:4].sum()-3)/9;expected.append([p,g,max(p-g,0)])
assert np.allclose(test.obs[['P_UCell','N_UCell','S_UCell']],expected,atol=1e-7);O=V/'GSE273952_frozen_programme_projection_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);(O/'private_spot_scores').mkdir(mode=0o700);full=np.full((n,len(allids)),np.nan,np.float32);files=[];genecount=None
with h5py.File(path,'r') as f:
 X=f['layers/counts'];genecount=int(X.attrs['shape'][1]);features=pd.DataFrame(index=['f'+str(i) for i in range(genecount)])
 for start in range(0,n,plan['spots_per_chunk']):
  assert shutil.disk_usage(V).free>(500+plan['output_budget_GiB'])*2**30;end=min(n,start+plan['spots_per_chunk']);ip=X['indptr'][start:end+1].astype(np.int64);x=csr_matrix((X['data'][int(ip[0]):int(ip[-1])],X['indices'][int(ip[0]):int(ip[-1])],ip-ip[0]),shape=(end-start,genecount));x.sum_duplicates();x.sort_indices();assert np.array_equal(np.asarray(x.sum(1)).ravel(),qc[start:end]);a=ad.AnnData(x,var=features.copy());uc.compute_ucell_scores(a,signatures=signatures,max_rank=1500,ties_method='average',missing_genes='skip',w_neg=1.0,chunk_size=plan['ucell_rank_batch'],n_jobs=1,device='cpu',suffix='_native');scores=np.full((end-start,len(allids)),np.nan,np.float32);support=np.zeros(scores.shape,np.uint8);frac=np.zeros(scores.shape,np.float32)
  for j,key in enumerate(allids):
   if scorable[j]:
    if key in signatures:scores[:,j]=a.obs[key+'_native'].to_numpy()
    else:assert j>=814 and len(indices[j])==0;scores[:,j]=0.0
   jj=indices[j]
   if len(jj):support[:,j]=np.asarray((x[:,jj]>0).sum(1)).ravel();frac[:,j]=np.asarray(x[:,jj].sum(1)).ravel()/qc[start:end]
  assert np.allclose(scores[:,-1],np.maximum(scores[:,-3]-scores[:,-2],0),atol=1e-7);full[start:end]=scores;rawlr=np.column_stack([np.asarray(x[:,lookup[g]].toarray()).ravel() for g in ['NAMPT','ITGA5','ITGB1','DDIT3']]);p=O/'private_spot_scores'/f'chunk_{start:08d}.npz';np.savez_compressed(p,spot_rows=np.arange(start,end),library_codes=Q['library_codes'][start:end],programme_ids=np.array(allids),scores=scores,detected_signature_genes=support,signature_RNA_fraction=frac,raw_NAMPT_ITGA5_ITGB1_DDIT3_counts=rawlr,RNA_total=qc[start:end],detected_genes=Q['detected_genes'][start:end]);files.append(dict(relative_file=str(p.relative_to(O)),sha256=sha(p),spots=end-start));(O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),spots_completed=end,total_spots=n,score_columns=len(allids)))+'\n');del a,x,scores,support,frac;gc.collect()
pd.DataFrame(coverage).to_csv(O/'public_all817_programme_measurement_scope.tsv',sep='\t',index=False);summary=[]
for j,key in enumerate(allids):
 val=full[:,j];ok=np.isfinite(val);v=val[ok];summary.append(dict(programme_id=key,valid_spots=int(ok.sum()),mean_score=float(v.astype(float).mean()) if len(v) else np.nan,q10=float(np.quantile(v,.1)) if len(v) else np.nan,median_score=float(np.median(v)) if len(v) else np.nan,q90=float(np.quantile(v,.9)) if len(v) else np.nan,interpretation='mixed_spot_RNA_programme_not_cell_identity_or_protein_activity'))
pd.DataFrame(summary).to_csv(O/'public_all817_cohort_RNA_score_summary.tsv',sep='\t',index=False);out=dict(status='ALL71398_MIXED_SPOTS_814_FROZEN_PROGRAMMES_PLUS3_FIXED_TF_RNA_COMPONENTS_SCORED_PENDING_INDEPENDENT_AUDIT',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),input_audit_sha256=sha(A/'inventory.json'),raw_input_sha256=iv['input_file_sha256'],dictionary_sha256=sha(dictionary),fixed_TF_gene_definition_sha256=sha(T),spots=n,native_genes=genecount,source_libraries=19,columns=817,original_programme_measurement_gates=int(sum(scorable[:814])),DDIT3_positive_targets=len(tfpos),DDIT3_negative_targets=len(tfneg),DDIT3_signed_UCell_weight_negative=1.0,DDIT3_signed_UCell_negative_scores_clipped_to_zero=True,synthetic_signed_scores_checked=12,case_level_clinical_outcome_values_used=False,old_deconvolution_results_used=False,spot_cell_identity_assigned=False,formal_state_or_spatial_biological_validation=False,private_outputs=files,public_outputs={p.name:sha(p) for p in O.glob('public_*.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['private_outputs','public_outputs']}))
