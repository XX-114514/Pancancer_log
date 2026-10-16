"""Independent all-row cache/partition recount, native gene-axis and QR design audit."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import qr as pivoted_qr
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();O=V/'state_after_DIALOGUE_all_LR_preflight_20261016_v1';iv=json.loads((O/'inventory.json').read_text())
for n,h in iv['files'].items():assert sha(O/n)==h
T=pd.read_csv(O/'all129472_LR_partition_coverage.tsv',sep='\t',keep_default_na=False);assert len(T)==129472 and not T.duplicated(['context','source','target','ligand','receptor']).any();A=pd.read_csv(O/'all28_cache_alignment_checks.tsv',sep='\t');E=R/'immune_same_name_family_source_review_v1/CM_edge_LR_extraction_v2';ei=json.loads((E/'inventory.json').read_text());assert sha(E/'inventory.json')==iv['LR_inventory_sha256'];F=R/'all_CM_state_native_input_preflight_20261014_v1';fs=pd.read_csv(F/'all_source_native_scope.tsv',sep='\t').set_index('gse_id');Q=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';qi=json.loads((Q/'inventory.json').read_text());qc={}
for info in qi['private_outputs']:
 p=Q/info['relative_file'];assert sha(p)==info['sha256'];z=np.load(p,allow_pickle=False)
 for k,r,d in zip(z['unit_index'],z['mean_full_RNA'],z['mean_full_detected_genes']):qc[int(k)]=(float(r),float(d))
checked=0;native_checks=0;rankchecks=0;pattern_cache={}
for a in A.itertuples():
 root=V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{a.context:04d}';ci=json.loads((root/'private_input_inventory.json').read_text());parts=pd.read_csv(root/'private_block_library_partition.tsv',sep='\t');p=E/'private_group_arrays'/(a.old_group+'.npz');assert sha(p)==ei['outputs'][str(p.relative_to(E))];z=np.load(p,allow_pickle=False);ri=pd.Index(z['samples'].astype(str)).get_indexer(parts.sample_uid.astype(str));assert (ri>=0).all();assert list(z['patient_blocks'][ri])==list(parts.effective_resampling_block);reverse=a.source==ci['endpoints'][1];assert a.source==ci['endpoints'][int(reverse)] and a.target==ci['endpoints'][int(not reverse)];unitS=parts['unit_index_B' if reverse else 'unit_index_A'];unitT=parts['unit_index_A' if reverse else 'unit_index_B'];qs=np.array([qc[int(k)] for k in unitS]);qt=np.array([qc[int(k)] for k in unitT]);X=np.column_stack([np.log1p(z['source_n'][ri]),np.log1p(z['target_n'][ri]),z['source_within_major'][ri],z['target_within_major'][ri],np.log1p(qs[:,0]),np.log1p(qt[:,0]),np.log1p(qs[:,1]),np.log1p(qt[:,1])]);train=parts.partition.eq('discovery').to_numpy();test=parts.partition.eq('holdout').to_numpy();values=z['scores'][ri,:,0];gp=AR/'20261010_uniform_expression_annotation_v1/by_gse'/a.gse_id/'genes.tsv';assert sha(gp)==fs.loc[a.gse_id,'native_gene_axis_sha256'];ng=set(pd.read_csv(gp,sep='\t').gene);native=np.array([set(str(l).split('_')+str(r).split('_'))<=ng for l,r in zip(z['ligands'],z['receptors'])]);assert np.array_equal(native,z['native_measured']);native_checks+=len(native);block=T[(T.context==a.context)&T.source.eq(a.source)&T.target.eq(a.target)].set_index(['ligand','receptor']);assert len(block)==4624
 for j,lr in enumerate(zip(z['ligands'],z['receptors'])):
  row=block.loc[tuple(lr)];good=np.isfinite(values[:,j]);it=np.flatnonzero(good&train);ih=np.flatnonzero(good&test);assert int(row.observed_discovery_blocks)==len(it) and int(row.observed_holdout_blocks)==len(ih) and bool(row.native_all_subunits_measured)==native[j];Y=X[it];rank=0;const=False
  if len(Y):
   sd=Y.std(0);const=bool((sd<=1e-12).any());Z=np.column_stack([np.ones(len(Y)),(Y-Y.mean(0))/np.where(sd>1e-12,sd,1)]);_,qr,_=pivoted_qr(Z,mode='economic',pivoting=True);rank=int((np.abs(np.diag(qr))>max(Z.shape)*np.finfo(float).eps*np.linalg.norm(Z,2)).sum())
  assert rank==int(row.discovery_design_rank),(a.context,lr,rank,row.discovery_design_rank);assert const==bool(row.constant_covariate);rankchecks+=1
  vt=len(it)>0 and np.ptp(values[it,j])>0;vh=len(ih)>0 and np.ptp(values[ih,j])>0;assert vt==bool(row.discovery_LR_variable) and vh==bool(row.holdout_LR_variable)
  expected='INPUT_ELIGIBLE_EFFECT_AND_MASKED_RECEIVER_SCORE_NOT_TESTED'
  # Reverse priority assignment independently reproduces the producer's first exclusion.
  if not vh:expected='CONSTANT_HELDOUT_LR'
  if len(ih)<10:expected='FEWER_THAN_10_OBSERVED_INTERNAL_HOLDOUT_BLOCKS'
  if not vt:expected='CONSTANT_DISCOVERY_LR'
  if const or rank!=9:expected='FULL_EIGHT_COVARIATE_DISCOVERY_DESIGN_NOT_IDENTIFIABLE'
  if len(it)/int(train.sum())<.8:expected='DISCOVERY_OUTPUT_COVERAGE_BELOW_80PCT'
  if len(it)<20:expected='FEWER_THAN_20_OBSERVED_DISCOVERY_BLOCKS'
  if not native[j]:expected='NATIVE_COMPLEX_SUBUNIT_UNMEASURED'
  assert row.status==expected,(a.context,lr,row.status,expected);checked+=1
assert checked==129472 and native_checks==129472;out=dict(status='PASS_ALL129472_CACHE_PARTITION_NATIVE_GENE_AND_QR_DESIGN_RECORDS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),records_checked=checked,native_complex_flags_checked=native_checks,independent_QR_design_checks=rankchecks,association_tests_performed=False,all_cases_reported_missing_not_zero=True);p=O/'independent_audit.json';assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
