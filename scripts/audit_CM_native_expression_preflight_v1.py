"""Separate sparse per-unit row projection of native expression and complexity.
Audits representation evidence, not physical binding or induced responses.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import pandas as pd,numpy as np
from scipy.sparse import load_npz

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8388608),b''):h.update(block)
 return h.hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'CM_native_expression_preflight_v1';U=AR/'20261010_uniform_expression_annotation_v1';B=AR/'20261010_expression_annotation_downstream_v1';L=P/'CM_edge_LR_extraction_v2';planpath=Path(os.environ['PANCANCER_NATIVE_PREFLIGHT_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['auditor_sha256'];inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath) and not (O/'independent_audit.json').exists()
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 c=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t',keep_default_na=False);d=pd.read_csv(B/'downstream_review_v2/number_to_annotation_dictionary.tsv',sep='\t',keep_default_na=False);name=dict(zip(d.node_id,d.subCluster));groups=pd.read_csv(L/'all68_directed_group_summary.tsv',sep='\t');reports=[];total=0;feature_groups=0
 for sid,gg in groups.groupby('gse_id',sort=True):
  root=U/'by_gse'/sid;units=pd.read_csv(root/'units.tsv',sep='\t',keep_default_na=False);assert np.array_equal(units.unit_row,np.arange(len(units)));genes=pd.read_csv(root/'genes.tsv',sep='\t',keep_default_na=False).gene;gi={g:i for i,g in enumerate(genes)};keys={}
  for unit in units.itertuples():
   if unit.candidate_eligible:keys[unit.sample_uid,name[c.iloc[unit.node_index].node_id]]=unit.unit_row
  logs=load_npz(root/'sum_log1p_CP10K.npz');positive=load_npz(root/'positive_counts.npz');cp=load_npz(root/'sum_CP10K.npz');badp=badm=obs=0
  for group in gg.itertuples():
   b=dict(np.load(L/'private_group_arrays'/(group.group_id+'.npz'),allow_pickle=False));f=dict(np.load(O/'private_group_features'/(group.group_id+'.npz'),allow_pickle=False));assert np.array_equal(f['samples'],b['samples']) and np.array_equal(f['ligands'],b['ligands']) and np.array_equal(f['receptors'],b['receptors']);sr=[keys[s,group.source] for s in b['samples']];tr=[keys[s,group.target] for s in b['samples']];sn=units.iloc[sr].n_cells.to_numpy();tn=units.iloc[tr].n_cells.to_numpy();assert np.array_equal(sn,b['source_n']) and np.array_equal(tn,b['target_n']);assert np.allclose(np.asarray(cp[sr].sum(axis=1)).ravel()/sn,10000,rtol=2e-6,atol=.02) and np.allclose(np.asarray(cp[tr].sum(axis=1)).ravel()/tn,10000,rtol=2e-6,atol=.02)
   for rows,n,label in [(sr,sn,'source'),(tr,tn,'target')]:
    assert np.allclose(np.asarray(positive[rows].sum(axis=1)).ravel()/n,f[label+'_mean_detected_genes'],atol=1e-10);assert np.allclose(np.asarray(logs[rows].sum(axis=1)).ravel()/n,f[label+'_mean_sum_log1p_CP10K'],rtol=1e-10,atol=1e-10)
   # Independent small sparse column caches, then explicit per-complex minima.
   sx={};tx={}
   for ligand,receptor in zip(b['ligands'],b['receptors']):
    for gene in set(ligand.split('_')+receptor.split('_')):
     if gene in gi and gene not in sx:sx[gene]=logs[sr,gi[gene]].toarray().ravel()/sn;tx[gene]=logs[tr,gi[gene]].toarray().ravel()/tn
   for j,(ligand,receptor) in enumerate(zip(b['ligands'],b['receptors'])):
    measured=all(g in sx for g in ligand.split('_')+receptor.split('_'));assert measured==bool(f['native_measured'][j])
    if not measured:assert np.isnan(f['native_meanlog_exprprod_and_lrmeans'][:,j]).all();continue
    sl=np.minimum.reduce([sx[g] for g in ligand.split('_')]);rr=np.minimum.reduce([tx[g] for g in receptor.split('_')]);native=np.column_stack([sl*rr,(sl+rr)/2]);assert np.allclose(native,f['native_meanlog_exprprod_and_lrmeans'][:,j],rtol=1e-10,atol=1e-10);mask=np.isfinite(b['scores'][:,j,0]);obs+=int(mask.sum());badp+=int((~np.isclose(native[mask,0],b['scores'][mask,j,0],rtol=2e-5,atol=2e-7)).sum());badm+=int((~np.isclose(native[mask,1],b['scores'][mask,j,1],rtol=2e-5,atol=2e-7)).sum())
   feature_groups+=1
  expected=next(row for row in inv['summary'] if row['gse_id']==sid);assert obs==expected['observed_original_LR_records'] and badp==expected['product_discrepancies'] and badm==expected['lrmeans_discrepancies'];total+=obs;reports.append(dict(gse_id=sid,observed_scores_checked=obs,exprprod_discrepancies=badp,lrmeans_discrepancies=badm,status='PASS_INDEPENDENT_NATIVE_PROJECTION'));print(json.dumps(reports[-1]),flush=True)
 assert total==381170 and feature_groups==68;report=dict(status='PASS_ALL68_NATIVE_GROUP_FEATURES_AND381170_SCORE_REPRESENTATION_COMPARISONS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),native_feature_groups_checked=68,original_scores_checked=total,source_reports=reports,biological_validation=False,fresh_raw_cell_recount_this_audit=False);(O/'independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
