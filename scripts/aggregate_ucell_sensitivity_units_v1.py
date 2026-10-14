"""Preserve all four score variants and QC per frozen sample-node, after full audit."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.sparse import csr_matrix

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);A=V/'all_source_native_UCell_audit_aggregation_20261014_v2';ai=json.loads((A/'inventory.json').read_text());assert ai['status']=='PASS_ALL4532933_IDENTITIES_SUPPORT_MASKS_AND_SAMPLED_MANUAL_UCELL_WITH_PRIVATE_UNIT_MEANS';S=V/'all_source_native_UCell_projection_20261014_v2';iv=json.loads((S/'inventory.json').read_text());assert ai['producer_inventory_sha256']==sha(S/'inventory.json');D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert ai['dictionary_sha256']==sha(D);recs=json.loads(D.read_text())['programmes'];F=R/'all_CM_state_native_input_preflight_20261014_v1';fi=json.loads((F/'inventory.json').read_text());assert sha(F/'private_all_frozen_sample_node_units.tsv.gz')==fi['private_unit_table_sha256'];units=pd.read_csv(F/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False)
 variants=['scores_maxrank1000','scores_maxrank1500','scores_maxrank3000','TPM_signature_scores_maxrank1500'];aggs={};seen=np.zeros(len(units),np.int64)
 for major in sorted(units.major.unique()):
  ui=units.loc[units.major.eq(major),'unit_index'].to_numpy(int);ids=np.array([r['programme_id'] for r in recs if major in r['allowed_majors']]);shape=(len(ui),len(ids));aggs[major]=dict(ui=ui,ids=ids,sums={k:np.zeros(shape) for k in variants},valid={k:np.zeros(shape,np.int64) for k in variants},RNA=np.zeros(len(ui)),detected=np.zeros(len(ui)))
 O=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';assert not O.exists();O.mkdir(mode=0o700);done=0
 for file in iv['private_score_files']:
  p=S/file['relative_file'];assert sha(p)==file['sha256'];z=dict(np.load(p,allow_pickle=False));uid=z['unit_index'];a=aggs[file['major']];loc=np.searchsorted(a['ui'],uid);assert np.array_equal(a['ui'][loc],uid) and np.array_equal(a['ids'],z['programme_ids']);assert units.iloc[uid].major.eq(file['major']).all() and units.iloc[uid].gse_id.eq(file['gse_id']).all();n=len(uid);inc=csr_matrix((np.ones(n),(loc,np.arange(n))),shape=(len(a['ui']),n));seen+=np.bincount(uid,minlength=len(units))
  for k in variants:
   y=z[k].astype(float);mask=np.isfinite(y);v=np.nan_to_num(y);summed=inc@v;check=np.zeros_like(summed);np.add.at(check,loc,v);assert np.allclose(summed,check,atol=1e-10,rtol=1e-12);a['sums'][k]+=summed;a['valid'][k]+=np.asarray(inc@mask.astype(np.int64),np.int64)
  for key,field in [('RNA','full_RNA_total'),('detected','full_detected_genes')]:
   summed=np.asarray(inc@z[field]).ravel();check=np.zeros_like(summed);np.add.at(check,loc,z[field]);assert np.array_equal(summed,check);a[key]+=summed
  done+=n;(O/'progress.json').write_text(json.dumps(dict(cells=done,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
 assert done==4532933 and np.array_equal(seen,units.n_cells.to_numpy());outputs=[]
 for i,(major,a) in enumerate(aggs.items()):
  old=next(x for x in ai['private_aggregate_files'] if x['major']==major);op=A/old['relative_file'];assert sha(op)==old['sha256'];prior=dict(np.load(op,allow_pickle=False));assert np.array_equal(prior['unit_index'],a['ui']) and np.array_equal(prior['programme_ids'],a['ids']);out=dict(unit_index=a['ui'],programme_ids=a['ids'],cells=seen[a['ui']],mean_full_RNA=a['RNA']/seen[a['ui']],mean_full_detected_genes=a['detected']/seen[a['ui']])
  for k in variants:out[k+'_mean']=np.divide(a['sums'][k],a['valid'][k],out=np.full_like(a['sums'][k],np.nan),where=a['valid'][k]>0);out[k+'_valid_cells']=a['valid'][k]
  assert np.array_equal(out['scores_maxrank1500_valid_cells'],prior['valid_score_cells']) and np.allclose(out['scores_maxrank1500_mean'],prior['programme_mean'],atol=1e-10,rtol=1e-12,equal_nan=True);p=O/f'private_major_{i:02}_sensitivity_units.npz';np.savez_compressed(p,**out);outputs.append(dict(major=major,relative_file=p.name,sha256=sha(p),units=len(a['ui'])))
 out=dict(status='COMPLETE_ALL4532933_UNIT_SCORE_SENSITIVITIES_QC_WITH_PRIMARY_AUDIT_AGREEMENT',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),score_inventory_sha256=sha(S/'inventory.json'),audit_inventory_sha256=sha(A/'inventory.json'),dictionary_sha256=sha(D),cells=done,units=len(units),variants=variants,all_chunks_two_aggregation_methods=True,primary1500_agrees_with_independent_audit=True,private_outputs=outputs,formal_states=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='private_outputs'}))
if __name__=='__main__':main()
