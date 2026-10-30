"""Independent whole-row readback of the lineage model-mask adapter."""
import argparse,datetime,hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def main(R):
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);end=P/'provenance/supervisor_terminal_v1.json'
 while not end.exists():time.sleep(5)
 assert load(end)['actual_child_exit_code']==0;pc=load(P/'config/spec_v1.json');H=Path(pc['parent_producer']);G=Path(pc['parent_graph']);t=load(P/'provenance/terminal_v1.json');bridge=P/'outputs/accepted_lineage_model_input_bridge_v1.json';assert sha(bridge)==t['bridge_SHA256'];b=load(bridge);assert sha(b['producer_manifest']['path'])==b['producer_manifest']['SHA256']==t['manifest_SHA256'];m=load(b['producer_manifest']['path']);original=load(b['original_identity_acceptance_reference']['path']);oldmanifest=load(original['producer_manifest']['path']);oldbase=Path(oldmanifest['index']['path']).parent;old={r['source']:r for r in load(oldmanifest['index']['path'])};assert sha(m['index']['path'])==m['index']['SHA256'];index=load(m['index']['path']);sources=load(H/'outputs/scope_and_common_genes_v1.json')['sources'];assert [r['source'] for r in index]==sources
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 labels=pd.read_csv(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz',sep='\t',usecols=['leiden_0_5']).iloc[:,0].to_numpy(np.int32);allowed=set(pc['core_parent_clusters'])|set(pc['whole_parent_control_clusters']);expected_positions=np.asarray([i for i,k in enumerate(labels) if int(k) in allowed],np.int64);np.testing.assert_array_equal(np.load(m['parent_fit_positions']['path']),expected_positions);assert sha(m['parent_fit_positions']['path'])==m['parent_fit_positions']['SHA256'];alln=0;records=[]
 for code,r in enumerate(index):
  source=r['source'];o=old[source];assert r['metadata']==str(oldbase/o['metadata']) and r['metadata_SHA256']==o['metadata_SHA256'];p=P/'outputs'/r['arrays'];assert sha(p)==r['arrays_SHA256'];assert sha(oldbase/o['arrays'])==o['arrays_SHA256']
  with np.load(oldbase/o['arrays'],allow_pickle=False) as f:original_arrays={k:f[k] for k in ['native_row','source_native_position','fit_pool_candidate_not_final_QC','source_qualification_code','native_total_RNA']}
  with np.load(p,allow_pickle=False) as f:new={k:f[k] for k in f.files}
  for field in ['native_row','source_native_position','source_qualification_code','native_total_RNA']:np.testing.assert_array_equal(new[field],original_arrays[field])
  np.testing.assert_array_equal(new['parent_fit_pool_candidate_not_final_QC'],original_arrays['fit_pool_candidate_not_final_QC']);selected=expected_positions[axes['source_code'][expected_positions]==code];expected_rows=axes['native_row'][selected];actual_rows=new['native_row'][new['fit_pool_candidate_not_final_QC']];np.testing.assert_array_equal(actual_rows,expected_rows);assert np.all(~new['fit_pool_candidate_not_final_QC']|original_arrays['fit_pool_candidate_not_final_QC']);np.testing.assert_array_equal(new['source_native_position'][new['fit_pool_candidate_not_final_QC']],axes['source_native_position'][selected]);assert len(actual_rows)==pc['expected_per_source'][source] and len(actual_rows)>100;alln+=len(actual_rows);records.append(dict(source=source,n_selected_fit=len(actual_rows),all_original_qualified_axes_and_parent_QC_mask_preserved=True))
 assert alln==pc['expected_selected_n']==t['n_selected_fit'];save(R/'outputs/accepted_lineage_model_input_scope_v1.json',dict(n_selected_fit=alln,n_sources=len(index),all_source_masks_independently_rebuilt_from_parent_graph_members=True,all_original_qualified_row_axes_RNA_and_parent_fit_masks_unchanged=True,whole_parents_and_entire_control_parent_included=True,bridge_path=str(bridge),bridge_SHA256=sha(bridge),parent_manifest_SHA256=sha(b['original_identity_acceptance_reference']['path']),source_checks=records,not_new_identity_or_QC_or_ABC_winner=True));save(R/'provenance/terminal_v1.json',dict(status='FULL_PARENT_LINEAGE_MODEL_INPUT_MASK_AND_RETAINED_AXES_INDEPENDENTLY_ACCEPTED_NOT_NEW_IDENTITY',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_selected_fit=alln,n_sources=len(index),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json')))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
