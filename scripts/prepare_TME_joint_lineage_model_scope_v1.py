"""Append-only model-input mask for full joint parent clusters and controls.

This bridge reuses the standard integration input protocol. It is not a new
cell identity release, QC decision, or deletion of any out-of-model record.
"""
import argparse,datetime,hashlib,json,sys,time
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
 cfg=load(R/'config/spec_v1.json');H=Path(cfg['parent_producer']);G=Path(cfg['parent_graph']);C=Path(cfg['coexpression_producer']);A=Path(cfg['coexpression_auditor']);wait=A/'provenance/supervisor_terminal_v1.json'
 while not wait.exists():time.sleep(5)
 for D in [H,G,C,A]:assert load(D/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'outputs/accepted_joint_cluster_coexpression_v1.json')['producer_terminal_SHA256']==sha(C/'provenance/terminal_v1.json');hc=load(H/'config/spec_v1.json');sys.path.insert(0,str(Path(hc['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);original_acceptance=Path(hc['working_acceptance']);assert sha(original_acceptance)==hc['working_acceptance_SHA256'];accepted=load(original_acceptance);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];original_manifest=load(accepted['producer_manifest']['path']);W=Path(original_manifest['index']['path']).parent;index=load(original_manifest['index']['path']);assert sha(original_manifest['index']['path'])==original_manifest['index']['SHA256']
 index=sorted(index,key=lambda r:r['source']);graph_axes=G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz';labels=pd.read_csv(graph_axes,sep='\t',usecols=['leiden_0_5']).iloc[:,0].to_numpy(np.int32)
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 sources=load(H/'outputs/scope_and_common_genes_v1.json')['sources'];assert [r['source'] for r in index]==sources;core=np.asarray(cfg['core_parent_clusters'],np.int32);control=np.asarray(cfg['whole_parent_control_clusters'],np.int32);assert not set(core)&set(control);selected=np.isin(labels,np.r_[core,control]);assert int(selected.sum())==cfg['expected_selected_n'];parent_positions=np.flatnonzero(selected);np.save(R/'outputs/private_parent_fit_positions_v1.npy',parent_positions);new_index=[];records=[]
 for code,record in enumerate(index):
  source=record['source'];p=W/record['arrays'];metadata=W/record['metadata'];assert sha(p)==record['arrays_SHA256'] and sha(metadata)==record['metadata_SHA256']
  with np.load(p,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','fit_pool_candidate_not_final_QC','source_qualification_code','native_total_RNA']}
  fit=np.flatnonzero(q['fit_pool_candidate_not_final_QC']);wh=np.flatnonzero(axes['source_code']==code);np.testing.assert_array_equal(q['native_row'][fit],axes['native_row'][wh]);np.testing.assert_array_equal(q['source_native_position'][fit],axes['source_native_position'][wh]);new_mask=np.zeros(len(q['native_row']),bool);new_mask[fit[selected[wh]]]=True;assert int(new_mask.sum())==cfg['expected_per_source'][source]
  old_mask=q['fit_pool_candidate_not_final_QC'];q['parent_fit_pool_candidate_not_final_QC']=old_mask;q['fit_pool_candidate_not_final_QC']=new_mask;out=R/'outputs'/(source+'_lineage_model_input_mask_v1.npz');np.savez_compressed(out,**q);new_index.append(dict(source=source,metadata=str(metadata),metadata_SHA256=record['metadata_SHA256'],arrays=out.name,arrays_SHA256=sha(out),protocol_role='model_row_mask_only_not_new_identity_or_QC'))
  records.append(dict(source=source,n_parent_qualified=len(new_mask),n_parent_fit=len(fit),n_selected_fit=int(new_mask.sum()),n_core_fit=int(np.isin(labels[wh],core).sum()),n_whole_parent_controls_fit=int(np.isin(labels[wh],control).sum()),all_original_qualified_rows_retained=True))
 ip=R/'outputs/lineage_model_input_index_v1.json';save(ip,new_index);mp=R/'outputs/lineage_model_input_manifest_v1.json';save(mp,dict(schema='joint_parent_lineage_model_mask_bridge_v1',index=dict(path=str(ip),SHA256=sha(ip)),original_identity_acceptance=dict(path=str(original_acceptance),SHA256=sha(original_acceptance)),core_parent_clusters=core.tolist(),whole_parent_control_clusters=control.tolist(),parent_graph_axes_SHA256=sha(graph_axes),parent_fit_axes_SHA256=sha(H/'outputs/private_cell_axes_v1.npz'),coexpression_acceptance_SHA256=sha(A/'outputs/accepted_joint_cluster_coexpression_v1.json'),parent_fit_positions=dict(path=str(R/'outputs/private_parent_fit_positions_v1.npy'),SHA256=sha(R/'outputs/private_parent_fit_positions_v1.npy')),source_scope_records=records,n_selected_fit=int(selected.sum()),all_selected_parent_cells_included_not_only_marker_positive=True,all_original_rows_and_identities_and_QC_retained=True,transductive_joint_cluster_refinement_not_ABC_winner=True))
 bridge=R/'outputs/accepted_lineage_model_input_bridge_v1.json';save(bridge,dict(schema='model_scope_adapter_not_new_working_identity_acceptance_v1',producer_manifest=dict(path=str(mp),SHA256=sha(mp)),original_identity_acceptance_reference=dict(path=str(original_acceptance),SHA256=sha(original_acceptance)),independent_mask_audit_still_required=True));storage_gate(R,cfg)
 save(R/'provenance/terminal_v1.json',dict(status='FULL_PARENT_LINEAGE_MODEL_INPUT_MASK_READY_REQUIRES_INDEPENDENT_MASK_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_selected_fit=int(selected.sum()),n_sources=len(index),bridge_SHA256=sha(bridge),manifest_SHA256=sha(mp)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
