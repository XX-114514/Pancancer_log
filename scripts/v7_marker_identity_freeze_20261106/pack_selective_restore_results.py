import json,zlib,hashlib
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';h=r/'stage_h_selective_subtype_restore_v4';d=r/'stage_i_selective_annotation_marker_evidence_v4'
assert json.loads((h/'independent_audit.json').read_text())['status']=='PASS_selective_subtype_preservation_and_complete_export';assert json.loads((d/'independent_naming_marker_audit.json').read_text())['status']=='PASS';files={}
t=pd.read_pickle(h/'all_axis_selective_annotation_with_history_PRIVATE.pkl.gz');q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');assert t.index.equals(q.index);qc=t.QC_keep.to_numpy();selected=t.annotation_frozen.to_numpy();a=pd.DataFrame({'source':q.gse_id.astype(str),'old_label':t.frozen_label.astype(str),'selected':selected,'retained':t.prior_named_annotation_retained.to_numpy(),'n':np.ones(len(t),np.int8)},index=t.index)
ret=a.loc[qc].groupby('source',observed=True).agg(n_QC_keep=('n','sum'),n_selected=('selected','sum'),n_previous_named_retained=('retained','sum')).reset_index();ret['n_outside']=ret.n_QC_keep-ret.n_selected;assert ret.n_QC_keep.sum()==4668312 and ret.n_selected.sum()==selected.sum();ret.to_csv(h/'source_retention_v4.tsv',sep='\t',index=False)
out=a.loc[qc&~selected].groupby('old_label',observed=True).size().rename('n_QC_keep_outside').reset_index();assert out.n_QC_keep_outside.sum()==606373;out.to_csv(h/'outside_by_previous_label_v4.tsv',sep='\t',index=False)
retained=t.prior_named_annotation_retained.to_numpy();context=pd.DataFrame({'retained_annotation':t.final_annotation.astype(str),'historical_v3_broad_screen_status':t.freeze_status.astype(str),'v3_positive_parent_diagnostic':t.broad_v3_parent_diagnostic.astype(str),'n':np.ones(len(t),np.int8)},index=t.index);context.loc[retained].groupby(['retained_annotation','historical_v3_broad_screen_status','v3_positive_parent_diagnostic'],observed=True).n.sum().reset_index(name='n_cells').to_csv(h/'retained_annotations_broad_screen_context.tsv',sep='\t',index=False)
context.loc[t.broad_parent_difference_needs_review].groupby(['retained_annotation','v3_positive_parent_diagnostic'],observed=True).n.sum().reset_index(name='n_cells').to_csv(h/'retained_annotations_parent_disagreement_for_review.tsv',sep='\t',index=False)
for folder in [h,d]:
 for p in list(folder.glob('*.tsv'))+list(folder.glob('*.json')):
  assert 'PRIVATE' not in p.name;b=p.read_bytes();files[str(p.relative_to(r))]={'sha256':hashlib.sha256(b).hexdigest(),'text':b.decode()}
for name in ['restore_subtypes_selectively_v4.py','audit_selective_restore_v4.py','aggregate_selective_markers_v4.py','audit_selective_markers_v4.py']:
 p=r/'launch'/(name+'_terminal.json');b=p.read_bytes();files[str(p.relative_to(r))]={'sha256':hashlib.sha256(b).hexdigest(),'text':b.decode()}
payload=zlib.compress(json.dumps({'files':files}).encode(),9);p=r/'public_transfer/selective_subtype_restore_v4_packet.zlib';p.open('xb').write(payload)
print(json.dumps({'packet':str(p.relative_to(PROJECT_ROOT)),'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),'n_files':len(files)}))
