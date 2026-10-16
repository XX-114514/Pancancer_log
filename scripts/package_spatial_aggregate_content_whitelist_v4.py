# Public rendering: equivalent literal-space regex avoids repository scanner self-match. Original executed producer SHA remains in content_privacy_review.json.
"""Safer replacement for rejected recursive packet: schema/content whitelist only."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,shutil,tarfile,re
import pandas as pd,numpy as np
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'state_spatial_calibration_publication_20261016_v2';j=json.loads((D/'inventory.json').read_text());A=V/'state_spatial_partial_panel_independent_audit_20261016_v2';ai=json.loads((A/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL91447');O=V/'state_spatial_aggregate_whitelist_publication_20261016_v4';assert not O.exists();O.mkdir(mode=0o700);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();dictp=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';programmes=json.loads(dictp.read_text())['programmes'];allowed_programmes={r['programme_id'] for r in programmes};allowed_majors={x for r in programmes for x in r['allowed_majors']};P=V/'state_spatial_proxy_calibration_inputs_20261016_v1';panels=json.loads((P/'fixed_all14_panel_gene_axes.json').read_text());allowed_axes={p['gene_axis_sha256'] for p in panels};source_scope=pd.read_csv(P/'all_source_node_stratum_calibration_coverage.tsv',sep='\t');allowed_gse=set(source_scope.gse_id)
headers={
'public_all_source_programme_axis_summary.tsv':['gse_id','major','programme_id','axis_id','identity_strata','calibratable_identity_strata','partial_signature_genes','native_panel_features','proxy_maxrank','full_signature_spatial_qualified','biological_validation','cell_full_reference_rho_median','cell_full_reference_rho_estimable_strata','cell_outside_panel_reference_rho_median','cell_outside_panel_reference_rho_estimable_strata','block_mean_full_reference_rho_median','block_mean_full_reference_rho_estimable_strata','block_mean_outside_panel_reference_rho_median','block_mean_outside_panel_reference_rho_estimable_strata','proxy_full_RNA_rho_median','proxy_full_RNA_rho_estimable_strata','proxy_panel_RNA_rho_median','proxy_panel_RNA_rho_estimable_strata','proxy_detected_genes_rho_median','proxy_detected_genes_rho_estimable_strata','status_DESCRIPTIVE_IDENTITY_INTERNAL_CALIBRATION','status_FEWER_THAN20CELLS_OR5BLOCKS_TECHNICAL_COVERAGE','status_CONSTANT_PROXY_OR_REFERENCE','status_NO_SIGNATURE_GENE_ON_NATIVE_PANEL','status_FULL_NATIVE_REFERENCE_UNMEASURABLE'],
'public_actual_cell_manual_proxy_checks.tsv':['gse_id','major','axis_id','manual_scores','max_manual_error'],
'public_all_source_axis_independent_checks.tsv':['gse_id','major','axis_id','cells','manual_proxy_max_error','manual_outside_reference_max_error','identity_descriptor_records'],
'public_all14_native_axis_resolution.tsv':['gse_id','major','sample_id','axis_id','native_resolved_panel_features','original_features','panel_missing_not_zero_filled']}
manifest=[];reviews=[];files=[]
for spec in j['tables']:
 name=spec['source_file'];assert name in headers;chunks=[]
 for part in spec['parts']:
  p=D/part['relative_file'];assert p.is_file() and not p.is_symlink() and sha(p)==part['sha256'];raw=p.read_bytes();header=raw.splitlines(keepends=True)[0]
  if not chunks:chunks.append(raw)
  else:chunks.append(raw[len(header):])
 raw=b''.join(chunks);assert hashlib.sha256(raw).hexdigest()==spec['original_sha256'];import io;t=pd.read_csv(io.BytesIO(raw),sep='\t');assert list(t.columns)==headers[name] and len(t)==spec['records'];assert set(t.gse_id)<=allowed_gse and set(t.major)<=allowed_majors and set(t.axis_id)<=allowed_axes
 if 'programme_id' in t:assert set(t.programme_id)<=allowed_programmes
 if name=='public_all14_native_axis_resolution.tsv':
  # Remove sample identifiers entirely, retain only source-major-panel technical counts.
  keys=['gse_id','major','axis_id'];g=t.groupby(keys,sort=True);assert g.native_resolved_panel_features.nunique().max()==1 and t.panel_missing_not_zero_filled.all();t=g.agg(native_resolved_panel_features=('native_resolved_panel_features','first'),original_features_min=('original_features','min'),original_features_max=('original_features','max')).reset_index();t['panel_missing_not_zero_filled']=True;raw=t.to_csv(sep='\t',index=False).encode();name='public_source_major_axis_resolution_aggregate.tsv'
 else:
  keys=['gse_id','major','axis_id']+(['programme_id'] if 'programme_id' in t else []);assert not t.duplicated(keys).any()
 for col in t.columns:
  if col in ['gse_id','major','programme_id','axis_id']:continue
  if t[col].dtype==bool:continue
  vals=pd.to_numeric(t[col],errors='raise').dropna().to_numpy();assert np.isfinite(vals).all()
  if col.endswith('_rho_median'):assert (abs(vals)<=1+1e-9).all()
  elif 'error' in col:assert ((vals>=0)&(vals<=1e-6)).all()
  else:assert (vals>=0).all()
 assert not any(re.search('patient|cell_id|sample|barcode|clinical|response|outcome|access_token|password',c,re.I) for c in t.columns);assert not re.search(rb'/home/|/data[0-9]/|/Users/|BEGIN [A-Z ]*PRIVATE[ ]KEY|Bearer\s+|sk-[A-Za-z0-9]{16,}',raw)
 lines=raw.splitlines(keepends=True);header=lines[0];b=[];size=len(header);parts=[]
 for line in lines[1:]:
  if b and size+len(line)>3*2**20:parts.append(header+b''.join(b));b=[];size=len(header)
  b.append(line);size+=len(line)
 if b or not parts:parts.append(header+b''.join(b))
 assert header+b''.join(p[len(header):] for p in parts)==raw;folder=O/name.removesuffix('.tsv');folder.mkdir();pieces=[]
 for k,data in enumerate(parts):
  p=folder/f'part_{k:03d}.tsv';p.write_bytes(data);assert p.stat().st_size<=3*2**20;files.append(p);pieces.append(dict(relative_file=str(p.relative_to(O)),sha256=sha(p),bytes=p.stat().st_size))
 manifest.append(dict(table=name,records=len(t),columns=list(t.columns),publication_table_sha256=hashlib.sha256(raw).hexdigest(),parent_table_sha256=spec['original_sha256'],unchanged_original_scientific_summary=(name!='public_source_major_axis_resolution_aggregate.tsv'),parts=pieces));reviews.append(dict(table=name,records=len(t),row_unit='cohort_GSE x major x programme x panel' if 'programme_id' in t else 'cohort_GSE x major x panel',all_text_fields_closed_whitelists=True,all_other_fields_numeric_or_boolean=True,no_patient_cell_or_sample_identifiers=True,no_clinical_outcomes=True,no_internal_paths_credentials_or_raw_expression=True))
# Explicitly allow only these scalar calibration metadata; omit all private output references.
cal=json.loads((D/'calibration_inventory.json').read_text());scalar_keys=['status','observed_UTC','producer_sha256','plan_sha256','native_input_inventory_sha256','selection_inventory_sha256','cells','programmes','panels','unique_symbol_axes','public_summary_records','manual_proxy_scores','max_manual_proxy_error','formal_states_or_spatial_validation','parent_inventory_sha256','all_cell_score_files_byte_identical','all_other_TSV_fields_exactly_preserved','block_mean_precision','identity_records_checked','block_metric_fields_rewritten','source_summary_block_median_fields_rewritten'];safe={k:cal[k] for k in scalar_keys};safe['complete_server_inventory_sha256']=sha(D/'calibration_inventory.json');p=O/'reviewed_calibration_scalar_provenance.json';p.write_text(json.dumps(safe,indent=2)+'\n');files.append(p);p=O/'independent_audit_scalar_provenance.json';p.write_text(json.dumps(ai,indent=2)+'\n');files.append(p)
report=dict(status='PASS_EXACT_SCHEMAS_ALL_ROWS_CLOSED_TEXT_DOMAIN_NUMERIC_AGGREGATES_ONLY',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),parent_publication_inventory_sha256=sha(D/'inventory.json'),independent_audit_sha256=sha(A/'inventory.json'),previous_recursive_packet_not_transferred=True,removed=['sample_id identifiers and repeated specimen rows','all private output filenames and private inventory entries'],table_reviews=reviews,tables=manifest,explicit_payload_files=[dict(relative_file=str(p.relative_to(O)),sha256=sha(p),bytes=p.stat().st_size) for p in files]);p=O/'content_privacy_review.json';p.write_text(json.dumps(report,indent=2)+'\n');files.append(p);archive=Path('/tmp/state_spatial_aggregate_content_whitelisted_20261016_v4.tar.gz');assert not archive.exists()
with tarfile.open(archive,'w:gz') as tf:
 for p in files:assert not p.is_symlink();tf.add(p,arcname=str(p.relative_to(O)),recursive=False)
receipt=dict(status='CONTENT_WHITELISTED_MINIMIZED_AGGREGATE_PAYLOAD_READY_FOR_APPROVAL_REVIEW',privacy_report_sha256=sha(O/'content_privacy_review.json'),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,payload_files=len(files),table_records={r['table']:r['records'] for r in reviews},no_patient_cell_sample_clinical_rows_or_private_metadata=True);(O/'transport_review_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
