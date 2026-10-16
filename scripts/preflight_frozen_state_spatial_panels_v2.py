"""Measurement-only bridge: all frozen814 programmes vs existing audited spatial axes."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,anndata as ad

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)=='6a3c4311d89a9a02d55a19f2d249c960c508b5973368c06c2e5c43bb605edec9';records=json.loads(D.read_text())['programmes'];assert len(records)==814;P=R/'immune_same_name_family_source_review_v1';M=P/'spatial_native_input_repair_v2/analysis_input_manifest.tsv';assert sha(M)=='6fda649c8b2bf4abe010ba8d53d50e495bae29ce8d30648675d9be9bf91dcf83';A=P/'spatial_contact_native_preflight_v3/independent_audit.json';audit=json.loads(A.read_text());assert audit['status']=='PASS_ALL65_NATIVE_VALUES_CONTROLS_QC_BACKGROUND_LINEAGE_OVERLAP_AUTHOR_AND_FOV_RECOUNT';I=P/'spatial_contact_native_preflight_v3/inventory.json';assert sha(I)==audit['producer_inventory_sha256'];prior=json.loads(I.read_text());files=pd.read_csv(M,sep='\t');assert len(files)==14;C=V/'original_CM_programme_pair_associations_20261015_v2';ci=json.loads((C/'inventory.json').read_text());dp=C/'all_programme_pair_associations.tsv';assert sha(dp)==ci['output_sha256'][dp.name];s=pd.read_csv(dp,sep='\t');s=s[s.BH_q<.05].copy();assert len(s)==36;out=[];pair=[];axes=[]
 for f in files.to_dict('records'):
  path=Path(f['analysis_path']);previous=next((x for x in prior['inputs'] if x['sample_id']==f['sample_id']),None);actual_sha=sha(path);assert previous is None or actual_sha==previous['sha256'];expected={'sha256':actual_sha};a=ad.read_h5ad(path,backed='r');assert a.var_names.is_unique;lookup={}
  for col,index in enumerate(a.var_names):
   names={str(index)}
   if 'gene_symbol' in a.var:names.add(str(a.var.iloc[col]['gene_symbol']))
   for gene in names:
    if gene and gene not in ['nan','None']:lookup.setdefault(gene,set()).add(col)
  rows={}
  for r in records:
   genes=r['primary_genes'];present=[g for g in genes if len(lookup.get(g,[]))==1];ambiguous=[g for g in genes if len(lookup.get(g,[]))>1];missing=[g for g in genes if g not in present and g not in ambiguous];gate=len(present)>=10 and len(present)/len(genes)>=.8
   # Independent count via feature-index union resolved per symbol.
   count=sum(1 for g in genes if g in lookup and len(lookup[g])==1);assert count==len(present)
   row=dict(dataset_id=f['dataset_id'],sample_id=f['sample_id'],spatial_cells=int(a.n_obs),spatial_features=int(a.n_vars),programme_id=r['programme_id'],method=r['method'],algorithm_stable_candidate=r['algorithm_stable_candidate'],allowed_majors=';'.join(r['allowed_majors']),signature_genes=len(genes),unambiguously_measured=len(present),measured_fraction=len(present)/len(genes),ambiguous_genes=';'.join(ambiguous),missing_genes=';'.join(missing),programme_measurement_gate=gate,original_exact_endpoint_identity_validated=False,spatial_score_computed=False,biological_validation=False);out.append(row);rows[r['programme_id']]=row
  for r in s.itertuples():
   aa=rows[r.programme_A];bb=rows[r.programme_B];pair.append(dict(sample_id=f['sample_id'],dataset_id=f['dataset_id'],source_association_context=r.context,source_pair_index=r.pair_index,endpoint_A=r.endpoint_A,endpoint_B=r.endpoint_B,programme_A=r.programme_A,programme_B=r.programme_B,measured_A=aa['unambiguously_measured'],measured_B=bb['unambiguously_measured'],gate_A=aa['programme_measurement_gate'],gate_B=bb['programme_measurement_gate'],both_measurement_gates=aa['programme_measurement_gate'] and bb['programme_measurement_gate'],exact_endpoint_identity_validated=False,inference='MEASUREMENT_ONLY_NO_SPATIAL_TEST'))
  axes.append(dict(dataset_id=f['dataset_id'],sample_id=f['sample_id'],cells=int(a.n_obs),features=int(a.n_vars),input_sha256=expected['sha256'],measurable_programmes=sum(x['programme_measurement_gate'] for x in rows.values())));a.file.close()
 O=V/'frozen_state_spatial_panel_preflight_20261016_v2';assert not O.exists();O.mkdir();pd.DataFrame(out).to_csv(O/'all814_by14_assay_programme_coverage.tsv',sep='\t',index=False);pd.DataFrame(pair).to_csv(O/'all36_by14_assay_pair_measurement.tsv',sep='\t',index=False);pd.DataFrame(axes).to_csv(O/'all14_spatial_axes.tsv',sep='\t',index=False);j=dict(status='COMPLETE_FROZEN_STATE_SPATIAL_MEASUREMENT_PREFLIGHT_NOT_SPATIAL_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),dictionary_sha256=sha(D),spatial_manifest_sha256=sha(M),prior_spatial_audit_sha256=sha(A),programmes=814,assays=14,coverage_records=len(out),pair_measurement_records=len(pair),gene_gate_min_fraction=.8,gene_gate_min_measured=10,fully_measurable_programme_assay_records=sum(r['programme_measurement_gate'] for r in out),both_programmes_measurable_pair_assay_records=sum(r['both_measurement_gates'] for r in pair),exact_endpoint_identity_validated=False,spatial_test_performed=False,additional9_count_QC_not_revalidated_here=True,output_sha256={p.name:sha(p) for p in O.glob('*.tsv')});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
if __name__=='__main__':main()
