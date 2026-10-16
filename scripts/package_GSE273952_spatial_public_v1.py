from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,tarfile,shutil
import pandas as pd,numpy as np
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'GSE273952_spatial_public_20261016_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();deps=['GSE273952_frozen_programme_independent_audit_20261016_v1','GSE273952_spatial_geometry_audit_20261016_v1','GSE273952_spatial_RNA_descriptors_audit_20261016_v1']
for tick in range(480):
 if all((V/n/'inventory.json').exists() for n in deps):break
 for name in deps+['GSE273952_frozen_programme_projection_20261016_v1','GSE273952_spatial_RNA_descriptors_20261016_v1']:
  e=V/(name+'.execution.json')
  if e.exists() and json.loads(e.read_text()).get('exit_code') not in [None,0]:raise RuntimeError('Upstream failure '+name)
 time.sleep(30)
else:raise TimeoutError('Four-hour publication wait')
for name in deps:assert json.loads((V/name/'inventory.json').read_text())['status'].startswith('PASS_ALL')
assert not O.exists();O.mkdir(mode=0o700);records=[];meta=[]
def copy(src,name):
 dst=O/name;shutil.copy2(src,dst);records.append(dict(file=name,bytes=dst.stat().st_size,sha256=sha(dst)))
for name in ['GSE273952_frozen_programme_projection_20261016_v1','GSE273952_spatial_neighbourhood_preflight_20261016_v1','GSE273952_spatial_RNA_descriptors_20261016_v1']+deps:
 src=V/name/'inventory.json';j=json.loads(src.read_text());j.pop('private_outputs',None);j['complete_server_inventory_sha256']=sha(src);dst=O/(name+'_scalar.json');dst.write_text(json.dumps(j,indent=2)+'\n');meta.append(dst);records.append(dict(file=dst.name,bytes=dst.stat().st_size,sha256=sha(dst)))
ids=set(r['programme_id'] for r in json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text())['programmes']);ids.update(['RNA_DDIT3_TARGETS_POS','RNA_DDIT3_TARGETS_NEG','RNA_DDIT3_TARGETS_SIGNED']);models={'quality_author_regions','plus_incomplete_legacy6_composition','plus_spatial_polynomial_gradient'};schemas=[]
# Explicit schema whitelists; no file-prefix heuristic or recursive private directory scan.
tables=[('GSE273952_frozen_programme_projection_20261016_v1','public_all817_programme_measurement_scope.tsv'),('GSE273952_frozen_programme_projection_20261016_v1','public_all817_cohort_RNA_score_summary.tsv'),('GSE273952_spatial_neighbourhood_preflight_20261016_v1','public_author_region_coverage_aggregate.tsv'),('GSE273952_spatial_RNA_descriptors_20261016_v1','public_all2451_programme_spatial_descriptors.tsv'),('GSE273952_spatial_RNA_descriptors_20261016_v1','public_all9_neighbour_NAMPT_target_RNA_descriptors.tsv')]
allowed_columns={'programme_id','role','signature_genes','measured_physical_genes','measurement_gate','allowed_single_cell_majors','spot_cell_identity_assigned','formal_state','valid_spots','mean_score','q10','median_score','q90','interpretation','annotation','region','spots','nonfolded_spots','libraries_with_region','model','endpoint','libraries_total','libraries_estimable','positive_libraries','negative_libraries','median','q25','q75','minimum','maximum','leave_one_library_median_min','leave_one_library_median_max','inference_pvalue_computed','cell_intrinsic_or_mechanism_validation'}
major_allowed=set()
for r in json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text())['programmes']:major_allowed.update(r['allowed_majors'])
region_allowed={'UNKNOWN_EMPTY','Folded area','Folded_area','NO_TLS','TLS_PFL','TLS_SFL','TLS_immature','TLS_sain','TLS_CD21NA_CD23NA','TLS_CD21NA_CD23neg','TLS_NA','NO_Tumor','NO_tumor','Tumor','Tumor zone'}
closed={'programme_id':ids,'endpoint':ids,'model':models,'role':{'frozen_single_cell_candidate_programme_on_mixed_spots','fixed_prior_target_RNA_proxy_not_original_ULM_or_protein'},'interpretation':{'mixed_spot_RNA_programme_not_cell_identity_or_protein_activity'},'annotation':{'TLS_primary','TLS_secondary','tumor_zone'},'region':region_allowed}
for source,name in tables:
 p=V/source/name;t=pd.read_csv(p,sep='\t');assert set(t.columns)<=allowed_columns
 for col in t.columns:
  if col in closed:assert set(t[col].dropna())<=closed[col]
  elif col=='allowed_single_cell_majors':
   for value in t[col].dropna():assert value=='not_cell_intrinsic_in_mixed_spots' or set(value.split(';'))<=major_allowed
  elif t[col].dtype==bool:pass
  else:vals=pd.to_numeric(t[col],errors='raise').dropna();assert np.isfinite(vals).all()
 copy(p,name);schemas.append(dict(file=name,rows=len(t),columns=list(t.columns),all_text_closed_domains=True,no_participant_spot_library_rows_or_clinical_outcomes=True))
review=dict(status='PASS_EXPLICIT6_SCALAR5_COHORT_OR_REGION_AGGREGATES_ALL_TEXT_CLOSED_DOMAINS_NUMERIC_FIELDS',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),files=records,table_schemas=schemas,private_spot_and_per_library_rows_not_included=True,raw_expression_or_clinical_data_not_included=True);rp=O/'public_content_review.json';rp.write_text(json.dumps(review,indent=2)+'\n');archive=Path('/tmp/GSE273952_spatial_public_content_whitelisted_20261016_v1.tar.gz');assert not archive.exists()
with tarfile.open(archive,'w:gz') as t:
 for r in records:t.add(O/r['file'],arcname=r['file'])
 t.add(rp,arcname=rp.name)
print(json.dumps(dict(archive_name=archive.name,bytes=archive.stat().st_size,sha256=sha(archive),explicit_files=len(records)+1,privacy_review_sha256=sha(rp))))
