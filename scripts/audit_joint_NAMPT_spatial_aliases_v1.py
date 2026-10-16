"""Fresh feature-table identifier audit; no expression zeros or biological negative claims."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,re
import numpy as np,pandas as pd,anndata as ad
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
D=V/'joint_NAMPT_DDIT3_spatial_measurement_20261016_v1';di=json.loads((D/'inventory.json').read_text());p=D/'public_all14_joint_axis_measurement.tsv';assert sha(p)==di['public_output_sha256'];table=pd.read_csv(p,sep='\t');M=R/'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv';mf=pd.read_csv(M,sep='\t');S=V/'state_spatial_proxy_calibration_inputs_20261016_v1';sp=S/'fixed_all14_panel_gene_axes.json';panels=json.loads(sp.read_text());alias={'NAMPT','PBEF','PBEF1','VF','VISFATIN','1110035O14RIK','ENSG00000105835','HGNC:30092','10135','P43490'};checks=[]
for row in table.itertuples():
 match=mf[mf.sample_id.eq(row.sample_id)];assert len(match)==1;path=Path(match.iloc[0].analysis_path);expected=next(x for x in panels if x['sample_id']==row.sample_id);assert sha(path)==expected['input_sha256'];a=ad.read_h5ad(path,backed='r');columns=[c for c in a.var.columns if re.search('gene|symbol|name|identifier|accession|target|probe|^id$',str(c),re.I) and not re.search('mean|count|total|variance',str(c),re.I)];fields={'feature_index':list(a.var_names)};fields.update({str(c):a.var[c].astype(str).tolist() for c in columns});hits=[]
 for col,values in fields.items():
  for j,text in enumerate(values):
   s=str(text).upper();tokens={s,re.sub(r'\.\d+$','',s)}|set(re.split('[,;| _/]+',s));tokens.update(re.sub(r'\.\d+$','',v) for v in list(tokens))
   if tokens&alias:hits.append((j,col,sorted(tokens&alias)))
 matched_features=len({x[0] for x in hits});assert bool(matched_features)==bool(row.NAMPT_measured),('ALIAS_MISMATCH_REQUIRES_REVIEW',row.sample_id,matched_features);checks.append(dict(sample_id=row.sample_id,features=int(a.n_vars),identifier_fields=';'.join(fields),NAMPT_alias_or_ID_feature_hits=matched_features,full_source_file_SHA_verified=True,missing_gene_not_expression_zero=True));a.file.close()
O=V/'joint_NAMPT_spatial_alias_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);pd.DataFrame(checks).to_csv(O/'public_all14_fresh_NAMPT_alias_checks.tsv',sep='\t',index=False);result=dict(status='PASS_ALL14_FRESH_FEATURE_TABLES_NAMPT_SYMBOL_ALIAS_ENSEMBL_HGNC_ENTREZ_UNIPROT_MATCHES_PREFLIGHT',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),measurement_inventory_sha256=sha(D/'inventory.json'),official_identifier_source='https://www.ncbi.nlm.nih.gov/gene/10135/',official_source_review_date='2026-10-16',identifiers_checked=sorted(alias),assays=14,assays_with_NAMPT_identifier=sum(x['NAMPT_alias_or_ID_feature_hits']>0 for x in checks),public_output_sha256=sha(O/'public_all14_fresh_NAMPT_alias_checks.tsv'),gene_absence_is_not_expression_negative=True,spatial_biological_validation=False);(O/'inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
