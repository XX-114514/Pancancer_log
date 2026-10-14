"""All24 audited scopes to complete candidate signatures, no formal states."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def frame(p):
 with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def markers(values,genes):
 ids=np.argsort(-values,kind='stable');ids=ids[values[ids]>0][:50];return np.asarray(genes)[ids].tolist()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path(os.environ['PANCANCER_ALL_STATE_DICTIONARY_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];S=V/'remaining_state_diagnostic_supervision_20261014_v2';sv=json.loads((S/'inventory.json').read_text());assert sv['status']=='COMPLETE_ALL6_BLOCKS_AND18_MAJORS_FIXED_NUMERIC_DIAGNOSTICS_NOT_STATES' and sv['completed_scopes']==24
accepted=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(accepted)==plan['accepted_input_inventory_sha256'];scope=json.loads(accepted.read_text());majors=sorted(r['major'] for r in scope['files']);records=[];parent_checks=[];correspondence=[];spectra_reports=[]
for i,major in enumerate(majors):
 if i==0:
  D=R/'all_lineage_cNMF_diagnostics_20261014_v1/major_00';A=R/'all_lineage_cNMF_diagnostic_audit_20261014_v2/major_00/independent_audit.json'
 else:D=V/'all_lineage_cNMF_diagnostics_20261014_v2'/('major_%02d'%i);A=V/'all_lineage_cNMF_diagnostic_audit_20261014_v3'/('major_%02d'%i)/'independent_audit.json'
 iv=json.loads((D/'inventory.json').read_text());a=json.loads(A.read_text());assert a['status'].startswith('PASS_ALL_MAJOR') and a['producer_inventory_sha256']==sha(D/'inventory.json');assert iv['major']==major and iv['candidates']==45 and iv['repetition_matches']==4500 and iv['crossK_comparisons']==650
 for name,h in iv['output_sha256'].items():assert sha(D/name)==h
 parent_checks.append(dict(role='cNMF',scope=major,diagnostic_inventory_sha256=sha(D/'inventory.json'),audit_sha256=sha(A)))
 table=pd.read_csv(D/'all45_candidates.tsv',sep='\t',keep_default_na=False);CX=pd.read_csv(D/'all650_crossK_comparisons.tsv',sep='\t');B=V/'all_lineage_cNMF_state_contrast_20261014_v2'/('major_%02d'%i)/'cnmf/cnmf_tmp';bases={}
 for k in plan['K']:
  cp=B/f'cnmf.gene_spectra_score.k_{k}.dt_0_5.df.npz';tp=B/f'cnmf.gene_spectra_tpm.k_{k}.dt_0_5.df.npz';bases[k]=(frame(cp),frame(tp));assert sha(cp)==iv['input_output_sha256'][str(cp.relative_to(B.parent.parent))];assert sha(tp)==iv['input_output_sha256'][str(tp.relative_to(B.parent.parent))]
 stable={(r.K,r.component):bool(r.repeat_and_crossK_candidate) for r in table.itertuples()}
 for row in table.itertuples():
  k=int(row.K);component=int(row.component);coef,tpm=bases[k];genes=coef.columns.to_numpy(str);primary=markers(np.maximum(coef.loc[component].values,0),genes);alternate=markers(tpm.loc[component].values,genes);assert primary==str(row.primary50_positive_coefficient_genes).split(';') or (not primary and row.primary50_positive_coefficient_genes=='')
  qualified=bool(row.repeat_and_crossK_candidate);records.append(dict(programme_id=row.candidate_id,method='cNMF',allowed_majors=[major],K=k,component=component,primary_genes=primary,TPM_sensitivity_genes=alternate,algorithm_stable_candidate=qualified,stability_main_candidate=qualified,native_training_supported_cells=int(row.informative_gene_support_cells),native_training_supported_fraction=float(row.informative_gene_support_fraction),signature_RNA_fraction_median=float(row.signature_RNA_fraction_median),source_r2=float(row.source_r2),source_node_r2=float(row.source_node_r2),source_node_increment=float(row.source_node_r2-row.source_r2),RNA_partial_rank_after_source_node=float(row.RNA_partial_rank_after_source_node),detected_partial_rank_after_source_node=float(row.detected_partial_rank_after_source_node),source_limited=bool(row.source_limited),gene_signature_origin='top50_genuinely_positive_native_fullRNA_TPM_regression_coefficients',no_zero_coefficients_padded=True,crossK_independent_state_count=False,formal_state_name_assigned=False,clinical_LR_p_used=False))
 for row in CX.itertuples():
  if row.positive_coef_cosine>=.8 and row.top50_jaccard>=.25 and row.source_residual_usage_rho>=.5:
   correspondence.append(dict(major=major,programme_A='ACN_'+major.replace(' ','_').replace('/','_')+f'_K{int(row.K_a)}_{int(row.component_a):02}',programme_B='ACN_'+major.replace(' ','_').replace('/','_')+f'_K{int(row.K_b)}_{int(row.component_b):02}',both_algorithm_stable=stable[int(row.K_a),int(row.component_a)] and stable[int(row.K_b),int(row.component_b)],coefficient_cosine=float(row.positive_coef_cosine),top_gene_jaccard=float(row.top50_jaccard),source_residual_usage_rho=float(row.source_residual_usage_rho),correspondence_is_not_identity_merge_or_unique_state_claim=True))
spplan=json.loads(Path(os.environ['PANCANCER_ALL_LINEAGE_SPECTRA_DIAGNOSTIC_PLAN']).read_text())
for block,maj in spplan['blocks'].items():
 if block=='lymphoid_T':D=R/'all_lineage_Spectra_diagnostics_20261014_v1'/block;A=R/'all_lineage_Spectra_diagnostic_audit_20261014_v1'/block/'independent_audit.json'
 else:D=V/'all_lineage_Spectra_diagnostics_20261014_v2'/block;A=V/'all_lineage_Spectra_diagnostic_audit_20261014_v2'/block/'independent_audit.json'
 iv=json.loads((D/'inventory.json').read_text());a=json.loads(A.read_text());assert a['status'].startswith('PASS_ALL_BLOCK') and a['producer_inventory_sha256']==sha(D/'inventory.json');parent_checks.append(dict(role='Spectra',scope=block,diagnostic_inventory_sha256=sha(D/'inventory.json'),audit_sha256=sha(A)))
 for name,h in iv['output_sha256'].items():assert sha(D/name)==h
 table=pd.read_csv(D/'all_candidate_diagnostics.tsv',sep='\t',keep_default_na=False);diag=pd.read_csv(D/'all_source_node_complexity_diagnostics.tsv',sep='\t');spectra_reports.append(table.assign(dictionary_input_block=block))
 # All factors stay in the audit registry; only prespecified stable/context candidates are signatures.
 for row in table.itertuples():
  if not row.stable_context_candidate:continue
  own=diag[diag.anchor_factor.eq(row.anchor_factor)];permitted=maj if row.scope=='block_global' else [row.scope];genes=str(row.primary50_median_L2_gene_weights).split(';');records.append(dict(programme_id=row.candidate_id,method='Spectra',allowed_majors=permitted,primary_genes=genes,TPM_sensitivity_genes=[],algorithm_stable_candidate=True,stability_main_candidate=bool(row.stable_shared_or_lineage_candidate),context_candidate=True,supported_majors=str(row.supported_majors).split(';'),model_scope=row.scope,gene_signature_origin='median_L2_weights_all3_matched_seeds',prior_descriptor=row.best_prior_descriptor,prior_descriptor_is_confirmed_state=False,source_r2=float(own.source_r2.max()),source_node_r2=float(own.source_node_r2.max()),native_training_supported_cells=int(own.informative_gene_support_cells.sum()),source_limited=bool(row.source_limited_block),complexity_review=bool(row.complexity_review),formal_state_name_assigned=False,clinical_LR_p_used=False))
panels=plan['diagnostic_panels'];reference=json.loads((R/'state_model_runtime_20261014_v1/resources/spectra_snapshot/src/Spectra/data/default_gene_sets.json').read_text())
for r in records:
 gs=set(r['primary_genes']);overlap={key:len(gs&set(value)) for key,value in panels.items()};r['diagnostic_panel_overlaps']=overlap;r['descriptive_programmes']=[key for key,n in overlap.items() if n>=3];r['foreign_RNA_sensitive']=any(overlap.get(key,0)>=3 for key in ['epithelial_RNA','stromal_RNA','mural_RNA','neural_RNA']);r['other_immune_lineage_review']=any((major in ['Macrophage','Monocyte','Myeloid','Dendritic cell','Neutrophil'] and overlap.get('T_core',0)>=3) or (major in ['T cell','T/NK lymphoid','Lymphoid cell'] and overlap.get('Myeloid_core',0)>=3) for major in r['allowed_majors']);r['ribosomal_top50_count']=sum(g.startswith(('RPL','RPS')) for g in gs);r['mitochondrial_top50_count']=sum(g.startswith('MT-') for g in gs);r['basic_expression_review_flag']=r['ribosomal_top50_count']>=10 or r['mitochondrial_top50_count']>=3;r['low_native_training_support_review']=r['native_training_supported_cells']<20;r['state_or_identity_or_other_RNA_unresolved']=True;r['formal_state_labels']=False;r['cross_method_agreement_assumed']=False
 best=[]
 for scope,sets in reference.items():
  for name,rr in sets.items():
   ss=set(rr);n=len(gs&ss)
   if n:best.append((n/len(gs|ss),n,scope+'::'+name))
 q=sorted(best,key=lambda x:(-x[0],-x[1],x[2]))[0] if best else (0,0,'none');r['best_reference_descriptor']=q[2];r['best_reference_jaccard']=q[0];r['best_reference_overlap']=q[1];r['reference_descriptor_is_confirmed_state']=False;r['state_naming_status']='CANDIDATE_ONLY_NATIVE_SOURCE_AND_IDENTITY_WITHIN_NODE_REVIEW_REQUIRED';r['gene_signature_missing_is_not_zero_expression']=True
assert len([r for r in records if r['method']=='cNMF'])==810 and len({r['programme_id'] for r in records})==len(records) and len(parent_checks)==24
O=V/'all_lineage_state_candidate_dictionary_20261014_v3';assert not O.exists();O.mkdir();out=dict(schema='all18_candidate_state_dictionary_v3',frozen_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),accepted_input_inventory_sha256=sha(accepted),supervision_inventory_sha256=sha(S/'inventory.json'),input_scopes=parent_checks,programmes=records,programme_count=len(records),not_unique_states=True,no_LR_or_clinical_p_used=True,formal_states=False);(O/'frozen_programme_dictionary.json').write_text(json.dumps(out,indent=2)+chr(10));rows=[{k:(';'.join(map(str,v)) if isinstance(v,list) else json.dumps(v) if isinstance(v,dict) else v) for k,v in r.items()} for r in records];pd.DataFrame(rows).to_csv(O/'all_candidate_programme_catalogue.tsv',sep='\t',index=False);pd.DataFrame(correspondence).to_csv(O/'crossK_correspondence_not_state_merge.tsv',sep='\t',index=False);pd.concat(spectra_reports,ignore_index=True).to_csv(O/'all_Spectra_factor_audit_registry.tsv',sep='\t',index=False)
inventory=dict(status='FROZEN_ALL18_CANDIDATE_SIGNATURES_PENDING_NATIVE_PROJECTION_NOT_STATES',frozen_UTC=out['frozen_UTC'],producer_sha256=sha(__file__),plan_sha256=sha(P),cNMF_candidates=810,Spectra_candidate_signatures=sum(r['method']=='Spectra' for r in records),all_Spectra_factors=len(pd.concat(spectra_reports)),algorithm_stable_candidates=sum(r['algorithm_stable_candidate'] for r in records),main_stability_candidates=sum(r['stability_main_candidate'] for r in records),dictionary_sha256=sha(O/'frozen_programme_dictionary.json'),public_output_sha256={p.name:sha(p) for p in O.glob('*.tsv')},formal_states=False,original_labels_changed=False,all9CM_complete=False,DIALOGUE_complete=False,spatial_complete=False);(O/'inventory.json').write_text(json.dumps(inventory,indent=2)+chr(10));print(json.dumps(inventory),flush=True)
