"""Independent full-cell recount and exact source/preparation/frequency checks."""
import json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,zarr
from anndata.io import read_elem
r=PROJECT_ROOT;run=r/'runs/20261106_v7_selective_CoVarNet_paper_v1';o=run/'inputs';saved=json.loads((o/'complete.json').read_text());f=pd.read_pickle(r/'runs/20261106_v7_publication_marker_annotation_v1/stage_h_selective_subtype_restore_v4/annotation_frozen_V7_selective_subtype_20261106_v4_PRIVATE.pkl.gz');g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');axis=pd.Index(read_elem(g['obs/_index']));pos=axis.get_indexer(f.index);assert (pos>=0).all();obs={k:np.asarray(read_elem(g['obs'][k])).astype(str)[pos] for k in ['gse_id','sampleID','sample_uid']};units=pd.read_csv(o/'sample_admission_PRIVATE.tsv',sep='\t',keep_default_na=False)
v=r.parent/'Project_v5/runs';core=pd.DataFrame(json.loads((v/'20261027_source_preanalysis_gap_gate_v7/outputs/source_qualification_work_units_core_v7.json').read_text())['units']);ctx=pd.read_csv(v/'20261030_TME_CM_technical_material_design_v1/outputs/private_analysis_unit_context_and_coverage_v1.tsv',sep='\t',keep_default_na=False).set_index('qualification_unit_key_SHA256');unitidx=pd.Index(units.sampleID);cell_unit=unitidx.get_indexer(obs['sample_uid']);assert (cell_unit>=0).all();n=np.bincount(cell_unit,minlength=len(units));assert np.array_equal(n,units.n_V4_cells.to_numpy());expected=[]
for i,z in units.iterrows():
 subset=core.loc[core.source.eq(z.gse_id)];hits=[];raw_ids=set(np.asarray(obs['sampleID'])[cell_unit==i])|{z.sampleID}
 for rec in subset.to_dict('records'):
  aliases={rec['qualification_unit_id']}
  if '::' in rec['qualification_unit_id']:aliases.add(rec['qualification_unit_id'].split('::')[-1])
  aliases|={x.split('|')[-1] for x in list(aliases) if '|' in x}
  raw_alias=set(raw_ids)
  for value in raw_ids:
   for prefix in [z.gse_id+'|',z.gse_id+'::']:
    if value.startswith(prefix):raw_alias.add(value[len(prefix):])
  if aliases&raw_alias:hits.append(rec)
 accept=False
 if len(hits)==1:
  rec=hits[0];key=rec['qualification_unit_key_sha256'];assert z.qualification_unit_key_sha256==key
  if key in ctx.index:
   a=ctx.loc[key];fresh=str(a.fresh_frozen).lower();accept=rec['source_qualification']=='来源资格通过' and ('fresh' in fresh or '新鲜' in fresh) and not ('冻' in fresh or 'frozen' in fresh) and a.RNA_modality=='scRNA' and n[i]>=50
 expected.append(accept)
assert np.array_equal(expected,units.admitted_before_non_epithelial.to_numpy(bool));take=np.asarray(expected)[cell_unit]&f.major_cell_type.astype(str).ne('Epithelial').to_numpy();assert int(take.sum())==saved['n_non_epithelial_cells'];cell_saved=pd.read_pickle(o/'selected_cell_annotation_and_sample_PRIVATE.pkl.gz');assert cell_saved.index.equals(f.index[take]);assert np.array_equal(cell_saved.subCluster.to_numpy(),f.loc[take,'final_annotation'].astype(str).to_numpy())
fq=pd.read_csv(o/'frequency_PRIVATE.tsv',sep='\t',index_col=0);node=pd.Index(fq.index);samples=pd.Index(fq.columns[1:]);allcounts=pd.read_csv(o/'counts_all_non_epithelial_PRIVATE.tsv',sep='\t',index_col=0);labels=pd.Index(allcounts.index);si=samples.get_indexer(obs['sample_uid'][take]);li=labels.get_indexer(f.loc[take,'final_annotation'].astype(str));assert (si>=0).all() and (li>=0).all();count=np.bincount(li*len(samples)+si,minlength=len(labels)*len(samples)).reshape(len(labels),len(samples));assert np.array_equal(count,allcounts.to_numpy())
catalog=pd.read_csv(o/'node_catalogue.tsv',sep='\t').set_index('subCluster');parents=catalog.loc[labels,'majorCluster'];den={p:count[parents.eq(p)].sum(axis=0) for p in parents.unique()};independent=np.vstack([np.divide(count[i],den[parents.iloc[i]],out=np.zeros(len(samples),float),where=den[parents.iloc[i]]>0) for i in range(len(labels))]);vary=np.ptp(independent,axis=1)>0;assert labels[vary].equals(node);assert np.max(np.abs(independent[vary]-fq.iloc[:,1:].to_numpy(float)))<1e-12
assert np.array_equal(fq.majorCluster.astype(str).to_numpy(),parents[vary].astype(str).to_numpy());assert np.isfinite(independent).all()
for p in parents.unique():
 ix=np.flatnonzero(parents.eq(p));total=independent[ix].sum(axis=0);assert np.max(np.abs(total-(den[p]>0)))<1e-12
for item in saved['files']:
 p=o/item['name'];assert p.stat().st_size==item['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
out={'status':'PASS_full_cell_source_admission_and_within_parent_frequency_audit','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':int(take.sum()),'n_samples':len(samples),'n_nodes':len(node),'n_sources':len(set(obs['gse_id'][take])),'all_sample_alias_source_status_fresh_scRNA_and_size_rules_recomputed':True,'all_original_cell_ids_labels_and_numerator_counts_exact':True,'all_parent_denominators_and_author_zero_convention_exact':True,'all_constant_row_decisions_exact':True,'all_artifact_SHA_verified':True,'platform_10x_not_newly_certified':True,'new_patient_independence_or_annotation_validation_claim':False};(o/'independent_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
