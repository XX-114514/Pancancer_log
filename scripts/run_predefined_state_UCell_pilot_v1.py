"""Official pyUCell scoring of an immutable development-only state reference.

No learned state parameters, no clinical data, no lockbox expression. Individual
arrays stay beside the source inputs. Scores are diagnostic until applicability,
patient/material independence and held-out readout checks are accepted.
"""
import argparse,datetime,hashlib,json,os,sys,time,traceback
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def cp(R,stage,**kw):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=now(),**kw))

def api_checks(uc):
 import numpy as np
 import pandas as pd
 from anndata import AnnData
 from scipy import sparse
 # Separate hand calculation: ranks of a,b in first cell are 1,2 -> score 1;
 # in second cell 3,2 -> score 1-(5-3)/(6-3)=1/3 at rmax=3.
 X=sparse.csr_matrix([[9.,5.,1.,0.],[1.,5.,9.,0.],[0.,0.,0.,0.]])
 a=AnnData(X,var=pd.DataFrame(index=['a','b','c','d']));uc.compute_ucell_scores(a,{'test':['a','b']},max_rank=3,n_jobs=1,chunk_size=1,device='cpu')
 np.testing.assert_allclose(a.obs.test_UCell,[1,1/3,0],atol=1e-6)
 b=AnnData(X,var=pd.DataFrame(index=['a','b','c','d']));uc.compute_ucell_scores(b,{'test':['a','b']},max_rank=3,n_jobs=1,chunk_size=3,device='cpu')
 np.testing.assert_array_equal(a.obs.test_UCell,b.obs.test_UCell)
 return {'n_actual_tests':2,'hand_calculation_and_chunk_invariance':True,'no_smoothing':True}

def main(R):
 start=time.monotonic();cfg=load(R/'config/spec_v1.json');sys.path.insert(0,str(R/'vendor'))
 import numpy as np
 import pandas as pd
 import scipy.sparse as sps
 import anndata
 import pyucell as uc
 import psutil
 assert uc.__version__=='0.8.0';assert psutil.virtual_memory().available>16*1024**3
 cat=load(R/'resources/predefined_state_catalog_v1.json');assert sha(R/'resources/predefined_state_catalog_v1.json')==cfg['catalog_SHA256']
 save(R/'provenance/API_validation_v1.json',dict(completed_UTC=now(),pyucell=uc.__version__,**api_checks(uc)))
 T=Path(cfg['existing_trial']);H=Path(cfg['producer']);tc=load(T/'config/spec_v1.json');hc=load(H/'config/spec_v1.json');scope=load(H/'outputs/scope_and_common_genes_v1.json');genes=scope['genes'];assert len(genes)==cat['n_universe_genes']
 for name in ['private_obs_discovery_cells_v1.tsv.gz','private_obs_other_development_cells_v1.tsv.gz']:
  assert sha(T/'outputs'/name)==cfg['observation_SHA256'][name]
 obs=pd.concat([pd.read_csv(T/'outputs/private_obs_discovery_cells_v1.tsv.gz',sep='\t',dtype={'sample':str,'subtype':str}).assign(role='discovery'),pd.read_csv(T/'outputs/private_obs_other_development_cells_v1.tsv.gz',sep='\t',dtype={'sample':str,'subtype':str}).assign(role='description_only')],ignore_index=True)
 assert len(obs)==776270 and not obs.duplicated(['source','row']).any()
 if cfg.get('smoke_sources'):
  obs=obs[obs.source.isin(cfg['smoke_sources'])].groupby('source',sort=True,group_keys=False).head(cfg['smoke_cells_per_source']).reset_index(drop=True)
 assert len(obs)==cfg['expected_cells']
 split=pd.read_csv(Path(tc['split_run'])/'outputs/source_split.tsv',sep='\t');lock=set(split.source[split.split=='lockbox']);assert not set(obs.source)&lock
 q=pd.DataFrame(tc['qualification_units']);q=q[q.source_qualification=='来源资格通过'];verified={s for s,g in q.groupby('source') if g.patient_block_basis_available.all()};up={(r.source,r.native_expression_unit):str(r.patient_id_source_scoped) for r in q.itertuples()}
 obs['patient_verified']=obs.source.isin(verified);obs['patient_id_source_scoped']=[up[(s,u)] if s in verified else '__unknown__' for s,u in zip(obs.source,obs['sample'])]
 J=Path(hc['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(T/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
 lineage={r['source']:r for r in load(H/'outputs/source_input_lineage_v1.json')}
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');meta=[r for r in cat['programs'] if r['scoreable_TNK_pilot']];sig={r['program_id']:r['common_universe_genes'] for r in meta};names=list(sig)
 np.save(R/'outputs/private_score_row_order_v1.npy',obs.row.to_numpy());obs.to_csv(R/'outputs/private_observations_v1.tsv.gz',sep='\t',index=False)
 save(R/'outputs/score_program_axis_v1.json',dict(program_ids=names,reference_metadata=meta,all_module_eligibility_pending=True))
 index=[];ncomplete=0
 for source,group in obs.groupby('source',sort=True):
  t=time.monotonic();r=lineage[source];names_native,columns=native_measured_genes(native,source);pos=dict(zip(names_native,columns));ph=np.asarray([pos[g] for g in genes]);z=native.store(source);rows=group.row.to_numpy();scorefile=R/'outputs'/('private_'+source+'_UCell_v1.npy');scores=np.lib.format.open_memmap(scorefile,mode='w+',dtype='float32',shape=(len(group),len(sig)));detected=np.empty(len(group),np.int32)
  for lo in range(0,len(group),cfg['read_block_size']):
   hi=min(lo+cfg['read_block_size'],len(group));ix=rows[lo:hi];nr=axes['native_row'][ix];x=gather_csr_rows(z[r['counts_key']],nr,max_span=8192).astype(np.float64).tocsr()
   if 'ambient_removed_counts' in z['layers']:x=(x-gather_nonnegative_float_delta(z['layers/ambient_removed_counts'],nr)).tocsr()
   x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();np.testing.assert_allclose(np.asarray(x.sum(axis=1)).ravel(),totals[ix],rtol=1e-10,atol=1e-7)
   x=x[:,ph].astype(np.float32);detected[lo:hi]=np.diff(x.indptr);a=anndata.AnnData(x,var=pd.DataFrame(index=genes));uc.compute_ucell_scores(a,signatures=sig,max_rank=cfg['max_rank'],missing_genes='skip',ties_method='average',chunk_size=cfg['score_chunk_size'],n_jobs=cfg['threads'],device='cpu');v=a.obs[[n+'_UCell' for n in names]].to_numpy(np.float32);assert np.isfinite(v).all() and v.min()>=-1e-6 and v.max()<=1.000001;scores[lo:hi]=v
  scores.flush();obs_s=group.reset_index(drop=True);obs_s['detected_common_genes']=detected
  obs_s.to_csv(R/'outputs'/('private_'+source+'_axis_v1.tsv.gz'),sep='\t',index=False)
  # Preserve sample/material units. Patient identities are sidecars, never substitute
  # one patient pooled over multiple materials for independent analysis samples.
  value=pd.DataFrame(np.asarray(scores),columns=names);unitkeys=['sample','subtype'];val=value.assign(sample=obs_s['sample'].to_numpy(),subtype=obs_s.subtype.to_numpy());gb=val.groupby(unitkeys,sort=True);count=gb.size();means=gb[names].mean();means.insert(0,'n_cells',count);means.to_csv(R/'outputs'/('private_'+source+'_unit_subtype_means_v1.tsv.gz'),sep='\t')
  for quant in [0.5,0.9,0.95]:
   vq=gb[names].quantile(quant);vq.insert(0,'n_cells',count);vq.to_csv(R/'outputs'/('private_'+source+'_unit_subtype_q'+str(quant).replace('.','')+'_v1.tsv.gz'),sep='\t')
  # source summaries do not replace patient-aware inference
  summ={'source':source,'n_cells':len(group),'n_expression_units':int(group['sample'].nunique()),'n_verified_patients':int(group.patient_id_source_scoped.nunique()) if source in verified else None,'role':group.role.iloc[0],'median_detected_common_genes':float(np.median(detected)),'program_mean':dict(zip(names,np.asarray(scores).mean(axis=0,dtype=np.float64).tolist())),'scores_SHA256':sha(scorefile),'wall_seconds':time.monotonic()-t}
  save(R/'outputs'/(source+'_summary_v1.json'),summ);index.append({k:v for k,v in summ.items() if k!='program_mean'});ncomplete+=len(group);cp(R,'source_scores_complete',source=source,n_complete_cells=ncomplete,n_complete_sources=len(index),n_programs=len(sig),source_wall_seconds=summ['wall_seconds']);native._fields.clear()
  if sum(p.stat().st_size for p in R.rglob('*') if p.is_file())>cfg['output_cap_bytes']:raise RuntimeError('output cap exceeded')
 save(R/'outputs/source_score_index_v1.json',index)
 save(R/'outputs/module_input_contract_v1.json',dict(status='STATE_SCORES_READY_FOR_READOUT_VALIDATION_NOT_MODULE_FIT',cell_axis='source/native_row via immutable row mapping',unit_axis='source/native_expression_unit/subtype; patient/material bridge retained separately',program_axis='frozen resource/program_id/common-universe coverage and applicable lineage',raw_score_method='official_pyucell_0.8.0_no_knn_smoothing',rank_universe_genes=len(genes),max_rank=cfg['max_rank'],global_scores='diagnostic shared cellular-process scores; compare within biological type and source',restricted_scores='diagnostic only until exact parent identity gate is accepted',missing_value_policy='NA_with_reason_not_zero_imputation',remaining_gates=['gene_symbol_alias_and_primary_provenance_review','biological_parent_applicability','patient_material_repeats','source_within_parent_split_half_and_downsampling','reserved_gene_or_orthogonal_validation','redundancy_and_module_input_selection','heldout_projection_under_frozen_reference'],no_new_CM_fit=True))
 save(R/'provenance/terminal_v1.json',dict(status='OFFICIAL_PREDEFINED_UCELL_DEVELOPMENT_SCORES_COMPLETE_REQUIRES_READOUT_VALIDATION',completed_UTC=now(),n_cells=ncomplete,n_sources=len(index),n_programs=len(sig),lockbox_expression_rows_read=0,wall_seconds=time.monotonic()-start,catalog_SHA256=cfg['catalog_SHA256'],index_SHA256=sha(R/'outputs/source_score_index_v1.json'),patient_rows_exported=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();R=Path(a.output)
 try:main(R)
 except BaseException:
  save(R/'provenance/failure_v1.json',dict(failed_UTC=now(),traceback=traceback.format_exc()));raise
