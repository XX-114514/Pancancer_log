"""All frozen breast CM/LR mappings, primary material, exact author patient blocks.
Uses frozen full-control-space math; four predeclared full test families.
No global patient/identity change and no independent validation claim.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.util
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';L=P/'CM_edge_LR_extraction_v2';A=P/'CM_LR_full_control_association_v1';M=P/'GSE161529_official_metadata_recon_v2';O=P/'GSE161529_primary_patient_sensitivity_v1';planpath=Path(os.environ['PANCANCER_BREAST_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists()
 for k,f in [('mapping_sha256',M/'private_all69_official_metadata_exact_mapping.tsv'),('metadata_audit_sha256',M/'independent_audit_v1.json'),('association_inventory_sha256',A/'inventory.json'),('association_audit_sha256',A/'independent_audit.json')]:assert sha(f)==plan[k]
 assert sha(L/'inventory.json')==plan['extraction_inventory_sha256'];li=json.loads((L/'inventory.json').read_text());assert sha(A/'all26994_full_control_CM_LR_associations.tsv.gz')==json.loads((A/'inventory.json').read_text())['outputs']['all26994_full_control_CM_LR_associations.tsv.gz']
 modulepath=R/'fullgene_receiver_target_code_v1/CM_LR_full_control_association_v1.py';assert sha(modulepath)==plan['frozen_model_math_sha256'];spec=importlib.util.spec_from_file_location('frozen_full_control_math',modulepath);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 m=pd.read_csv(M/'private_all69_official_metadata_exact_mapping.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID');old=pd.read_csv(A/'all26994_full_control_CM_LR_associations.tsv.gz',sep='\t');base=old[(old.gse_id=='GSE161529')&(old.variant=='endpoint_excluded_library')];assert len(base)==1178 and base.case_id.nunique()==426;rows=[];O.mkdir();(O/'private_resampling_examples').mkdir();start=datetime.now(timezone.utc).isoformat();saved=0
 for gid,g in base.groupby('group_id'):
  file=L/'private_group_arrays'/(gid+'.npz');assert sha(file)==li['outputs'][str(file.relative_to(L))];b=dict(np.load(file,allow_pickle=False));meta=m.loc[b['samples']];primary=meta.official_material_role.eq('tumor_mammary_unresolved_repeat_specimens').to_numpy();assert meta.cell_population[primary].eq('Total').all();index={(l,r):i for i,(l,r) in enumerate(zip(b['ligands'],b['receptors']))};cmidx={c:i for i,c in enumerate(b['CMs'])};cache={}
  for r in g.itertuples():
   j=index[r.ligand,r.receptor];mask=primary&np.isfinite(b['scores'][:,j,0]);nendpoint=int(primary.sum());nobserved=int(mask.sum());blocks=('GSE161529:tumor_author:'+meta.author_patient_token[mask]).to_numpy(dtype=str);Y=np.log1p(b['scores'][mask,j,0]);controls=np.column_stack([b['source_within_major'][mask],b['target_within_major'][mask],np.log1p(b['source_n'][mask]),np.log1p(b['target_n'][mask])]);cat=meta.author_cancer_type[mask];levels=sorted(set(cat));dummy=np.column_stack([cat.eq(k).to_numpy(dtype=float) for k in levels[1:]]) if len(levels)>1 else np.empty((nobserved,0))
   for exposure in ['endpoint_excluded_library','original_H_library']:
    for model in ['required_controls','required_controls_plus_exact_author_cancer_subtype']:
     variant=exposure+'__'+model;activity=b['endpoint_excluded_activity' if exposure=='endpoint_excluded_library' else 'original_H'][mask,cmidx[r.cm]];raw=np.column_stack([activity,controls,dummy]) if model.endswith('subtype') else np.column_stack([activity,controls]);key=(gid,np.packbits(mask).tobytes(),r.cm,variant);info={k:getattr(r,k) for k in ['case_id','group_id','pair_id','gse_id','cancerType','source','target','ligand','receptor','cm','binding_status']};info.update(variant=variant,n_endpoint_primary_libraries=nendpoint,n_observed_primary_libraries=nobserved,n_author_patient_blocks=len(set(blocks)),observed_fraction=nobserved/nendpoint if nendpoint else None,exact_author_cancer_categories=';'.join(levels),original_material_mixed_result_preserved=True,global_registry_unchanged=True,independent_validation=False)
     if not b['native_measured'][j]:result=dict(status='NATIVE_SUBUNIT_NOT_MEASURED')
     elif nendpoint<20 or nobserved<20 or len(set(blocks))<10 or nobserved/nendpoint<.8:result=dict(status='INSUFFICIENT_PRIMARY_PATIENT_COVERAGE')
     else:
      if key not in cache:cache[key]=module.setup(raw,blocks,('breast_primary_sensitivity_v1',)+key)
      result=module.evaluate(cache[key],Y)
     rows.append(info|result)
     if saved<8 and result['status']=='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION':
      np.savez_compressed(O/'private_resampling_examples'/('example'+str(saved+1)+'.npz'),raw_predictors=raw,raw_response=Y,patient_blocks=blocks,case_id=r.case_id,cm=r.cm,variant=variant,resampling_seed=result['resampling_seed']);saved+=1
  print('completed_group',gid,'models',len(rows),flush=True)
 d=pd.DataFrame(rows);assert len(d)==4712
 for v,g in d.groupby('variant'):
  assert len(g)==1178;p=g.residual_block_permutation_p.fillna(1);q=multipletests(p,method='fdr_bh')[1];d.loc[g.index,'full1178_BH_q']=np.where(g.residual_block_permutation_p.notna(),q,np.nan)
 file=O/'all4712_primary_author_patient_models.tsv.gz';d.to_csv(file,sep='\t',index=False,na_rep='NA');summary={v:dict(models=len(g),status_counts=g.status.value_counts().to_dict(),q05=int(g.full1178_BH_q.lt(.05).sum())) for v,g in d.groupby('variant')};inv=dict(status='COMPLETE_SENSITIVITY_PENDING_AUDIT_NOT_INDEPENDENT_VALIDATION',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),summary=summary,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n');print(json.dumps(inv),flush=True)
if __name__=='__main__':main()
