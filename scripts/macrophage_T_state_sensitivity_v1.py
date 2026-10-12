"""Coverage-selected Mac<->T primary-tumor bridge, all LR/TFs, exploratory.
Context-specific full nuisance blocks; shared reviewed CRC patients remain clustered.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import os,json,hashlib,importlib.util
import numpy as np,pandas as pd
from scipy.linalg import block_diag
from scipy.sparse import load_npz
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';D=P/'major_LR_TF_design_preflight_v1';O=P/'macrophage_T_state_sensitivity_v1';planp=Path(os.environ['PANCANCER_MAC_T_STATE_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists();corep=Path(os.environ['PANCANCER_CR2_CORE']);assert sha(corep)==plan['core_sha256'];spec=importlib.util.spec_from_file_location('core',corep);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
 for f,h in [(N/'inventory.json',plan['native_inventory_sha256']),(T/'inventory.json',plan['TF_inventory_sha256']),(D/'inventory.json',plan['design_inventory_sha256']),(D/'independent_audit.json',plan['design_audit_sha256'])]:assert sha(f)==h
 ni=json.loads((N/'inventory.json').read_text());ti=json.loads((T/'inventory.json').read_text());di=json.loads((D/'inventory.json').read_text());counts_path=AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv';assert sha(counts_path)==plan['counts_sha256'];total=pd.read_csv(counts_path,sep='\t').groupby('sampleID').n_cells.sum();resourcep=Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']);assert sha(resourcep)==plan['LR_sha256'];resource=pd.read_csv(resourcep,sep='\t');scope=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t',keep_default_na=False).set_index('context_id');data={};started=datetime.now(timezone.utc).isoformat()
 for cid in plan['contexts']:
  row=scope.loc[cid];af=D/'private_context_arrays'/(cid+'.npz');bf=T/'private_context_arrays'/(cid+'.npz');nf=N/'private_group_arrays'/(row.group_id+'.npz')
  for f,i,root in [(af,di,D),(bf,ti,T),(nf,ni,N)]:assert sha(f)==i['outputs'][str(f.relative_to(root))]
  a=np.load(af);b=np.load(bf);nat=np.load(nf);mask=nat['corrected_material_for_next_inference']==row.material;samples=a['samples'];assert np.array_equal(samples,b['samples']) and np.array_equal(samples,nat['samples'][mask]);sn=nat['source_n'][mask];tn=nat['target_n'][mask];base=np.column_stack((sn/total.loc[samples],tn/total.loc[samples],np.log1p(sn),np.log1p(tn),nat['source_mean_detected_nonLR_genes'][mask],nat['target_mean_detected_nonLR_genes'][mask]));data[cid]=dict(row=row,a=a,b=b,L=nat['ligand_meanlog'][mask],Rec=nat['receptor_meanlog'][mask],base=base,blocks=b['blocks'],samples=samples)
 # Predeclared full-family sensitivity: sender HSF1 program (LR genes excluded)
 # and receiver-intrinsic same ligand RNA. These are covariates, not state labels.
 source_cache={}
 for cid,d in data.items():
  sid=d['row'].gse_id
  if sid not in source_cache:
   uf=N/'private_source_major_arrays'/(sid+'_units.tsv');gf=N/'private_source_major_arrays'/(sid+'_axes.npz');mf=N/'private_source_major_arrays'/(sid+'_meanlog.npz');tf=T/'private_source_TF_arrays'/(sid+'.npz')
   for f,inv,root in [(uf,ni,N),(gf,ni,N),(mf,ni,N),(tf,ti,T)]:assert sha(f)==inv['outputs'][str(f.relative_to(root))]
   u=pd.read_csv(uf,sep='\t');assert np.array_equal(u.unit_row,np.arange(len(u)));axis=np.load(gf);scores=np.load(tf);hix=np.flatnonzero(scores['TFs']=='HSF1');assert len(hix)==1 and scores['primary_target_counts'][hix[0]]>=20;lookup={(r.sampleID,r.majorCluster):int(r.unit_row) for r in u.itertuples()};assert len(lookup)==len(u);source_cache[sid]=(lookup,{str(g):i for i,g in enumerate(axis['genes'])},load_npz(mf),scores['primary_scores'][:,hix[0]])
  lookup,genes,mean,hscore=source_cache[sid];si=np.array([lookup[(sample,d['row'].source)] for sample in d['samples']]);ri=np.array([lookup[(sample,d['row'].target)] for sample in d['samples']]);d['sender_HSF1']=hscore[si];assert np.isfinite(d['sender_HSF1']).all();d['receiver_ligand']=np.full(d['L'].shape,np.nan)
  for li,ligand in enumerate(resource.ligand):
   components=str(ligand).split('_')
   if all(g in genes for g in components):d['receiver_ligand'][:,li]=mean[ri][:,[genes[g] for g in components]].toarray().min(axis=1)
 O.mkdir();(O/'private_pooled_input_cases').mkdir();rows=[];component=[];sources=[];leave=[];coverage=[];saved=0;fixture_counts=defaultdict(int)
 for direction,cids in plan['direction_contexts'].items():
  for variant in ['abundance_complexity','subtype_mixture','full_gene_sensitivity']:
   for li in range(4624):
    byTF=defaultdict(list)
    for cid in cids:
     d=data[cid];a=d['a'];ds=a['mixture_design_status' if variant=='subtype_mixture' else 'primary_design_status'][li]
     if ds!=0:continue
     for tf in a['candidate_TF_index'][a['candidate_LR_index']==li]:byTF[int(tf)].append(cid)
    patterns=defaultdict(list)
    for tf,ids in byTF.items():
     if len({data[c]['row'].gse_id for c in ids})>=2:patterns[tuple(ids)].append(tf)
    result=dict(direction=direction,variant=variant,LR_index=li,ligand=resource.iloc[li].ligand,receptor=resource.iloc[li].receptor,eligible_TFs=len([v for vv in patterns.values() for v in vv]),status='coverage_fewer2sources',joint_p_exploratory=np.nan,lead_TF='',lead_beta=np.nan,lead_lower=np.nan,lead_upper=np.nan,lead_df=np.nan,lead_sources=np.nan,lead_same_sign_sources=np.nan,lead_LOSO_same_sign=False,lead_max_source_information_share=np.nan,family_size=9248,independent_validation=False);allTF=[]
    for ids,TFlist in patterns.items():
     blocks=np.concatenate([data[c]['blocks'] for c in ids]);samples=np.concatenate([data[c]['samples'] for c in ids]);assert len(set(samples))==len(samples);gse=np.concatenate([np.full(len(data[c]['samples']),data[c]['row'].gse_id) for c in ids]);ZZ=[];xx=[];yy=[]
     for cid in ids:
      d=data[cid];n=len(d['samples']);assert np.isfinite(d['receiver_ligand'][:,li]).all();Z=np.column_stack((np.ones(n),core.standardize(np.column_stack((d['base'],d['Rec'][:,li],d['sender_HSF1'],d['receiver_ligand'][:,li])))))
      if variant=='subtype_mixture':Z=np.column_stack((Z,core.standardize(d['a']['source_target_mixture_PCs'])))
      ZZ.append(Z);xx.append(core.standardize(d['L'][:,li]));yy.append(core.standardize(d['b']['sensitivity_TF_scores' if variant=='full_gene_sensitivity' else 'primary_TF_scores'][:,TFlist]))
     Z=block_diag(*ZZ);x=np.concatenate(xx);Y=np.vstack(yy);prep=core.prepare_cr2(Z,x,blocks,min_df=plan['min_df']);effect=core.estimate_cr2(prep,Y) if 'coef_vector' in prep else None;per=[];loso=[];information={};xr=prep['coef_vector'] if 'coef_vector' in prep else np.zeros(len(x))
     for sid in sorted(set(gse)):
      mask=gse==sid;info=float((xr[mask]@xr[mask])/(xr@xr)) if xr@xr>0 else 0;information[sid]=info;pp=core.prepare_cr2(Z[mask],x[mask],blocks[mask],min_df=0)
      if 'coef_vector' in pp:
       ee=core.estimate_cr2(pp,Y[mask]);per.append((sid,ee));
       for j,tf in enumerate(TFlist):sources.append(dict(direction=direction,variant=variant,LR_index=li,TF=str(data[ids[0]]['b']['TFs'][tf]),gse_id=sid,beta=float(ee['beta'][j]),se=float(ee['se'][j]),df=float(pp['df']),lower=float(ee['lower'][j]),upper=float(ee['upper'][j]),n_libraries=int(mask.sum()),n_patient_or_sample_blocks=len(set(blocks[mask])),information_share=info))
      keep=~mask;lp=core.prepare_cr2(Z[keep],x[keep],blocks[keep],min_df=0)
      if 'coef_vector' in lp:
       le=core.estimate_cr2(lp,Y[keep]);loso.append((sid,le))
       for j,tf in enumerate(TFlist):leave.append(dict(direction=direction,variant=variant,LR_index=li,TF=str(data[ids[0]]['b']['TFs'][tf]),omitted_gse=sid,beta=float(le['beta'][j]),se=float(le['se'][j]),df=float(lp['df'])))
     for j,tf in enumerate(TFlist):
      name=str(data[ids[0]]['b']['TFs'][tf]);p=float(effect['p'][j]) if effect is not None and prep['status']=='inference_estimable' else np.nan;beta=float(effect['beta'][j]) if effect is not None else np.nan;item=dict(direction=direction,variant=variant,LR_index=li,TF=name,status=prep['status'],beta=beta,se=float(effect['se'][j]) if effect is not None else np.nan,p_exploratory=p,df=prep.get('df',np.nan),lower=float(effect['lower'][j]) if effect is not None else np.nan,upper=float(effect['upper'][j]) if effect is not None else np.nan,n_libraries=len(samples),n_analysis_blocks=len(set(blocks)),n_contexts=len(ids),n_GSEs=len(set(gse)),same_sign_source_effects=sum(np.sign(ee['beta'][j])==np.sign(beta) for sid,ee in per),source_effects_estimated=len(per),all_LOSO_same_sign=bool(loso) and all(np.sign(ee['beta'][j])==np.sign(beta) for sid,ee in loso),max_source_information_share=max(information.values()),contexts=';'.join(ids));component.append(item);allTF.append(item)
     coverage.append(dict(direction=direction,variant=variant,LR_index=li,TFs=len(TFlist),contexts=';'.join(ids),n_GSEs=len(set(gse)),n_libraries=len(samples),blocks=len(set(blocks))))
     if fixture_counts[(direction,variant)]<4 and prep['status']=='inference_estimable':
      np.savez_compressed(O/'private_pooled_input_cases'/('POOL'+str(saved).zfill(3)+'.npz'),Z=Z,x=x,Y=Y,blocks=blocks,gse=gse,TFs=np.array([str(data[ids[0]]['b']['TFs'][tf]) for tf in TFlist]),LR_index=li,direction=direction,variant=variant);saved+=1;fixture_counts[(direction,variant)]+=1
    if allTF:
     valid=[a for a in allTF if np.isfinite(a['p_exploratory'])]
     if valid:
      leader=min(valid,key=lambda a:a['p_exploratory']);result.update(status='pooled_exploratory_estimated',joint_p_exploratory=min(1,len(allTF)*leader['p_exploratory']),lead_TF=leader['TF'],lead_beta=leader['beta'],lead_lower=leader['lower'],lead_upper=leader['upper'],lead_df=leader['df'],lead_sources=leader['n_GSEs'],lead_same_sign_sources=leader['same_sign_source_effects'],lead_LOSO_same_sign=leader['all_LOSO_same_sign'],lead_max_source_information_share=leader['max_source_information_share'])
     else:result['status']='pooled_inference_insufficient'
    rows.append(result)
   print(json.dumps(dict(direction=direction,variant=variant,LR_rows=len(rows),TF_components=len(component))),flush=True)
 out=pd.DataFrame(rows);assert len(out)==3*9248
 for v,g in out.groupby('variant'):
  assert len(g)==9248;p=g.joint_p_exploratory.fillna(1).to_numpy();q=core.bh(p);valid=g.status.eq('pooled_exploratory_estimated').to_numpy();out.loc[g.index,'joint_BH_exploratory']=np.where(valid,q,np.nan);out.loc[g.index,'joint_BY_exploratory']=np.where(valid,np.minimum(1,q*sum(1/i for i in range(1,9249))),np.nan)
 out.to_csv(O/'all27744_pooled_LR_models.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(component).to_csv(O/'all_pooled_TF_components.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(sources).to_csv(O/'all_source_TF_effects.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(leave).to_csv(O/'all_LOSO_TF_effects.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(coverage).to_csv(O/'all_pattern_coverage.tsv',sep='\t',index=False);summary=out.groupby('variant').agg(models=('LR_index','size'),estimated=('status',lambda v:int((v=='pooled_exploratory_estimated').sum())),BH05=('joint_BH_exploratory',lambda v:int((v<.05).sum())),BY05=('joint_BY_exploratory',lambda v:int((v<.05).sum()))).reset_index();summary.to_csv(O/'summary.tsv',sep='\t',index=False);j=dict(status='MULTISOURCE_EXPLORATORY_POOL_COMPLETE_PENDING_NUMERIC_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),LR_rows=len(out),TF_components=len(component),summary=summary.to_dict('records'),no_causal_binding_or_independent_validation=True,intrinsic_program_controls=['sender_HSF1_LR_excluded','receiver_same_ligand_native_meanlog'],state_labels_assigned=False,outputs={str(f.relative_to(O)):sha(f) for f in O.rglob('*') if f.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}))
if __name__=='__main__':main()
