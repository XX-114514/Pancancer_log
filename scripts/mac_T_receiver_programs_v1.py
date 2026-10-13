"""Fixed external receiver programs across the complete original Mac/T LR family.
Same-data exploration; no new cell identities or receptor-specific causal claim.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.util
import numpy as np,pandas as pd
from scipy.sparse import load_npz
from scipy.linalg import block_diag

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';D=P/'major_LR_TF_design_preflight_v1';O=P/'mac_T_receiver_programs_v1';planp=Path(os.environ['PANCANCER_RECEIVER_PROGRAM_PLAN']);plan=json.loads(planp.read_text());assert not O.exists() and sha(__file__)==plan['producer_sha256'];corep=Path(os.environ['PANCANCER_CR2_CORE']);assert sha(corep)==plan['core_sha256'];spec=importlib.util.spec_from_file_location('core',corep);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
 inv={}
 for name,root,key in [('N',N,'native_inventory_sha256'),('T',T,'TF_inventory_sha256'),('D',D,'design_inventory_sha256')]:assert sha(root/'inventory.json')==plan[key];inv[name]=json.loads((root/'inventory.json').read_text())
 dictionary=Path(os.environ['PANCANCER_PROGRAM_DICTIONARY']);assert sha(dictionary)==plan['dictionary_sha256'];programs=json.loads(dictionary.read_text())['programs'];resourcep=Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']);assert sha(resourcep)==plan['LR_sha256'];resource=pd.read_csv(resourcep,sep='\t');LR=set(g for col in ['ligand','receptor'] for v in resource[col] for g in v.split('_'));counts=AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv';assert sha(counts)==plan['counts_sha256'];totals=pd.read_csv(counts,sep='\t').groupby('sampleID').n_cells.sum();scope=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t').set_index('context_id');source_cache={};data={};started=datetime.now(timezone.utc).isoformat()
 def checked(file,which,root):assert sha(file)==inv[which]['outputs'][str(file.relative_to(root))]
 for cid in plan['contexts']:
  row=scope.loc[cid];sid=row.gse_id
  if sid not in source_cache:
   uf=N/'private_source_major_arrays'/(sid+'_units.tsv');gf=N/'private_source_major_arrays'/(sid+'_axes.npz');mf=N/'private_source_major_arrays'/(sid+'_meanlog.npz');pf=N/'private_source_major_arrays'/(sid+'_pct_positive.npz');tf=T/'private_source_TF_arrays'/(sid+'.npz')
   for f in [uf,gf,mf,pf]:checked(f,'N',N)
   checked(tf,'T',T);u=pd.read_csv(uf,sep='\t');assert np.array_equal(u.unit_row,np.arange(len(u)));lookup={(r.sampleID,r.majorCluster):r.unit_row for r in u.itertuples()};assert len(lookup)==len(u);genes=np.load(gf)['genes'];t=np.load(tf);source_cache[sid]=dict(lookup=lookup,genes={str(g):i for i,g in enumerate(genes)},mean=load_npz(mf),pct=load_npz(pf),TF=t,HSF=t['primary_scores'][:,np.flatnonzero(t['TFs']=='HSF1')[0]])
  af=D/'private_context_arrays'/(cid+'.npz');bf=T/'private_context_arrays'/(cid+'.npz');nf=N/'private_group_arrays'/(row.group_id+'.npz');checked(af,'D',D);checked(bf,'T',T);checked(nf,'N',N);a=np.load(af);b=np.load(bf);native=np.load(nf);mask=native['corrected_material_for_next_inference']==row.material;samples=a['samples'];assert np.array_equal(samples,b['samples']) and np.array_equal(samples,native['samples'][mask]);src=source_cache[sid];si=np.array([src['lookup'][(s,row.source)] for s in samples]);ri=np.array([src['lookup'][(s,row.target)] for s in samples]);sn=native['source_n'][mask];tn=native['target_n'][mask];base=np.column_stack((sn/totals.loc[samples],tn/totals.loc[samples],np.log1p(sn),np.log1p(tn),native['source_mean_detected_nonLR_genes'][mask],native['target_mean_detected_nonLR_genes'][mask]));data[cid]=dict(row=row,a=a,b=b,src=src,si=si,ri=ri,mean=src['mean'][ri],pct=src['pct'][ri],L=native['ligand_meanlog'][mask],Rec=native['receptor_meanlog'][mask],HSF=src['HSF'][si],base=base,samples=samples,blocks=b['blocks'])
 O.mkdir();(O/'private_numeric_fixtures').mkdir();rows=[];sources=[];leave=[];coverage=[];fixture_counts={};fixture_number=0
 variants=['interface_excluded','global_LR_excluded','IRF1_targets_excluded','subtype_mixture']
 for direction,cids in plan['direction_contexts'].items():
  for variant in variants:
   for li in range(4624):
    ligand=str(resource.iloc[li].ligand);receptor=str(resource.iloc[li].receptor);interface=set((ligand+'_'+receptor).split('_'));ids=[c for c in cids if data[c]['a']['mixture_design_status' if variant=='subtype_mixture' else 'primary_design_status'][li]==0]
    result=[dict(direction=direction,variant=variant,LR_index=li,ligand=ligand,receptor=receptor,program=pg['program'],status='fewer_two_original_design_sources',beta=np.nan,se=np.nan,lower=np.nan,upper=np.nan,p_exploratory=np.nan,df=np.nan,n_GSEs=0,n_libraries=0,n_analysis_blocks=0,same_sign_sources=0,source_effects_estimated=0,LOSO_same_sign=False,max_source_information_share=np.nan,score_genes=0) for pg in programs]
    if len({data[c]['row'].gse_id for c in ids})<2:rows.extend(result);continue
    genesets=[]
    for pg in programs:
     key='LR_excluded_genes' if variant=='global_LR_excluded' else 'IRF1_disjoint_genes' if variant=='IRF1_targets_excluded' else 'shared_primary_score_genes';gs=sorted(set(pg[key])-interface);genesets.append(gs)
    valid_indices=[];byProgram={}
    for pi,gs in enumerate(genesets):
     if len(gs)<plan['min_score_genes']:result[pi]['status']='program_genes_insufficient';continue
     ready=[]
     for cid in ids:
      d=data[cid];idx=[d['src']['genes'][g] for g in gs];detect=(d['pct'][:,idx].toarray()>=.01).sum(axis=1);score=np.asarray(d['mean'][:,idx].mean(axis=1)).ravel()
      if (detect>=plan['min_detected_genes']).mean()>=plan['min_library_detected_fraction'] and score.std()>1e-12:ready.append(cid)
     if len({data[c]['row'].gse_id for c in ready})<2:result[pi]['status']='program_detected_coverage_insufficient';continue
     byProgram.setdefault(tuple(ready),[]).append(pi)
    for ids2,pis in byProgram.items():
     ZZ=[];xx=[];yy=[];bb=[];ss=[];names=[]
     for cid in ids2:
      d=data[cid];idx=[d['src']['genes'][g] for g in ligand.split('_')];own=d['mean'][:,idx].toarray().min(axis=1);Z=np.column_stack((np.ones(len(d['samples'])),core.standardize(np.column_stack((d['base'],d['Rec'][:,li],d['HSF'],own)))))
      if variant=='subtype_mixture':Z=np.column_stack((Z,core.standardize(d['a']['source_target_mixture_PCs'])))
      Y=np.column_stack([np.asarray(d['mean'][:,[d['src']['genes'][g] for g in genesets[pi]]].mean(axis=1)).ravel() for pi in pis]);ZZ.append(Z);xx.append(core.standardize(d['L'][:,li]));yy.append(core.standardize(Y));bb.append(d['blocks']);ss.append(np.full(len(d['samples']),d['row'].gse_id));names.extend(d['samples'])
     assert len(set(names))==len(names);Z=block_diag(*ZZ);x=np.concatenate(xx);Y=np.vstack(yy);blocks=np.concatenate(bb);gse=np.concatenate(ss);prep=core.prepare_cr2(Z,x,blocks,min_df=plan['min_df']);effect=core.estimate_cr2(prep,Y) if 'coef_vector' in prep else None;per={};ls={};info={}
     for sid in sorted(set(gse)):
      m=gse==sid;pp=core.prepare_cr2(Z[m],x[m],blocks[m],min_df=0);lp=core.prepare_cr2(Z[~m],x[~m],blocks[~m],min_df=0)
      if 'coef_vector' in pp:per[sid]=core.estimate_cr2(pp,Y[m])
      if 'coef_vector' in lp:ls[sid]=core.estimate_cr2(lp,Y[~m])
      c=prep.get('coef_vector',np.zeros(len(x)));info[sid]=float(c[m]@c[m]/(c@c)) if c@c>0 else 0
      for j,pi in enumerate(pis):
       if sid in per:
        e=per[sid];sources.append(dict(direction=direction,variant=variant,LR_index=li,program=programs[pi]['program'],gse_id=sid,beta=e['beta'][j],se=e['se'][j],lower=e['lower'][j],upper=e['upper'][j],df=pp['df'],information_share=info[sid]))
       if sid in ls:leave.append(dict(direction=direction,variant=variant,LR_index=li,program=programs[pi]['program'],omitted_gse=sid,beta=ls[sid]['beta'][j],se=ls[sid]['se'][j],df=lp['df']))
     for j,pi in enumerate(pis):
      r=result[pi];beta=float(effect['beta'][j]) if effect is not None else np.nan;r.update(status=prep['status'],beta=beta,se=float(effect['se'][j]) if effect is not None else np.nan,lower=float(effect['lower'][j]) if effect is not None else np.nan,upper=float(effect['upper'][j]) if effect is not None else np.nan,p_exploratory=float(effect['p'][j]) if effect is not None and prep['status']=='inference_estimable' else np.nan,df=prep.get('df',np.nan),n_GSEs=len(set(gse)),n_libraries=len(x),n_analysis_blocks=len(set(blocks)),score_genes=len(genesets[pi]),same_sign_sources=sum(np.sign(e['beta'][j])==np.sign(beta) for e in per.values()),source_effects_estimated=len(per),LOSO_same_sign=bool(ls) and all(np.sign(e['beta'][j])==np.sign(beta) for e in ls.values()),max_source_information_share=max(info.values()))
     key=(direction,variant);count=fixture_counts.get(key,0)
     if count<plan['fixture_groups_per_direction_variant'] and prep['status']=='inference_estimable':
      np.savez_compressed(O/'private_numeric_fixtures'/('CASE'+str(fixture_number).zfill(3)+'.npz'),Z=Z,x=x,Y=Y,blocks=blocks,gse=gse,programs=np.array([programs[pi]['program'] for pi in pis]),direction=direction,variant=variant,LR_index=li);fixture_number+=1;fixture_counts[key]=count+1
     if (ligand,receptor) in [tuple(pair) for pair in plan['targeted_numeric_audit_interfaces']] and prep['status']=='inference_estimable':
      group_hash=hashlib.sha256(';'.join(ids2).encode()).hexdigest()[:10];np.savez_compressed(O/'private_numeric_fixtures'/('EXTRA_'+direction+'_'+variant+'_LR'+str(li)+'_'+group_hash+'.npz'),Z=Z,x=x,Y=Y,blocks=blocks,gse=gse,programs=np.array([programs[pi]['program'] for pi in pis]),direction=direction,variant=variant,LR_index=li)
    rows.extend(result)
   print(json.dumps(dict(direction=direction,variant=variant,rows=len(rows))),flush=True)
 out=pd.DataFrame(rows);assert len(out)==4*4*9248
 for variant,g in out.groupby('variant'):
  assert len(g)==36992;q=core.bh(g.p_exploratory.fillna(1));valid=g.status.eq('inference_estimable');out.loc[g.index,'BH_exploratory']=np.where(valid,q,np.nan);out.loc[g.index,'BY_exploratory']=np.where(valid,np.minimum(1,q*np.sum(1/np.arange(1,len(g)+1))),np.nan)
 out.to_csv(O/'all147968_receiver_program_models.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(sources).to_csv(O/'all_source_effects.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(leave).to_csv(O/'all_LOSO_effects.tsv.gz',sep='\t',index=False,na_rep='NA');summary=out.groupby(['variant','program']).agg(models=('LR_index','size'),estimated=('status',lambda s:int(s.eq('inference_estimable').sum())),BH05=('BH_exploratory',lambda s:int((s<.05).sum())),BY05=('BY_exploratory',lambda s:int((s<.05).sum()))).reset_index();summary.to_csv(O/'summary.tsv',sep='\t',index=False);j=dict(status='PRODUCED_PENDING_INDEPENDENT_NUMERIC_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),rows=len(out),family_per_variant=36992,programs_are_external_reference_not_independent_cohort=True,state_labels_assigned=False,summary=summary.to_dict('records'),outputs={str(f.relative_to(O)):sha(f) for f in O.rglob('*') if f.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}),flush=True)
if __name__=='__main__':main()
