"""Finite post-result IFNG exploration on all three frozen T/Mac CM edges.
Sample matrices stay private. Pseudobulk gene ranks are not UCell scores.
"""
import os,json,hashlib,importlib.util,gc
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.sparse import load_npz
from scipy.linalg import block_diag
from scipy.stats import rankdata

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);C=Path(os.environ['PANCANCER_CODE_ROOT']);P=R/'immune_same_name_family_source_review_v1';U=AR/'20261010_uniform_expression_annotation_v1';B=AR/'20261010_expression_annotation_downstream_v1';N=P/'major_pair_native_assembly_v1';D=P/'major_LR_TF_design_preflight_v1';O=P/'CM_fine_IFNG_finite_round_v1';pp=C/'CM_fine_IFNG_finite_round_plan_v1.json';plan=json.loads(pp.read_text());assert not O.exists() and sha(__file__)==plan['producer_sha256'];start=datetime.now(timezone.utc).isoformat()
 spec=importlib.util.spec_from_file_location('core',C/'cluster_CR2_core_v1.py');core=importlib.util.module_from_spec(spec);assert sha(C/'cluster_CR2_core_v1.py')==plan['core_sha256'];spec.loader.exec_module(core)
 assert sha(N/'inventory.json')==plan['native_inventory_sha256'] and sha(D/'inventory.json')==plan['design_inventory_sha256'];ni=json.loads((N/'inventory.json').read_text());di=json.loads((D/'inventory.json').read_text())
 for name,h in plan['native_cache']['root_files'].items():assert sha(U/name)==h
 cp=B/'inputs/sample_major_subcluster_counts.tsv';catp=B/'downstream_review_v2/number_to_annotation_dictionary.tsv';assert sha(cp)==plan['counts_sha256'] and sha(catp)==plan['current_input_sha256']['catalogue'];counts=pd.read_csv(cp,sep='\t');maj=counts.groupby(['sampleID','majorCluster']).n_cells.sum();catalog=pd.read_csv(catp,sep='\t').set_index('node_id');original=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t')
 dp=C/'program_reference_coverage_v1.json';assert sha(dp)==plan['dictionary_sha256'];programs=json.loads(dp.read_text())['programs'];assert len(programs)==4
 lr=C/'frozen_LIANA171_consensus4624_v1.tsv';assert sha(lr)==plan['LR_sha256'];resource=pd.read_csv(lr,sep='\t');lrgenes=set(g for v in resource.to_numpy().ravel() for g in v.split('_'));scope=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t').set_index('context_id')
 O.mkdir();(O/'private_fixtures').mkdir();data={};source_coverage=[]
 wanted=set(v for pair in plan['pairs'] for v in [pair['sender'],pair['receiver']])
 for sid in sorted(set(scope.loc[plan['contexts'],'gse_id'])):
  root=U/'by_gse'/sid;assert sha(root/'manifest.json')==plan['native_cache']['sources'][sid]['manifest_sha256'];manifest=json.loads((root/'manifest.json').read_text())
  for name,h in manifest['output_sha256'].items():assert sha(root/name)==h
  units=pd.read_csv(root/'units.tsv',sep='\t');assert np.array_equal(units.unit_row,np.arange(len(units)));units['node_id']=original.node_id.to_numpy()[units.node_index];units['subCluster']=units.node_id.map(catalog.subCluster);units['majorCluster']=units.node_id.map(catalog.majorCluster);take=units[units.candidate_eligible & units.subCluster.isin(wanted)].copy();lookup={(r.sample_uid,r.subCluster):i for i,r in enumerate(take.itertuples())};assert len(lookup)==len(take)
  expected=counts.set_index(['sampleID','subCluster']).n_cells;assert all(int(expected.loc[s,c])==int(take.iloc[i].n_cells) for (s,c),i in lookup.items())
  genes=pd.read_csv(root/'genes.tsv',sep='\t').gene.to_numpy();assert len(set(genes))==len(genes);gi={g:i for i,g in enumerate(genes)};logs=load_npz(root/'sum_log1p_CP10K.npz')[take.unit_row].tocsr();positive=load_npz(root/'positive_counts.npz')[take.unit_row].tocsr();nc=take.n_cells.to_numpy();mean=logs.multiply(1/nc[:,None]).tocsr();pct=positive.multiply(1/nc[:,None]).tocsr();nonlr=~pd.Index(genes).isin(lrgenes);complexity=np.asarray(positive[:,nonlr].sum(1)).ravel()/nc;dense=mean.toarray();ranks=rankdata(dense,axis=1,method='average')/len(genes)
  program_scores={};program_detect={}
  for pg in programs:
   gs=sorted(set(pg['shared_primary_score_genes'])-{'IFNG','IFNGR1','IFNGR2'});assert set(gs).issubset(gi);ix=[gi[g] for g in gs];program_scores[pg['program']]={'meanlog':dense[:,ix].mean(1),'pseudobulk_gene_rank':ranks[:,ix].mean(1)};program_detect[pg['program']]=np.sum(pct[:,ix].toarray()>=.01,axis=1);source_coverage.append(dict(gse_id=sid,program=pg['program'],score_genes=len(gs),native_genes=len(genes),selected_fine_units=len(take)))
  np.savez_compressed(O/('private_score_example_'+sid+'.npz'),unit_row=int(take.iloc[0].unit_row),n_cells=nc[0],programs=np.array([pg['program'] for pg in programs]),meanlog=np.array([program_scores[pg['program']]['meanlog'][0] for pg in programs]),rank=np.array([program_scores[pg['program']]['pseudobulk_gene_rank'][0] for pg in programs]),detected=np.array([program_detect[pg['program']][0] for pg in programs]))
  data[sid]=dict(units=take,lookup=lookup,nc=nc,complexity=complexity,IFNG=dense[:,gi['IFNG']],R=dense[:,[gi['IFNGR1'],gi['IFNGR2']]].min(1),scores=program_scores,detected=program_detect,Lpct=pct[:,gi['IFNG']].toarray().ravel(),Rpct=pct[:,[gi['IFNGR1'],gi['IFNGR2']]].toarray().min(1));del dense,ranks,logs,positive,mean,pct;gc.collect()
 rows=[];coverage=[];source_effects=[];loso=[];fixture=0
 for pair in plan['pairs']:
  pairdata={}
  for cid in plan['contexts']:
   r=scope.loc[cid];sid=r.gse_id;src=data[sid];af=D/'private_context_arrays'/(cid+'.npz');nf=N/'private_group_arrays'/(r.group_id+'.npz');assert sha(af)==di['outputs'][str(af.relative_to(D))] and sha(nf)==ni['outputs'][str(nf.relative_to(N))];a=np.load(af);n=np.load(nf);mask=n['corrected_material_for_next_inference']==r.material;assert np.array_equal(a['samples'],n['samples'][mask]);samples=a['samples'];blocks=n['corrected_blocks_for_next_inference'][mask];reviewed=n['corrected_identity_reviewed_for_next_inference'][mask]
   present=np.array([(s,pair['sender']) in src['lookup'] and (s,pair['receiver']) in src['lookup'] for s in samples]);samp=samples[present];si=np.array([src['lookup'][s,pair['sender']] for s in samp],dtype=int);ri=np.array([src['lookup'][s,pair['receiver']] for s in samp],dtype=int);enough=(src['nc'][si]>=20)&(src['nc'][ri]>=20);si=si[enough];ri=ri[enough];samp=samp[enough];bb=blocks[present][enough];rv=reviewed[present][enough];ok=len(samp)>=20 and len(set(bb))>=10
   q=dict(pair_id=pair['pair_id'],CM=pair['CM'],context_id=cid,gse_id=sid,cancer_label=r.cancerType,material=r.material,sender=pair['sender'],receiver=pair['receiver'],n_parent_libraries=len(samples),n_fine_both_captured=int(present.sum()),n_fine_both_ge20=len(samp),n_blocks=len(set(bb)),n_reviewed_blocks=len(set(bb[rv])),coverage_ready=ok);coverage.append(q)
   if not ok:continue
   ms=np.array([maj.loc[s,'T cell'] for s in samp]);mr=np.array([maj.loc[s,'Macrophage'] for s in samp]);base=np.column_stack((src['nc'][si]/ms,src['nc'][ri]/mr,np.log1p(src['nc'][si]),np.log1p(src['nc'][ri]),src['complexity'][si],src['complexity'][ri],src['R'][ri]));pairdata[cid]=dict(sid=sid,base=base,x=src['IFNG'][si],own=src['IFNG'][ri],si=si,ri=ri,src=src,blocks=bb,n=len(samp),RNA_ready=(src['Lpct'][si]>=.01).mean()>=.5 and (src['Rpct'][ri]>=.01).mean()>=.8)
  for score in plan['score_variants']:
   for control in plan['control_variants']:
    for pg in programs:
     name=pg['program'];ready={}
     for cid,d in pairdata.items():
      s=d['src'];y=s['scores'][name][score][d['ri']];sendergamma=s['scores'][programs[0]['program']][score][d['si']];det=s['detected'][name][d['ri']];cols=d['base'] if control==plan['control_variants'][0] else np.column_stack((d['base'],sendergamma,d['own']));Z=np.column_stack((np.ones(d['n']),core.standardize(cols)));x=core.standardize(d['x']);Y=core.standardize(y)[:,None];X=np.column_stack((Z,x));reason='ready'
      if not d['RNA_ready']:reason='endpoint_RNA_detection_insufficient'
      elif (det>=10).mean()<.8:reason='receiver_program_detection_insufficient'
      elif np.std(y)<=1e-12:reason='receiver_program_constant'
      elif np.linalg.matrix_rank(X,tol=1e-10)<X.shape[1]:reason='design_collinear'
      elif np.linalg.cond(X)>1e8:reason='design_condition_insufficient'
      if reason=='ready':ready[cid]=dict(Z=Z,x=x,Y=Y,blocks=d['blocks'],sid=d['sid'])
      d.setdefault('reasons',{})[(score,control,name)]=reason
     for test_id in plan['contexts']+['POOLED']:
      row=dict(pair_id=pair['pair_id'],CM=pair['CM'],sender=pair['sender'],receiver=pair['receiver'],score=score,control=control,program=name,test_id=test_id,status='fine_coverage_insufficient',beta=np.nan,se=np.nan,df=np.nan,p_exploratory=np.nan,lower=np.nan,upper=np.nan,n_libraries=0,n_blocks=0,n_GSEs=0)
      ids=[test_id] if test_id in ready else []
      if test_id=='POOLED':
       ids=list(ready);row['status']='fewer_two_estimable_sources'
       if len({ready[c]['sid'] for c in ids})<2:ids=[]
      elif test_id in pairdata:row['status']=pairdata[test_id]['reasons'][(score,control,name)]
      if ids:
       Z=block_diag(*[ready[c]['Z'] for c in ids]);x=np.concatenate([ready[c]['x'] for c in ids]);Y=np.vstack([ready[c]['Y'] for c in ids]);bb=np.concatenate([ready[c]['blocks'] for c in ids]);gse=np.concatenate([np.full(len(ready[c]['x']),ready[c]['sid']) for c in ids]);prep=core.prepare_cr2(Z,x,bb,min_df=4);row.update(status=prep['status'],n_libraries=len(x),n_blocks=len(set(bb)),n_GSEs=len(set(gse)))
       if 'coef_vector' in prep:
        e=core.estimate_cr2(prep,Y);row.update({k:float(e[k][0]) for k in ['beta','se','lower','upper']});row['df']=prep['df'];row['p_exploratory']=float(e['p'][0]) if prep['status']=='inference_estimable' else np.nan
        np.savez_compressed(O/'private_fixtures'/('CASE'+str(fixture).zfill(3)+'.npz'),Z=Z,x=x,Y=Y,blocks=bb,gse=gse,row_index=len(rows));fixture+=1
        if test_id=='POOLED':
         c=prep['coef_vector'];infos=[];signs=[];leaves=[]
         for sid in sorted(set(gse)):
          m=gse==sid;info=float(c[m]@c[m]/(c@c));infos.append(info)
          for mask_,table,label in [(m,source_effects,'gse_id'),(~m,loso,'omitted_gse')]:
           zz=Z[mask_];keep=np.any(abs(zz)>1e-12,axis=0);ep=core.prepare_cr2(zz[:,keep],x[mask_],bb[mask_],min_df=0)
           if 'coef_vector' in ep:
            ee=core.estimate_cr2(ep,Y[mask_]);b=float(ee['beta'][0]);table.append({k:row[k] for k in ['pair_id','score','control','program']}|{label:sid,'beta':b,'df':ep['df'],'information_share':info});(signs if label=='gse_id' else leaves).append(b)
         row.update(same_sign_sources=sum(np.sign(b)==np.sign(row['beta']) for b in signs),LOSO_same_sign=bool(leaves) and all(np.sign(b)==np.sign(row['beta']) for b in leaves),max_source_information_share=max(infos))
      rows.append(row)
 out=pd.DataFrame(rows);assert len(out)==336
 for _,g in out.groupby(['score','control']):
  assert len(g)==84;q=core.bh(g.p_exploratory.fillna(1));valid=g.status.eq('inference_estimable');out.loc[g.index,'BH_exploratory']=np.where(valid,q,np.nan);out.loc[g.index,'BY_exploratory']=np.where(valid,np.minimum(1,q*sum(1/np.arange(1,85))),np.nan)
 out.to_csv(O/'all336_fine_IFNG_models.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(coverage).to_csv(O/'all18_fine_context_coverage.tsv',sep='\t',index=False);pd.DataFrame(source_coverage).to_csv(O/'program_gene_coverage.tsv',sep='\t',index=False);pd.DataFrame(source_effects).to_csv(O/'pooled_source_effects.tsv',sep='\t',index=False);pd.DataFrame(loso).to_csv(O/'pooled_LOSO_effects.tsv',sep='\t',index=False)
 summary=out.groupby(['score','control']).agg(models=('status','size'),estimated=('status',lambda s:int(s.eq('inference_estimable').sum())),BH05=('BH_exploratory',lambda s:int((s<.05).sum())),BY05=('BY_exploratory',lambda s:int((s<.05).sum()))).reset_index();summary.to_csv(O/'summary.tsv',sep='\t',index=False)
 j=dict(status='PRODUCED_PENDING_NUMERIC_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(pp),models=len(out),fixtures=fixture,summary=summary.to_dict('records'),outputs={str(f.relative_to(O)):sha(f) for f in O.rglob('*') if f.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}),flush=True)
if __name__=='__main__':main()
