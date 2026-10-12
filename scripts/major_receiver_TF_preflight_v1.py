"""All frozen major units: signed CollecTRI ULM expression footprints.
Preflight only. No sample vectors exported; no ligand-response p-values.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd
from scipy.sparse import load_npz,csr_matrix

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()

def ulm_sparse(X,W):
 """Centered correlation, then the official decoupler1.9.2 t transform."""
 n=X.shape[1];sx=np.asarray(X.sum(1)).ravel();sy=np.asarray(W.sum(0)).ravel()
 xx=np.asarray(X.multiply(X).sum(1)).ravel()-sx*sx/n
 yy=np.asarray(W.multiply(W).sum(0)).ravel()-sy*sy/n
 cross=np.asarray((X@W).toarray())-sx[:,None]*sy[None,:]/n
 r=cross/np.sqrt(xx[:,None]*yy[None,:]);r=np.clip(r,-1,1)
 return r*np.sqrt((n-2)/((1-r+1e-16)*(1+r+1e-16)))

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';O=P/'major_receiver_TF_preflight_v1';D=Path(os.environ['PANCANCER_TF_REFERENCE']);planpath=Path(os.environ['PANCANCER_TF_PLAN']);plan=json.loads(planpath.read_text());assert not O.exists() and sha(__file__)==plan['producer_sha256'];assert sha(N/'inventory.json')==plan['major_inventory_sha256'];assert sha(N/'independent_audit.json')==plan['major_audit_sha256']
 inv=json.loads((N/'inventory.json').read_text());audit=json.loads((N/'independent_audit.json').read_text());assert audit['case_records']==1276224 and audit['status'].startswith('PASS_')
 for n,h in plan['reference_files'].items():assert sha(D/n)==h
 raw=pd.read_csv(D/'human_collectri_raw.tsv',sep='\t');net=pd.read_csv(D/'human_collectri_signed_complexes.tsv',sep='\t');assert len(net)==40634 and net.source.nunique()==1163 and not net.duplicated(['source','target']).any()
 resource=pd.read_csv(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']),sep='\t');assert sha(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']))==plan['LR_resource_sha256'];lrgenes={g for row in resource.itertuples(index=False,name=None) for v in row for g in v.split('_')};TFs=np.sort(net.source.unique());tc={v:i for i,v in enumerate(TFs)};total_targets=net.groupby('source').size().to_dict()
 O.mkdir();(O/'private_source_TF_arrays').mkdir();(O/'private_context_arrays').mkdir();source_summary=[];coverage_rows=[];sources={};started=datetime.now(timezone.utc).isoformat()
 for sid in sorted(r['gse_id'] for r in inv['summary']):
  for suf in ['units.tsv','meanlog.npz','pct_positive.npz','axes.npz']:
   f=N/'private_source_major_arrays'/(sid+'_'+suf);assert sha(f)==inv['outputs'][str(f.relative_to(N))]
  units=pd.read_csv(N/'private_source_major_arrays'/(sid+'_units.tsv'),sep='\t',keep_default_na=False);axis=np.load(N/'private_source_major_arrays'/(sid+'_axes.npz'));genes=axis['genes'];X=load_npz(N/'private_source_major_arrays'/(sid+'_meanlog.npz')).astype(float);pct=load_npz(N/'private_source_major_arrays'/(sid+'_pct_positive.npz')).astype(float);assert X.shape==(len(units),len(genes)) and len(set(genes))==len(genes)
  # Official extraction removes genes empty across supplied units; never structural NA->zero.
  active=np.asarray(X.getnnz(0)).ravel()>0;all_scores=[];all_counts=[];all_detected=[];universes=[]
  for variant,keep in [('nonLR_primary',active & np.array([g not in lrgenes for g in genes])),('full_gene_sensitivity',active)]:
   ids=np.where(keep)[0];gg=genes[ids];cc={g:i for i,g in enumerate(gg)};nn=net[net.target.isin(cc)];wi=np.array([cc[g] for g in nn.target]);wj=np.array([tc[g] for g in nn.source]);W=csr_matrix((nn.weight.to_numpy(float),(wi,wj)),shape=(len(gg),len(TFs)));counts=np.asarray((W!=0).sum(0)).ravel().astype(int);eligible=counts>=5
   Y=ulm_sparse(X[:,ids],W[:,eligible]);scores=np.full((len(units),len(TFs)),np.nan);scores[:,eligible]=Y
   # >=1% cells detecting a target; all targets weighted by gene identity, not sign here.
   detected=(pct[:,ids]>=.01).astype(np.int32)@(W!=0).astype(np.int32);detected=detected.toarray()
   all_scores.append(scores);all_counts.append(counts);all_detected.append(detected);universes.append(ids)
   for j,tf in enumerate(TFs):coverage_rows.append(dict(gse_id=sid,variant=variant,TF=tf,all_resource_targets=int(total_targets[tf]),measured_nonempty_targets=int(counts[j]),min5_scored=bool(eligible[j]),min20_target_gate=bool(counts[j]>=20),selected_major_units=len(units),gene_universe=len(ids),median_detected_targets_ge1pct=float(np.median(detected[:,j]))))
  np.savez_compressed(O/'private_source_TF_arrays'/(sid+'.npz'),TFs=TFs,primary_scores=all_scores[0],sensitivity_scores=all_scores[1],primary_target_counts=all_counts[0],sensitivity_target_counts=all_counts[1],primary_detected_target_counts=all_detected[0],sensitivity_detected_target_counts=all_detected[1],primary_gene_indices=universes[0],sensitivity_gene_indices=universes[1])
  sources[sid]=(units,all_scores,all_counts,all_detected);source_summary.append(dict(gse_id=sid,major_units=len(units),native_genes=len(genes),nonempty_genes=int(active.sum()),nonLR_universe=len(universes[0]),TFs_primary_ge5=int((all_counts[0]>=5).sum()),TFs_primary_ge20=int((all_counts[0]>=20).sum()),TFs_full_ge20=int((all_counts[1]>=20).sum())));print(json.dumps(source_summary[-1]),flush=True);del X,pct;gc.collect()
 groups=pd.read_csv(N/'all276_major_direction_contexts.tsv',sep='\t',keep_default_na=False);contexts=[];tfrows=[];lrrows=[]
 for base in groups.to_dict('records'):
  sid=base['gse_id'];f=N/'private_group_arrays'/(base['group_id']+'.npz');assert sha(f)==inv['outputs'][str(f.relative_to(N))];a=np.load(f);units,ss,cc,dd=sources[sid];index={(r.sampleID,r.majorCluster):i for r in units.itertuples(index=False)};targets=np.array([index[(s,base['target'])] for s in a['samples']]);materials=a['corrected_material_for_next_inference'];blocks=a['corrected_blocks_for_next_inference'];reviewed=a['corrected_identity_reviewed_for_next_inference'];assert len(targets)==len(blocks)
  for ix,material in enumerate(sorted(set(materials))):
   mask=materials==material;cid=base['group_id']+'_M'+str(ix+1).zfill(2);samples=a['samples'][mask];bb=blocks[mask];rv=reviewed[mask];Nlib=int(mask.sum());G=int(len(set(bb)));passed=Nlib>=20 and G>=10;ti=targets[mask]
   keys={k:base[k] for k in ['group_id','gse_id','stratum','cancerType','sample_class','source','target']};keys.update(context_id=cid,material=str(material),n_libraries=Nlib,n_analysis_blocks=G,n_reviewed_blocks=int(len(set(bb[rv]))),coverage_ge20libs_ge10blocks=passed)
   primary=ss[0][ti];sens=ss[1][ti];det=dd[0][ti];support=(det>=10).sum(0);variable=np.array([np.ptp(x[np.isfinite(x)])>1e-10 if np.isfinite(x).sum()>1 else False for x in primary.T]);tfready=(cc[0]>=20)&(support>=.8*Nlib)&variable&passed
   l=a['ligand_meanlog'][mask];r=a['receptor_meanlog'][mask];lp=a['ligand_min_positive_fraction'][mask];rp=a['receptor_min_positive_fraction'][mask];measured=a['native_measured'];ligand_var=np.ptp(np.nan_to_num(l),axis=0)>1e-12
   # Prospective RNA support: >=1% cells per end, >=50% ligand/80% receptor libraries. Not binding.
   lrsupport=measured&ligand_var&((lp>=.01).sum(0)>=.5*Nlib)&((rp>=.01).sum(0)>=.8*Nlib)&passed
   out=keys|dict(TFs_primary_ready=int(tfready.sum()),LR_RNA_support_ready=int(lrsupport.sum()),both_program_and_LR_ready=bool(tfready.any() and lrsupport.any()))
   contexts.append(out)
   np.savez_compressed(O/'private_context_arrays'/(cid+'.npz'),samples=samples,blocks=bb,identity_reviewed=rv,primary_TF_scores=primary,sensitivity_TF_scores=sens,TFs=TFs,TF_ready=tfready,LR_RNA_support_ready=lrsupport,original_group_id=base['group_id'],material=material)
   for j,tf in enumerate(TFs):tfrows.append(keys|dict(TF=tf,measured_targets=int(cc[0][j]),libraries_detected_at_least10targets=int(support[j]),variable=bool(variable[j]),TF_ready=bool(tfready[j])))
   lrrows.append(keys|dict(native_resource_pairs=int(measured.sum()),LR_RNA_support_ready=int(lrsupport.sum())))
 pd.DataFrame(coverage_rows).to_csv(O/'all_source_TF_coverage.tsv',sep='\t',index=False);pd.DataFrame(source_summary).to_csv(O/'source_TF_summary.tsv',sep='\t',index=False);pd.DataFrame(contexts).to_csv(O/'all_corrected_context_readiness.tsv',sep='\t',index=False);pd.DataFrame(tfrows).to_csv(O/'private_all_context_TF_readiness.tsv.gz',sep='\t',index=False);pd.DataFrame(lrrows).to_csv(O/'all_corrected_context_LR_support.tsv',sep='\t',index=False)
 out=dict(status='COMPLETE_TF_PREFLIGHT_PENDING_INDEPENDENT_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),corrected_contexts=len(contexts),context_TF_rows=len(tfrows),qualified_contexts=int(sum(x['coverage_ge20libs_ge10blocks'] for x in contexts)),contexts_with_both_TF_and_LR_support=int(sum(x['both_program_and_LR_ready'] for x in contexts)),source_summary=source_summary,TF_network_sha256=sha(D/'human_collectri_signed_complexes.tsv'),no_new_association_p=True,not_receptor_specific_response_or_protein_binding=True,outputs={str(f.relative_to(O)):sha(f) for f in O.rglob('*') if f.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='outputs'}),flush=True)
if __name__=='__main__':main()
