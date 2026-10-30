"""Source-balanced descriptive marker effects for accepted joint clusters.

This supplements official all-cell Wilcoxon; it is neither patient-level DE
nor an independent identity truth set. No final label is written.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def main(R):
 started=time.monotonic();cfg=load(R/'config/spec_v1.json');M=Path(cfg['marker_producer']);A=Path(cfg['marker_auditor']);H=Path(load(M/'config/spec_v1.json')['producer']);hc=load(H/'config/spec_v1.json');sys.path.insert(0,str(Path(hc['joined'])/'scripts'))
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 import psutil
 storage_gate(R,cfg);assert psutil.virtual_memory().available>=cfg['minimum_available_RAM_bytes']
 accepted=load(A/'outputs/all_native_profiles_acceptance_v1.json');contract=load(M/'outputs/marker_input_contract_v1.json');idx=load(M/'outputs/native_profiles_index_v1.json');catalog=load(M/'outputs/cluster_membership_catalog_v1.json');assert [x['source'] for x in idx]==[x['source'] for x in accepted['source_checks']]
 genes=contract['common_gene_names'];ng=contract['n_clusters'];N=contract['n_fit'];assert len(catalog)==ng
 source_n=[];effect=[];detect_effect=[];own=[];detect=[]
 for r,a in zip(idx,accepted['source_checks']):
  p=M/'outputs'/r['path'];assert sha(p)==r['SHA256']==a['source_profile_SHA256']
  with np.load(p,allow_pickle=False) as f:
   names=f['native_gene_names'].tolist();by={g:i for i,g in enumerate(names)};cols=[by[g] for g in genes];n=f['n_cells'].astype(np.int64);log=f['logCP10K_sum'][:,cols];det=f['detected'][:,cols]
  assert n.sum()==r['n_cells'];total=int(n.sum());other=total-n;valid=(n>0)&(other>0)
  own_mean=np.divide(log,n[:,None],out=np.full_like(log,np.nan),where=n[:,None]>0);own_det=np.divide(det,n[:,None],out=np.full_like(det,np.nan),where=n[:,None]>0)
  other_mean=np.divide(log.sum(axis=0)-log,other[:,None],out=np.full_like(log,np.nan),where=other[:,None]>0);other_det=np.divide(det.sum(axis=0)-det,other[:,None],out=np.full_like(det,np.nan),where=other[:,None]>0)
  e=own_mean-other_mean;d=own_det-other_det;e[~valid]=np.nan;d[~valid]=np.nan;effect.append(e);detect_effect.append(d);own.append(own_mean);detect.append(own_det);source_n.append(n)
 effect=np.stack(effect);detect_effect=np.stack(detect_effect);own=np.stack(own);detect=np.stack(detect);n=np.stack(source_n);assert n.sum()==N
 rows=[];tops=[];scope=[];cut=int(cfg['minimum_source_group_and_comparator_cells']);assert cut==30
 for k in range(ng):
  valid=(n[:,k]>0)&(n.sum(axis=1)-n[:,k]>0);adequate=(n[:,k]>=cut)&(n.sum(axis=1)-n[:,k]>=cut);scope.append(dict(cluster=str(k),n_cells=int(n[:,k].sum()),n_observed_sources=int((n[:,k]>0).sum()),n_comparable_sources=int(valid.sum()),n_sources_with_minimum_cells=int(adequate.sum()),minimum_group_and_comparator_cells=cut,membership_SHA256=catalog[k]['exact_membership_SHA256']))
  # Empty source/group is missing, never a zero effect or negative marker.
  for j,gene in enumerate(genes):
   row=dict(cluster=str(k),gene=gene,n_comparable_sources=int(valid.sum()),n_sources_with_minimum_cells=int(adequate.sum()),n_positive_sources=int((effect[valid,k,j]>0).sum()),n_negative_sources=int((effect[valid,k,j]<0).sum()),mean_within_source_logCP10K_difference=float(np.mean(effect[valid,k,j])) if valid.any() else None,median_within_source_logCP10K_difference=float(np.median(effect[valid,k,j])) if valid.any() else None,mean_within_source_detection_difference=float(np.mean(detect_effect[valid,k,j])) if valid.any() else None,n_positive_adequate_sources=int((effect[adequate,k,j]>0).sum()),mean_adequate_source_logCP10K_difference=float(np.mean(effect[adequate,k,j])) if adequate.any() else None,mean_adequate_source_detection_difference=float(np.mean(detect_effect[adequate,k,j])) if adequate.any() else None,source_equal_mean_own_detection=float(np.mean(detect[valid,k,j])) if valid.any() else None,cell_weighted_own_mean_logCP10K=float(np.nansum(own[:,k,j]*n[:,k])/n[:,k].sum()))
   rows.append(row)
  candidates=pd.DataFrame(rows[-len(genes):]);order='mean_adequate_source_logCP10K_difference' if adequate.any() else 'mean_within_source_logCP10K_difference';candidates=candidates.sort_values([order,'gene'],ascending=[False,True],kind='stable');tops.extend(candidates.head(cfg['top_per_cluster']).to_dict('records'))
 pd.DataFrame(rows).to_csv(R/'outputs/all_common_gene_within_source_effects_v1.tsv.gz',sep='\t',index=False);pd.DataFrame(tops).to_csv(R/'outputs/source_balanced_top_candidates_v1.tsv',sep='\t',index=False)
 np.savez_compressed(R/'outputs/private_source_cluster_common_gene_effects_v1.npz',sources=np.asarray([r['source'] for r in idx]),genes=np.asarray(genes),n_cells=n,within_source_logCP10K_difference=effect,within_source_detection_difference=detect_effect)
 save(R/'outputs/source_group_coverage_v1.json',scope);storage_gate(R,cfg)
 save(R/'outputs/source_effect_contract_v1.json',dict(n_fit=N,n_sources=len(idx),n_clusters=ng,n_common_genes=len(genes),native_profiles_acceptance_SHA256=sha(A/'outputs/all_native_profiles_acceptance_v1.json'),all_native_profile_input_SHA256=[dict(source=r['source'],SHA256=r['SHA256']) for r in idx],source_equal_effect_not_patient_DEG=True,cell_level_Wilcoxon_not_replaced=True,source_and_cancer_confounding_not_removed=True,minimum_cell_cutoff_inherited_from_existing_cluster_review_policy=cut,no_P_values_or_final_identity_assigned=True))
 save(R/'provenance/terminal_v1.json',dict(status='WITHIN_SOURCE_JOINT_CLUSTER_EFFECTS_AND_COVERAGE_COMPLETE_REQUIRES_COHERENT_MARKERS_AND_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,n_clusters=ng,n_common_genes=len(genes),wall_seconds=time.monotonic()-started,contract_SHA256=sha(R/'outputs/source_effect_contract_v1.json')))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
