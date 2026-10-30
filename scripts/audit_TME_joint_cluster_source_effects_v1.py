"""Replay all descriptive source effects and candidate ordering from profiles."""
import argparse,datetime,hashlib,json,time
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
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['producer']);terminal=P/'provenance/supervisor_terminal_v1.json'
 while not terminal.exists():time.sleep(5)
 assert load(terminal)['actual_child_exit_code']==0;pc=load(P/'config/spec_v1.json');M=Path(pc['marker_producer']);c=load(P/'outputs/source_effect_contract_v1.json');assert sha(P/'outputs/source_effect_contract_v1.json')==load(P/'provenance/terminal_v1.json')['contract_SHA256'];genes=load(M/'outputs/marker_input_contract_v1.json')['common_gene_names'];index=load(M/'outputs/native_profiles_index_v1.json');ng=c['n_clusters'];cut=pc['minimum_source_group_and_comparator_cells']
 with np.load(P/'outputs/private_source_cluster_common_gene_effects_v1.npz',allow_pickle=False) as f:saved={k:f[k] for k in f.files}
 assert saved['sources'].tolist()==[r['source'] for r in index] and saved['genes'].tolist()==genes
 total_log=np.zeros((ng,len(genes)));alln=[];allown=[]
 for i,r in enumerate(index):
  p=M/'outputs'/r['path'];assert sha(p)==r['SHA256']
  with np.load(p,allow_pickle=False) as f:
   names=f['native_gene_names'].tolist();by={g:j for j,g in enumerate(names)};columns=[by[g] for g in genes];n=f['n_cells'];lg=f['logCP10K_sum'][:,columns];det=f['detected'][:,columns]
  np.testing.assert_array_equal(n,saved['n_cells'][i]);expected=np.full((ng,len(genes)),np.nan);edet=expected.copy();denom=int(n.sum());alln.append(n);own=np.full_like(expected,np.nan);total_log+=lg
  for k in np.flatnonzero(n):
   own[k]=det[k]/n[k]
   if denom==n[k]:continue
   rest=np.delete(lg,k,axis=0).sum(axis=0);rest_det=np.delete(det,k,axis=0).sum(axis=0);expected[k]=lg[k]/n[k]-rest/(denom-n[k]);edet[k]=det[k]/n[k]-rest_det/(denom-n[k])
  np.testing.assert_allclose(expected,saved['within_source_logCP10K_difference'][i],rtol=1e-9,atol=1e-10,equal_nan=True);np.testing.assert_allclose(edet,saved['within_source_detection_difference'][i],rtol=1e-9,atol=1e-10,equal_nan=True);allown.append(own)
 n=np.stack(alln);own=np.stack(allown);df=pd.read_csv(P/'outputs/all_common_gene_within_source_effects_v1.tsv.gz',sep='\t',dtype={'cluster':str,'gene':str});assert len(df)==ng*len(genes) and not df[['cluster','gene']].duplicated().any();top=pd.read_csv(P/'outputs/source_balanced_top_candidates_v1.tsv',sep='\t',dtype={'cluster':str,'gene':str});assert len(top)==ng*pc['top_per_cluster'];scope=load(P/'outputs/source_group_coverage_v1.json');assert len(scope)==ng;positive_case=negative_case=unobserved_case=False
 for k in range(ng):
  rows=df[df.cluster==str(k)].set_index('gene').loc[genes];valid=(n[:,k]>0)&(n.sum(axis=1)>n[:,k]);adequate=(n[:,k]>=cut)&(n.sum(axis=1)-n[:,k]>=cut);e=saved['within_source_logCP10K_difference'][valid,k];d=saved['within_source_detection_difference'][valid,k];ea=saved['within_source_logCP10K_difference'][adequate,k];da=saved['within_source_detection_difference'][adequate,k]
  assert (rows.n_comparable_sources==valid.sum()).all() and (rows.n_sources_with_minimum_cells==adequate.sum()).all()
  expected={'n_positive_sources':(e>0).sum(axis=0),'n_negative_sources':(e<0).sum(axis=0),'mean_within_source_logCP10K_difference':e.mean(axis=0),'median_within_source_logCP10K_difference':np.median(e,axis=0),'mean_within_source_detection_difference':d.mean(axis=0),'n_positive_adequate_sources':(ea>0).sum(axis=0),'mean_adequate_source_logCP10K_difference':ea.mean(axis=0) if adequate.any() else np.full(len(genes),np.nan),'mean_adequate_source_detection_difference':da.mean(axis=0) if adequate.any() else np.full(len(genes),np.nan),'source_equal_mean_own_detection':own[valid,k].mean(axis=0),'cell_weighted_own_mean_logCP10K':total_log[k]/n[:,k].sum()}
  for name,v in expected.items():np.testing.assert_allclose(rows[name],v,rtol=1e-9,atol=1e-10,equal_nan=True)
  order='mean_adequate_source_logCP10K_difference' if adequate.any() else 'mean_within_source_logCP10K_difference';expected_top=rows.reset_index().sort_values([order,'gene'],ascending=[False,True],kind='stable').head(pc['top_per_cluster']);actual_top=top[top.cluster==str(k)];assert actual_top.gene.tolist()==expected_top.gene.tolist();assert scope[k]['n_cells']==int(n[:,k].sum()) and scope[k]['n_comparable_sources']==int(valid.sum()) and scope[k]['n_sources_with_minimum_cells']==int(adequate.sum())
  positive_case|=bool((e>0).any());negative_case|=bool((e<0).any())
  if (n[:,k]==0).any():assert np.isnan(saved['within_source_logCP10K_difference'][n[:,k]==0,k]).all();unobserved_case=True
 assert positive_case and negative_case and unobserved_case
 save(R/'outputs/accepted_within_source_joint_cluster_effects_v1.json',dict(n_fit=c['n_fit'],all_source_gene_and_cluster_effects_replayed=True,all_source_balanced_summary_columns_replayed=True,all_candidate_orders_reproduced=True,positive_and_negative_effect_cases=True,unobserved_source_cluster_missing_not_negative=True,producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),source_effects_not_independent_patient_DEG_or_final_labels=True))
 save(R/'provenance/terminal_v1.json',dict(status='ALL_WITHIN_SOURCE_JOINT_CLUSTER_EFFECTS_AND_RANKING_INDEPENDENTLY_REPLAYED_REQUIRES_BIOLOGICAL_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=c['n_fit']))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
