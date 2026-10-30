"""Independent full-source marker-profile replay and bounded rank cases."""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
def load(p):return json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def check(R,stage,**v):save(R/'checkpoints'/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**v))
def await_file(p,P):
 while not p.exists():
  t=P/'provenance/supervisor_terminal_v1.json'
  if t.exists():assert load(t)['actual_child_exit_code']==0,'Producer failed; preserve failed version'
  time.sleep(5)
def read_rows(csr,rows,integer):
 unique,inverse=np.unique(rows,return_inverse=True);parts=[];start=0
 while start<len(unique):
  end=int(np.searchsorted(unique,unique[start]+1024));lo=int(unique[start]);hi=int(unique[end-1])+1;ptr=np.asarray(csr['indptr'][lo:hi+1],np.int64);a=int(ptr[0]);b=int(ptr[-1]);data=np.asarray(csr['data'][a:b],np.float64);ix=np.asarray(csr['indices'][a:b],np.int64)
  assert np.isfinite(data).all() and (data>=0).all()
  if integer:assert np.equal(data,np.floor(data)).all()
  x=sp.csr_matrix((data,ix,ptr-a),shape=(hi-lo,int(csr.attrs['shape'][1])));x.sum_duplicates();x.eliminate_zeros();parts.append(x[unique[start:end]-lo]);start=end
 return sp.vstack(parts,format='csr')[inverse]
def main(R):
 cfg=load(R/'config/spec_v1.json');P=Path(cfg['marker_producer']);await_file(P/'outputs/marker_input_contract_v1.json',P);mc=load(P/'config/spec_v1.json');H=Path(mc['producer']);G=Path(mc['graph']);hc=load(H/'config/spec_v1.json');J=Path(hc['joined']);sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from anndata._io.specs import read_elem
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);contract=load(P/'outputs/marker_input_contract_v1.json');genes=contract['common_gene_names'];ng=contract['n_clusters'];N=contract['n_fit'];lineage=load(H/'outputs/source_input_lineage_v1.json');labels=pd.read_csv(G/'outputs/private_full_graph_cluster_axes_v1.tsv.gz',sep='\t',usecols=['leiden_0_5']).iloc[:,0].to_numpy(np.int32)
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');aggregate_log=np.zeros((ng,len(genes)));aggregate_detect=np.zeros((ng,len(genes)));aggregate_n=np.zeros(ng,np.int64);nonestimable=None;source_checks=[]
 wa=load(hc['working_acceptance']);wm=load(wa['producer_manifest']['path']);W=Path(wm['index']['path']).parent;workidx={r['source']:r for r in load(wm['index']['path'])};crosswalk=[]
 for code,r in enumerate(lineage):
  source=r['source'];record=P/'outputs'/('native_cluster_profiles_'+source+'.json');await_file(record,P);m=load(record);p=P/'outputs'/m['path'];assert sha(p)==m['SHA256'];wh=np.flatnonzero(axes['source_code']==code);z=native.store(source);var=read_elem(z['var']);lookup=native.feature_lookup(source);names=[];cols=[]
  for j,idx in enumerate(var.index):
   name=str(var.iloc[j]['gene_symbol']) if 'gene_symbol' in var else str(idx)
   if lookup.get(name,set())=={j}:names.append(name);cols.append(j)
  with np.load(p,allow_pickle=False) as f:saved={k:f[k] for k in f.files}
  assert saved['native_gene_names'].tolist()==names;np.testing.assert_array_equal(saved['native_columns'],cols);n=np.bincount(labels[wh],minlength=ng);np.testing.assert_array_equal(saved['n_cells'],n);expected={k:np.zeros((ng,len(names)),np.float64) for k in ['CP10K_sum','logCP10K_sum','detected']};delta=z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None
  for lo in range(0,len(wh),1024):
   take=wh[lo:lo+1024];raw=read_rows(z[r['counts_key']],axes['native_row'][take],True);np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),axes['native_total_RNA'][take]);x=raw
   if delta is not None:x=(x-read_rows(delta,axes['native_row'][take],False)).tocsr()
   x.eliminate_zeros();assert x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();RNA=np.asarray(x.sum(axis=1)).ravel();np.testing.assert_allclose(RNA,totals[take],rtol=1e-10,atol=1e-7);coo=x[:,cols].tocoo();which=labels[take][coo.row];v=coo.data*10000/RNA[coo.row];np.add.at(expected['CP10K_sum'],(which,coo.col),v);np.add.at(expected['logCP10K_sum'],(which,coo.col),np.log1p(v));np.add.at(expected['detected'],(which,coo.col),1)
  for k in expected:
   if k=='detected':np.testing.assert_array_equal(expected[k],saved[k])
   else:np.testing.assert_allclose(expected[k],saved[k],rtol=2e-10,atol=1e-7)
  by={g:i for i,g in enumerate(names)};positions=[by[g] for g in genes];aggregate_log+=expected['logCP10K_sum'][:,positions];aggregate_detect+=expected['detected'][:,positions];aggregate_n+=n
  if nonestimable is None and (n==0).any():
   k=int(np.flatnonzero(n==0)[0]);assert all(not expected[a][k].any() for a in expected);nonestimable=dict(source=source,cluster=k,status='SOURCE_DID_NOT_CAPTURE_THIS_GROUP_NOT_A_NEGATIVE_MARKER_RESULT')
  wr=workidx[source];meta=load(W/wr['metadata'])
  with np.load(W/wr['arrays'],allow_pickle=False) as f:
   fit=np.flatnonzero(f['fit_pool_candidate_not_final_QC']);np.testing.assert_array_equal(f['native_row'][fit],axes['native_row'][wh]);major=np.asarray(meta['candidate_dictionaries']['consensus_major'])[f['consensus_major_code'][fit]];fine=np.asarray(meta['candidate_dictionaries']['consensus_celltype'])[f['consensus_celltype_code'][fit]];oldgroup=f['annotation_group_code'][fit]
  tab=pd.DataFrame(dict(cluster=labels[wh],previous_major=major,previous_celltype=fine,previous_source_group=oldgroup));tab=tab.groupby(list(tab.columns),observed=True).size().reset_index(name='n_cells');tab.insert(0,'source',source);crosswalk.extend(tab.to_dict('records'));source_checks.append(dict(source=source,n_cells=len(wh),all_native_genes_and_all_groups_replayed=True,source_profile_SHA256=m['SHA256']));save(R/'outputs'/('source_profile_replay_'+source+'.json'),source_checks[-1]);check(R,'independent_all_native_source_profile_accepted',source=source,n_cells=len(wh));native._fields.clear()
 assert aggregate_n.sum()==N and nonestimable is not None;save(R/'outputs/all_native_profiles_acceptance_v1.json',dict(n_sources=len(source_checks),n_fit=N,all_native_source_gene_axes_and_all_group_sufficient_statistics_replayed=True,source_checks=source_checks,nonestimable_source_group_case=nonestimable,not_patient_independent_inference=True))
 # Derive positive/negative arithmetic probes from full means, not clinical effects.
 effect=aggregate_log[0]/aggregate_n[0]-(aggregate_log.sum(axis=0)-aggregate_log[0])/(N-aggregate_n[0]);chosen=[int(np.argmax(effect)),int(np.argmin(effect))];assert effect[chosen[0]]>0 and effect[chosen[1]]<0;vectors=np.zeros((N,2),np.float32)
 for code,r in enumerate(lineage):
  source=r['source'];wh=np.flatnonzero(axes['source_code']==code);z=native.store(source);var=read_elem(z['var']);lookup=native.feature_lookup(source);columns=[]
  for g in [genes[j] for j in chosen]:assert len(lookup[g])==1;columns.append(next(iter(lookup[g])))
  delta=z['layers/ambient_removed_counts'] if 'ambient_removed_counts' in z['layers'] else None
  for lo in range(0,len(wh),1024):
   take=wh[lo:lo+1024];x=read_rows(z[r['counts_key']],axes['native_row'][take],True)
   if delta is not None:x=(x-read_rows(delta,axes['native_row'][take],False)).tocsr()
   v=x[:,columns].toarray();v[v<0]=0;vectors[take]=np.log1p(v*(10000/totals[take])[:,None]).astype(np.float32)
  native._fields.clear()
 await_file(P/'provenance/supervisor_terminal_v1.json',P);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;df=pd.read_csv(P/'outputs/all_cell_common_gene_Wilcoxon_v1.tsv.gz',sep='\t',dtype={'group':str,'names':str},keep_default_na=False);assert len(df)==ng*len(genes)
 from statsmodels.stats.multitest import multipletests
 for name,g in df.groupby('group',observed=True):np.testing.assert_allclose(multipletests(g.pvals.to_numpy(),method='fdr_bh')[1],g.pvals_adj,rtol=1e-10,atol=1e-12)
 from scipy.stats import mannwhitneyu,norm
 cases=[];inside=labels==0;n=int(inside.sum());rest=N-n
 for col,j in enumerate(chosen):
  v=vectors[:,col];u=mannwhitneyu(v[inside],v[~inside],method='asymptotic',use_continuity=False);_,ties=np.unique(v,return_counts=True);ties=ties.astype(np.float64);variance=n*float(rest)*(N+1-np.sum(ties**3-ties)/(N*(N-1)))/12;z=(u.statistic-n*float(rest)/2)/np.sqrt(variance);expected_p=2*norm.sf(abs(z));row=df[(df.group=='0')&(df.names==genes[j])];assert len(row)==1;np.testing.assert_allclose(float(row.scores.iloc[0]),z,rtol=1e-5,atol=1e-4);np.testing.assert_allclose(float(row.pvals.iloc[0]),expected_p,rtol=1e-7,atol=1e-12);cases.append(dict(direction='positive' if col==0 else 'negative',cluster='0',gene=genes[j],independent_Mann_Whitney_and_tie_corrected_score_pass=True))
 expected=pd.DataFrame(crosswalk);actual=pd.read_csv(P/'outputs/previous_working_labels_to_new_clusters_v1.tsv',sep='\t',dtype={'source':str,'previous_major':str,'previous_celltype':str},keep_default_na=False);keys=['source','cluster','previous_major','previous_celltype','previous_source_group'];pd.testing.assert_frame_equal(actual.sort_values(keys).reset_index(drop=True),expected[actual.columns].sort_values(keys).reset_index(drop=True),check_dtype=False)
 save(R/'outputs/full_cluster_marker_acceptance_v1.json',dict(all_native_profiles_accepted=True,all_common_gene_BH_values_independently_recomputed=True,positive_negative_rank_cases=cases,nonestimable_source_group_case=nonestimable,all_previous_label_crosswalk_counts_reproduced=True,n_fit=N,n_sources=len(source_checks),cell_level_marker_screen_not_patient_DEG=True,no_final_identity_promoted=True))
 save(R/'provenance/terminal_v1.json',dict(status='FULL_NATIVE_CLUSTER_PROFILES_ALL_BH_AND_INDEPENDENT_POSITIVE_NEGATIVE_NONESTIMABLE_CASES_ACCEPTED_REQUIRES_BIOLOGICAL_MARKER_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_fit=N,n_sources=len(source_checks)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
