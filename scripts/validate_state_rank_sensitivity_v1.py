"""Official pyUCell rank and measured-universe sensitivity on frozen stratified rows."""
import sys,json,hashlib,datetime,traceback,time
from pathlib import Path
import numpy as np
import pandas as pd

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def main(R,P):
 start=time.monotonic();cfg=load(R/'config.json');assert sha(P/'provenance/terminal_v1.json')==cfg['producer_terminal_SHA256'];pc=load(P/'config/spec_v1.json');H=Path(pc['producer']);T=Path(pc['existing_trial']);hc=load(H/'config/spec_v1.json');J=Path(hc['joined'])
 sys.path[:0]=[str(P/'vendor'),str(J/'scripts'),str(T/'scripts')]
 import pyucell,anndata
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
 assert pyucell.__version__=='0.8.0'
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);genes=load(H/'outputs/scope_and_common_genes_v1.json')['genes'];lineage={x['source']:x for x in load(H/'outputs/source_input_lineage_v1.json')}
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 totals=np.load(H/'outputs/private_analysis_RNA_totals_v1.npy',mmap_mode='r');axis=load(P/'outputs/score_program_axis_v1.json');meta=axis['reference_metadata'];names=axis['program_ids'];sig={r['program_id']:r['common_universe_genes'] for r in meta};summary=[];n=0
 for sourceinfo in load(P/'outputs/source_score_index_v1.json'):
  s=sourceinfo['source'];o=pd.read_csv(P/'outputs'/('private_'+s+'_axis_v1.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str});assert o.role.nunique()==1 and o.role.iloc[0] in ['discovery','description_only']
  # Round-robin material/subtype strata; stable hash ordering, independent of scores.
  strata=[]
  for key,ix in o.groupby(['sample','subtype'],sort=True).indices.items():
   ordered=sorted(ix,key=lambda i:hashlib.sha256((s+'|'+str(int(o.row.iloc[i]))+'|20261030').encode()).digest());strata.append(ordered[:cfg['max_per_unit_subtype']])
  chosen=[a[k] for k in range(cfg['max_per_unit_subtype']) for a in strata if len(a)>k][:cfg['max_per_source']];chosen=np.asarray(chosen,int);g=o.iloc[chosen].copy();rows=g.row.to_numpy();r=lineage[s];ng,columns=native_measured_genes(native,s);pos=dict(zip(ng,columns));assert len(pos)==len(ng);z=native.store(s);x=gather_csr_rows(z[r['counts_key']],axes['native_row'][rows],max_span=8192).astype(np.float64).tocsr()
  if 'ambient_removed_counts' in z['layers']:x=(x-gather_nonnegative_float_delta(z['layers/ambient_removed_counts'],axes['native_row'][rows])).tocsr()
  x.eliminate_zeros();assert np.isfinite(x.data).all() and x.data.min(initial=0)>=-1e-8;x.data[x.data<0]=0;x.eliminate_zeros();np.testing.assert_allclose(np.asarray(x.sum(1)).ravel(),totals[rows],rtol=1e-10,atol=1e-7)
  common=x[:,[pos[t] for t in genes]].astype(np.float32);full=x[:,columns].astype(np.float32);scores={}
  for universe,matrix,gn in [('common',common,genes),('native',full,ng)]:
   for rank in cfg['ranks']:
    a=anndata.AnnData(matrix,var=pd.DataFrame(index=gn));pyucell.compute_ucell_scores(a,sig,max_rank=rank,missing_genes='skip',ties_method='average',chunk_size=256,n_jobs=2,device='cpu');values=a.obs[[k+'_UCell' for k in names]].to_numpy(np.float32);assert np.isfinite(values).all() and values.min()>=0 and values.max()<=1.000001;scores[universe+'_'+str(rank)]=values
  baseline=np.load(P/'outputs'/('private_'+s+'_UCell_v1.npy'),mmap_mode='r')[chosen];np.testing.assert_allclose(scores['common_1500'],baseline,atol=1e-7,rtol=1e-6)
  g.to_csv(R/'outputs'/('private_'+s+'_probe_axis.tsv.gz'),sep='\t',index=False);np.savez_compressed(R/'outputs'/('private_'+s+'_scores.npz'),**scores)
  for label,values in scores.items():
   diff=values-baseline
   for j,name in enumerate(names):summary.append(dict(source=s,role=o.role.iloc[0],program_id=name,comparison=label,n_cells=len(g),mean_absolute_difference=float(np.mean(np.abs(diff[:,j]))),mean_signed_difference=float(diff[:,j].mean()),n_native_measured_genes=len(ng)))
  n+=len(g);save(R/'checkpoints'/(s+'.json'),dict(completed_UTC=now(),source=s,n_sampled_cells=n));native._fields.clear()
 pd.DataFrame(summary).to_csv(R/'outputs/aggregate_source_program_rank_sensitivity.tsv',sep='\t',index=False)
 save(R/'outputs/summary.json',dict(completed_UTC=now(),n_cells=n,n_sources=37,n_programs=126,n_settings=6,official_pyucell_version=pyucell.__version__,signature_genes_fixed_to_scored_common_subset=True,native_comparison_changes_rank_universe_only=True,not_full_native_signature_rescue=True,not_patient_or_biological_validation=True,no_module_fit=True,wall_seconds=time.monotonic()-start));save(R/'provenance/terminal.json',dict(status='RANK_SENSITIVITY_COMPLETE',actual_exit_code=0,**load(R/'outputs/summary.json')))
if __name__=='__main__':
 R=Path(sys.argv[1]);P=Path(sys.argv[2])
 try:main(R,P)
 except BaseException:save(R/'provenance/failure.json',dict(failed_UTC=now(),traceback=traceback.format_exc()));raise
