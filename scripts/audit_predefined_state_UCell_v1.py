"""Independent axis, all-value and group-mean checks plus native scalar-score probes."""
from pathlib import Path
import sys,json,hashlib,time,datetime,traceback

def save(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(A,R):
 started=time.monotonic()
 while not (R/'provenance/supervisor_terminal_v1.json').exists():
  if time.monotonic()-started>11000:raise TimeoutError('producer did not finish within audit wait budget')
  time.sleep(5)
 assert load(R/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 terminal=load(R/'provenance/terminal_v1.json');cfg=load(R/'config/spec_v1.json')
 import numpy as np,pandas as pd
 from scipy.stats import rankdata
 H=Path(cfg['producer']);tc=load(Path(cfg['existing_trial'])/'config/spec_v1.json');J=Path(load(H/'config/spec_v1.json')['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(Path(cfg['existing_trial'])/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_annotation_review_v1 import native_measured_genes
 from native_panel_aggregate_v1 import gather_csr_rows
 from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);lineage={r['source']:r for r in load(H/'outputs/source_input_lineage_v1.json')};scope=load(H/'outputs/scope_and_common_genes_v1.json');genes=scope['genes'];gpos={g:i for i,g in enumerate(genes)};programs=load(R/'outputs/score_program_axis_v1.json');names=programs['program_ids'];meta=programs['reference_metadata']
 with np.load(H/'outputs/private_cell_axes_v1.npz',allow_pickle=False) as f:axes={k:f[k] for k in f.files}
 oldobs=pd.concat([pd.read_csv(Path(cfg['existing_trial'])/'outputs'/name,sep='\t',dtype={'sample':str,'subtype':str}) for name in ['private_obs_discovery_cells_v1.tsv.gz','private_obs_other_development_cells_v1.tsv.gz']],ignore_index=True)
 obs=pd.read_csv(R/'outputs/private_observations_v1.tsv.gz',sep='\t',dtype={'sample':str,'subtype':str});assert len(obs)==terminal['n_cells']==len(oldobs)==776270
 pd.testing.assert_frame_equal(obs[['row','source','sample','subtype']],oldobs[['row','source','sample','subtype']],check_dtype=False);assert not obs.duplicated(['source','row']).any();assert np.array_equal(obs.row.to_numpy(),np.load(R/'outputs/private_score_row_order_v1.npy',allow_pickle=False))
 split=pd.read_csv(Path(tc['split_run'])/'outputs/source_split.tsv',sep='\t');assert not set(obs.source)&set(split.source[split.split=='lockbox'])
 cat=load(R/'resources/predefined_state_catalog_v1.json');assert sha(R/'resources/predefined_state_catalog_v1.json')==cfg['catalog_SHA256'];assert names==[r['program_id'] for r in cat['programs'] if r['scoreable_TNK_pilot']]
 total=0;manual=0;maxerr=0.;checked_means=0;rows=[]
 for source,group in obs.groupby('source',sort=True):
  axis=pd.read_csv(R/'outputs'/('private_'+source+'_axis_v1.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str});pd.testing.assert_frame_equal(axis[['row','source','sample','subtype']].reset_index(drop=True),group[['row','source','sample','subtype']].reset_index(drop=True),check_dtype=False)
  sf=R/'outputs'/('private_'+source+'_UCell_v1.npy');y=np.load(sf,mmap_mode='r',allow_pickle=False);assert y.shape==(len(group),len(names));assert np.isfinite(y).all() and y.min()>=-1e-6 and y.max()<=1.000001
  summ=load(R/'outputs'/(source+'_summary_v1.json'));assert summ['scores_SHA256']==sha(sf);total+=y.size
  # An independent sorted-segment sum reproduces every persisted group mean.
  stored=pd.read_csv(R/'outputs'/('private_'+source+'_unit_subtype_means_v1.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str});mapping={(str(r['sample']),str(r['subtype'])):i for i,r in stored.iterrows()};codes=np.array([mapping[(str(s),str(t))] for s,t in zip(axis['sample'],axis.subtype)]);counts=np.bincount(codes,minlength=len(stored));assert np.array_equal(counts,stored.n_cells.to_numpy());order=np.argsort(codes,kind='stable');starts=np.r_[0,np.cumsum(counts)[:-1]];means=np.add.reduceat(np.asarray(y)[order].astype(np.float64),starts,axis=0)/counts[:,None];np.testing.assert_allclose(means,stored[names].to_numpy(),rtol=2e-5,atol=2e-6);checked_means+=means.size
  probe=np.unique(np.linspace(0,len(axis)-1,min(3,len(axis)),dtype=int));r=lineage[source];z=native.store(source);nat_names,col=native_measured_genes(native,source);physical=dict(zip(nat_names,col));ph=np.asarray([physical[g] for g in genes]);nr=axes['native_row'][axis.row.to_numpy()[probe]];x=gather_csr_rows(z[r['counts_key']],nr,max_span=8192).astype(np.float64).tocsr()
  if 'ambient_removed_counts' in z['layers']:x=(x-gather_nonnegative_float_delta(z['layers/ambient_removed_counts'],nr)).tocsr()
  x.eliminate_zeros();x.data[x.data<0]=0;x=x[:,ph].toarray().astype(np.float32)
  for k,local in enumerate(probe):
   vec=x[k];positive=np.flatnonzero(vec>0);ranks=np.full(len(vec),cfg['max_rank'],dtype=np.float64);pr=rankdata(-vec[positive],method='average');ranks[positive]=np.minimum(pr,cfg['max_rank'])
   for j,m in enumerate(meta):
    n=len(m['common_universe_genes']);rank_sum=sum(ranks[gpos[g]] for g in m['common_universe_genes']);base=n*(n+1)/2.;v=1-(rank_sum-base)/(n*cfg['max_rank']-base);err=abs(v-float(y[local,j]));maxerr=max(maxerr,err);assert err<3e-6;manual+=1
  rows.append(dict(source=source,n_cells=len(axis),n_values=int(y.size),n_means_checked=int(means.size),manual_score_values=len(probe)*len(meta)));native._fields.clear()
  save(A/'checkpoints'/(source+'_accepted_v1.json'),rows[-1])
 assert sum(r['n_cells'] for r in rows)==776270
 save(A/'outputs/acceptance_v1.json',dict(status='OFFICIAL_PREDEFINED_UCELL_AXES_ALL_VALUES_GROUP_MEANS_AND_NATIVE_SCORE_PROBES_ACCEPTED_NOT_BIOLOGICAL_STATE_VALIDATION',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),producer_terminal_SHA256=sha(R/'provenance/terminal_v1.json'),n_cells=776270,n_sources=len(rows),n_programs=len(meta),all_score_values_range_and_shape_checked=total,all_group_mean_values_independently_recomputed=checked_means,native_manual_score_values=manual,max_native_manual_error=maxerr,quantile_files_not_independently_recomputed=True,all_values_independently_rescored=False,patient_material_independence_not_certified=True,scientific_module_eligibility_unchanged=True,wall_seconds=time.monotonic()-started))
 save(A/'provenance/terminal_v1.json',dict(status='UCELL_NUMERICAL_READBACK_ACCEPTED_REQUIRES_BIOLOGICAL_AND_READOUT_REVIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(A/'outputs/acceptance_v1.json')))

if __name__=='__main__':
 A=Path(sys.argv[1]);R=Path(sys.argv[2]);rc=0
 try:main(A,R)
 except BaseException:
  rc=1;save(A/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()))
 save(A/'provenance/supervisor_terminal_v1.json',dict(actual_child_exit_code=rc,science_terminal_exists=(A/'provenance/terminal_v1.json').exists(),completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()));sys.exit(rc)
