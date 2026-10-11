from pathlib import Path
import importlib.util,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
S=Path('/tmp/fullgene_receiver_target_readiness_v3.py');spec=importlib.util.spec_from_file_location('fixture_producer_v5',S);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
W=Path('/tmp/fullgene_original_node_fixture_v5');assert not W.exists();W.mkdir()
for n in ['private_source_unit_arrays','private_source_unit_catalogs']:(W/n).mkdir()
raw=np.array([[9,1,0],[2,0,2],[0,3,1],[3,4,0],[1,0,1],[0,7,0],[4,1,0]],np.int32);x=csr_matrix(raw);h5=W/'source.h5ad'
with h5py.File(h5,'w') as f:
 v=f.create_group('var');v.attrs['_index']='gene';v.create_dataset('gene',data=np.asarray(['A','B','R'],dtype=h5py.string_dtype()))
 a=f.create_group('X');a.attrs['shape']=x.shape;a.attrs['encoding-type']='csr_matrix'
 for key,val in [('data',x.data),('indices',x.indices),('indptr',x.indptr)]:a.create_dataset(key,data=val)
rows=np.array([1,2,4,6]);subset=raw[rows];labels=np.array(['T_cell_IL7R_28','T_cell_IL7R_28','T_cell_IL7R_36','T_cell_IL7R_36']);sc=np.array([0,0,1,0]);cache=W/'cache.npz'
np.savez_compressed(cache,subCluster=labels,samples=np.array(['fixture_unit_A','fixture_unit_B']),family_code=np.zeros(4,int),families=np.array(['T_cell_IL7R']),sample_code=sc,source_row=rows,genes=np.array(['A','R']),raw_counts=subset[:,[0,2]],full_RNA_total=subset.sum(axis=1),native_detected_genes=(subset>0).sum(axis=1))
rec=dict(selected_cache=str(cache),native_H5=str(h5),cache_sha256=m.sha(cache),native_H5_sha256=m.sha(h5),native_source_cells=7,selected_cells=4,gse_id='SYNTHETIC_FIXTURE')
sm,genes,bg=m.aggregate_source(rec,None,W);b=dict(np.load(W/'private_source_unit_arrays/SYNTHETIC_FIXTURE.npz',allow_pickle=False));assert sm['selected_cells']==4 and sm['unit_rows']==4
assert b['n_cells'].tolist()==[2,0,1,1] and np.isnan(b['mean_CP10K'][1]).all()
for ni,node in enumerate(b['nodes']):
 for si in range(2):
  ii=ni*2+si;mask=(labels==node)&(sc==si);v=subset[mask];n=len(v)
  assert np.array_equal(b['sum_raw_counts'][ii],v.sum(axis=0)) and np.array_equal(b['positive_cell_counts'][ii],(v>0).sum(axis=0))
  if n:
   norm=v/v.sum(axis=1)[:,None]*10000
   assert np.allclose(b['mean_CP10K'][ii],norm.mean(axis=0)) and np.allclose(b['mean_log1p_CP10K'][ii],np.log1p(norm).mean(axis=0))
cat=pd.read_csv(W/'private_source_unit_catalogs/SYNTHETIC_FIXTURE.tsv',sep='\t');assert not cat['background_diagnostics_eligible'].any() and cat['detected_genes1pct_n20'].isna().all() and cat['detected_genes10pct_n20'].isna().all()
assert np.isclose(b['mean_CP10K'][0,0],2500) and b['mean_CP10K'][0,1]==3750
lr=pd.DataFrame({'from':['L','CCL5','IL7R'],'to':['R','R','R'],'database':['fixture']*3,'source':['fixture']*3});tg=pd.DataFrame({'L':[1,2,2,9]},index=['A','Z','Y','R']);tf=pd.DataFrame({'L':[1,0]},index=['TF1','TF2'])
for lig,rec,expected in [('L','R','EXACT_SINGLETON_LR_PRIOR'),('CCL5_IL7R','R','COMPLEX_NOT_ASSEMBLED_BY_FLAT_PRIOR'),('UNKNOWN','R','NO_EXACT_LR_PRIOR')]:
 row=type('Row',(),dict(ligand=lig,receptor=rec));v=m.prior_fields(row,lr,tg,tf);assert v['LR_status']==expected and not v['receptor_specific_targets_available']
 if lig=='L':assert v['top200_after_exclusions']=='Y;Z;A' and v['nonzero_ligand_tf_prior_row_symbols']==1
 if lig=='UNKNOWN':assert not v['exact_ligand_target_column'] and v['nonzero_target_genes']==0
row=type('Row',(),dict(model_key='synthetic_key',gse_id='synthetic_discovery',ligand='L',receptor='R',receiver_family='fixture_receiver',LR_status='EXACT_SINGLETON_LR_PRIOR',exact_ligand_target_column=True,top200_after_exclusions='A;B'))
v=m.coverage_fields(row,'SYNTHETIC_SOURCE',{'A','B'},{});assert v['target_coverage_status']=='RECEIVER_FAMILY_NOT_CAPTURED' and v['targets_receiver_detected1pct'] is None and v['receiver_background_genes1pct'] is None and v['targets_in_native_axis']==2
v=m.coverage_fields(row,'SYNTHETIC_SOURCE',{'A','B'},{'fixture_receiver':dict(n_cells=20,genes1={'A'},genes10=set())});assert v['targets_receiver_detected1pct']==1 and v['targets_receiver_detected10pct']==0
row.exact_ligand_target_column=False;row.top200_after_exclusions='';v=m.coverage_fields(row,'SYNTHETIC_SOURCE',{'A','B'},{'fixture_receiver':dict(n_cells=20,genes1={'A'},genes10=set())});assert v['targets_in_native_axis'] is None and v['targets_receiver_detected1pct'] is None and v['receiver_background_genes1pct']==1
print(json.dumps(dict(status='PASS_RAW_LOG_FULL_DENOMINATOR_EMPTY_UNIT_ORIGINAL_NUMBERED_IDENTITIES_PRIOR_DIRECTION_COMPLEX_MISSING_TIES_AND_ABSENT_RECEIVER_NA',producer_sha256=m.sha(S))))
