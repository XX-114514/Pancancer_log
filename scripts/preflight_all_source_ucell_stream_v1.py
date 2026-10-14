"""Read-only full native score storage/streaming feasibility; no new RNA copies."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,shutil
import numpy as np,pandas as pd,zarr

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path(os.environ['PANCANCER_UCELL_STREAM_PREFLIGHT_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];F=R/'all_CM_state_native_input_preflight_20261014_v1';iv=json.loads((F/'inventory.json').read_text());assert sha(F/'inventory.json')==plan['scope_inventory_sha256'];U=AR/'20261010_uniform_expression_annotation_v1';assert sha(U/'input_contract.json')==plan['canonical_contract_sha256'];contract=json.loads((U/'input_contract.json').read_text());assert contract['raw_layer']=='X'
for name,key in [('private_canonical_row_selection.npz','private_row_selection_sha256'),('private_all_frozen_sample_node_units.tsv.gz','private_unit_table_sha256')]:assert sha(F/name)==iv[key]
units=pd.read_csv(F/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False)
with np.load(F/'private_canonical_row_selection.npz',allow_pickle=False) as z:rows=z['global_rows'];uid=z['unit_index']
assert len(rows)==4532933 and np.all(np.diff(rows)>0) and len(uid)==len(rows);assert np.array_equal(units.unit_index,np.arange(len(units)));assert np.array_equal(np.bincount(uid,minlength=len(units)),units.n_cells)
g=zarr.open_group(contract['core_zarr'],mode='r');indptr=np.asarray(g['X']['indptr'][:],dtype=np.int64);assert len(indptr)==contract['raw_shape'][0]+1;nnz=indptr[rows+1]-indptr[rows];source=units.gse_id.to_numpy()[uid];major=units.major.to_numpy()[uid];sources=pd.read_csv(F/'all_source_native_scope.tsv',sep='\t').set_index('gse_id');results=[]
for sid in sorted(set(source)):
 p=U/'by_gse'/sid/'genes.tsv';assert sha(p)==sources.loc[sid,'native_gene_axis_sha256'];genes=pd.read_csv(p,sep='\t');assert genes.gene.is_unique and genes.core_column.is_unique;assert genes.core_column.between(0,contract['raw_shape'][1]-1).all()
 for m in sorted(set(major[source==sid])):
  ix=np.flatnonzero((source==sid)&(major==m));n=len(ix);nchunk=(n+plan['output_cells_per_chunk']-1)//plan['output_cells_per_chunk'];rowbytes=int(nnz[ix].sum())*8+(n+1)*8
  results.append(dict(gse_id=sid,major=m,cells=n,native_genes=len(genes),canonical_nnz_upper_bound=int(nnz[ix].sum()),source_major_CSR_bytes_upper_bound=rowbytes,output_chunks=nchunk,max_individual_row_nnz=int(nnz[ix].max()),primary_signature_capacity=plan['programme_capacity_per_major'],score_sidecars_bytes_upper_bound=n*plan['programme_capacity_per_major']*plan['bytes_per_cell_programme'],native_dense_rank_batch_bytes=len(genes)*plan['ucell_rank_batch_cells']*8*plan['dense_rank_workspace_multiplier']))
t=pd.DataFrame(results);bytes_out=int(t.score_sidecars_bytes_upper_bound.sum());required=bytes_out+plan['metadata_and_aggregate_budget_GiB']*2**30;free=shutil.disk_usage(V).free;fits=required<=plan['total_output_budget_GiB']*2**30 and free-required>=plan['reserve_free_GiB']*2**30;assert t.cells.sum()==4532933 and t.gse_id.nunique()==52 and t.major.nunique()==18
O=V/'all_source_ucell_stream_preflight_20261014_v1';assert not O.exists();O.mkdir();t.to_csv(O/'all_source_major_stream_budget.tsv',sep='\t',index=False);out=dict(status='PASS_NATIVE_STREAM_CAPACITY_PENDING_COMPLETE_DICTIONARY' if fits else 'INSUFFICIENT_NATIVE_STREAM_CAPACITY',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),sources=52,majors=18,cells=len(rows),source_major_groups=len(t),output_chunks=int(t.output_chunks.sum()),max_native_genes=int(t.native_genes.max()),max_dense_rank_workspace_bytes=int(t.native_dense_rank_batch_bytes.max()),score_sidecars_bytes_upper_bound=bytes_out,required_bytes_with_metadata=required,volume_free_bytes=free,capacity_pass=bool(fits),programme_capacity_per_major=plan['programme_capacity_per_major'],capacity_must_be_rechecked_with_actual_dictionary=True,native_RNA_copied=False,actual_scoring_started=False,formal_states=False,output_sha256={'all_source_major_stream_budget.tsv':sha(O/'all_source_major_stream_budget.tsv')});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);assert fits
