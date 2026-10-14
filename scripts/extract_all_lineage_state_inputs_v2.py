"""Limited original-major sparse discovery counts, native total RNA denominators."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,re,subprocess,gc
import numpy as np,pandas as pd,anndata as ad,zarr
from scipy.sparse import csr_matrix,vstack
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def free():
 t=subprocess.check_output(['/usr/bin/quota','-g'],text=True);m=re.search(r'^\s*(\d+)\s+0\s+(\d+)\s+\d+\s+0\s+0\s*$',t,re.M);assert m
 return (int(m[2])-int(m[1]))*1024
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=Path(os.environ['PANCANCER_ALL_LINEAGE_INPUT_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256']
D=R/'all_lineage_state_discovery_row_selection_20261014_v1';A=R/'all_lineage_state_discovery_row_audit_20261014_v1/independent_audit.json'
assert sha(D/'inventory.json')==plan['selection_inventory_sha256'] and sha(A)==plan['selection_audit_sha256']
assert json.loads(A.read_text())['status'].startswith('PASS_ALL100766')
assert free()>=(plan['reserve_GiB']+plan['stage_budget_GiB'])*2**30
U=AR/'20261010_uniform_expression_annotation_v1';assert sha(U/'input_contract.json')==plan['input_contract_sha256']
c=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(c['core_zarr'],mode='r');x=g['X'];shape=c['raw_shape'];genes=json.loads((D/'fixed_common_training_genes_by_major.json').read_text())
md=pd.read_csv(D/'private_discovery_metadata.tsv.gz',sep='\t',keep_default_na=False);z=np.load(D/'private_discovery_row_indices.npz');assert np.array_equal(md.core_position,z['canonical_global_rows'])
O=R/'all_lineage_state_singlecell_inputs_20261014_v2';assert not O.exists();O.mkdir();(O/'private_major_counts').mkdir();start=datetime.now(timezone.utc).isoformat();summaries=[];files=[]
for mi,major in enumerate(sorted(genes)):
 assert free()>=(plan['reserve_GiB']+1)*2**30
 frame=md[md.major.eq(major)].copy().reset_index(drop=True);rows=frame.core_position.to_numpy();assert (np.diff(rows)>0).all();common=genes[major];mapping=None;allowed={}
 for sid in sorted(frame.gse_id.unique()):
  gt=pd.read_csv(U/'by_gse'/sid/'genes.tsv',sep='\t').set_index('gene');assert gt.index.is_unique
  cols=gt.loc[common,'core_column'].to_numpy();assert mapping is None or np.array_equal(mapping,cols);mapping=cols
  mask=np.zeros(shape[1],bool);mask[gt.core_column.to_numpy()]=True;allowed[sid]=mask
 chunks=[];totals=[];detect=[];begin=0
 while begin<len(rows):
  end=min(begin+plan['selected_cell_chunk'],len(rows));end=min(end,int(np.searchsorted(rows,rows[begin]+plan['max_core_row_span'],side='left')));assert end>begin
  rr=rows[begin:end];lo=int(rr[0]);hi=int(rr[-1])+1;p=x['indptr'][lo:hi+1].astype(np.int64);a,b=int(p[0]),int(p[-1]);full=csr_matrix((x['data'][a:b],x['indices'][a:b],p-a),shape=(hi-lo,shape[1]))[rr-lo].tocsr()
  assert np.isfinite(full.data).all() and (full.data>=0).all() and np.array_equal(full.data,np.floor(full.data)) and full.data.max(initial=0)<np.iinfo(np.int32).max
  full.sort_indices();assert full.has_canonical_format
  for sid in sorted(set(frame.gse_id.iloc[begin:end])):
   sm=frame.gse_id.iloc[begin:end].eq(sid).to_numpy();subset=full[sm];assert allowed[sid][subset.indices].all(),sid
  total=np.asarray(full.sum(axis=1)).ravel();det=np.asarray((full>0).sum(axis=1)).ravel();assert (total>0).all()
  chunks.append(full[:,mapping].astype(np.int32));totals.extend(total.tolist());detect.extend(det.tolist());begin=end
 raw=vstack(chunks,format='csr');del chunks
 frame['full_RNA_total']=np.asarray(totals,np.int64);frame['full_detected_genes']=np.asarray(detect,np.int32)
 key=g['obs'].attrs['_index'];frame.index=pd.Index(g['obs'][key].get_orthogonal_selection((rows,)).astype(str),name='cell_id');assert frame.index.is_unique
 frame['analysis_block']=frame.effective_resampling_block;frame['patient_identity_reviewed']=frame.unit_kind.eq('reviewed_patient')
 obj=ad.AnnData(X=raw,obs=frame,var=pd.DataFrame(index=pd.Index(common,name='gene')));obj.uns['normalization_contract']='Raw original counts; full_RNA_total native denominator, never common-axis sum; source-specific absence not zero-imputed'
 obj.uns['original_major_frozen']=True;name='major_%02d.h5ad'%mi;path=O/'private_major_counts'/name;obj.write_h5ad(path,compression='gzip',compression_opts=4)
 r=dict(major=major,cells=obj.n_obs,genes=obj.n_vars,sources=frame.gse_id.nunique(),libraries=frame.sample_uid.nunique(),inherited_blocks=frame.analysis_block.nunique(),native_nnz_bound=int(sum(detect)),common_axis_nnz=int(raw.nnz),source_limited=frame.gse_id.nunique()<2,bytes=path.stat().st_size,sha256=sha(path),private_relative_file='private_major_counts/'+name)
 summaries.append(r);files.append(dict(major=major,path=str(path),sha256=r['sha256']))
 (O/'progress.json').write_text(json.dumps(dict(status='EXTRACTING_ORIGINAL_MAJOR_COUNTS_NOT_MODEL_FIT',started_UTC=start,completed_majors=len(summaries),summary=summaries),indent=2)+'\n');print(json.dumps(r),flush=True);del obj,raw;gc.collect()
assert sum(r['cells'] for r in summaries)==100766
pd.DataFrame([{k:v for k,v in r.items() if k!='private_relative_file'} for r in summaries]).to_csv(O/'all18_discovery_count_scope.tsv',sep='\t',index=False)
out=dict(status='EXTRACTED_ALL100766_COUNTS_PENDING_NATIVE_READBACK_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),selection_inventory_sha256=sha(D/'inventory.json'),selection_audit_sha256=sha(A),cells=100766,majors=18,native_denominators=True,canonical_counts_unchanged=True,original_labels_unchanged=True,formal_state_labels=False,files=files,public_scope_sha256=sha(O/'all18_discovery_count_scope.tsv'),free_GiB=free()/2**30)
(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='files'}),flush=True)
