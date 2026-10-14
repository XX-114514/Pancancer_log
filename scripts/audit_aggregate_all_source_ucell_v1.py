"""Full identity/support audit, prespecified manual scores, private sample-node means."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,zarr
from scipy.sparse import csr_matrix,vstack

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def tied_ranks(row,cap):
 out=np.full(len(row),float(cap));ix=np.flatnonzero(row>0);ix=ix[np.argsort(-row[ix],kind='stable')];a=0
 while a<len(ix):
  b=a+1
  while b<len(ix) and row[ix[b]]==row[ix[a]]:b+=1
  out[ix[a:b]]=min((a+1+b)/2,cap);a=b
 return out
def score(ranks,ix,cap):
 n=len(ix);return 1-(sum(float(ranks[i]) for i in ix)-n*(n+1)/2)/(n*cap-n*(n+1)/2)
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path('/tmp/all_source_ucell_score_audit_plan_v1.json');plan=json.loads(P.read_text());assert sha(__file__)==plan['auditor_sha256'];S=V/'all_source_native_UCell_projection_20261014_v1';iv=json.loads((S/'inventory.json').read_text());assert iv['status']=='COMPLETE_ALL4532933_NATIVE_UCELL_PENDING_INDEPENDENT_AUDIT';D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)==iv['dictionary_sha256'];records=json.loads(D.read_text())['programmes'];F=R/'all_CM_state_native_input_preflight_20261014_v1';fi=json.loads((F/'inventory.json').read_text());assert sha(F/'inventory.json')==iv['scope_inventory_sha256'];assert sha(F/'private_canonical_row_selection.npz')==fi['private_row_selection_sha256'];assert sha(F/'private_all_frozen_sample_node_units.tsv.gz')==fi['private_unit_table_sha256'];units=pd.read_csv(F/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False)
 with np.load(F/'private_canonical_row_selection.npz',allow_pickle=False) as z:expected_rows=z['global_rows'];expected_units=z['unit_index']
 seen=np.zeros(len(expected_rows),bool);U=AR/'20261010_uniform_expression_annotation_v1';assert sha(U/'input_contract.json')==plan['canonical_contract_sha256'];contract=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');x=g['X'];key=g['obs'].attrs['_index'];shape=contract['raw_shape'];src=pd.read_csv(F/'all_source_native_scope.tsv',sep='\t').set_index('gse_id');majors=sorted(units.major.unique());O=V/'all_source_native_UCell_audit_aggregation_20261014_v1';assert not O.exists();O.mkdir(mode=0o700);aggs={};completed=[];manual_count=0;quantity_count=0;maxerr=0.;rng=np.random.default_rng(20261014);gene_cache={}
 for major in majors:
  ui=units.loc[units.major.eq(major),'unit_index'].to_numpy(int);rr=[r for r in records if major in r['allowed_majors']];n,p=len(ui),len(rr);aggs[major]=dict(unit_index=ui,programmes=rr,sums=np.zeros((n,p)),valid=np.zeros((n,p),np.int64),support=np.zeros((n,p),np.int64),RNA_fraction_sum=np.zeros((n,p)),cells=np.zeros(n,np.int64))
 for file in iv['private_score_files']:
  p=S/file['relative_file'];assert sha(p)==file['sha256'];z=dict(np.load(p,allow_pickle=False));rr=z['global_rows'];uid=z['unit_index'];where=np.searchsorted(expected_rows,rr);assert (where<len(expected_rows)).all() and np.array_equal(expected_rows[where],rr) and np.array_equal(expected_units[where],uid) and not seen[where].any();seen[where]=True;md=units.iloc[uid];major=file['major'];sid=file['gse_id'];assert md.major.eq(major).all() and md.gse_id.eq(sid).all();assert np.array_equal(g['obs'][key].get_orthogonal_selection((rr,)).astype(str),z['cell_ids']);a=aggs[major];recs=a['programmes'];assert np.array_equal(z['programme_ids'],np.array([r['programme_id'] for r in recs]));n,pcount=len(rr),len(recs)
  if sid not in gene_cache:
   gp=U/'by_gse'/sid/'genes.tsv';assert sha(gp)==src.loc[sid,'native_gene_axis_sha256'];gene_cache[sid]=pd.read_csv(gp,sep='\t')
  gt=gene_cache[sid];genes=pd.Index(gt.gene.astype(str));cols=gt.core_column.to_numpy(int);allowed=np.zeros(shape[1],bool);allowed[cols]=True;primary=[];alternate=[];gate=[];agate=[];ar=[];ac=[]
  for j,r in enumerate(recs):
   ix=genes.get_indexer(r['primary_genes']);present=ix>=0;primary.append(ix[present]);gate.append(present.sum()>=10 and present.mean()>=.8);an=r['TPM_sensitivity_genes'];ai=genes.get_indexer(an);ap=ai>=0;alternate.append(ai[ap]);agate.append(bool(len(an) and ap.sum()>=10 and ap.mean()>=.8));ar.extend(cols[ix[present]].tolist());ac.extend([j]*int(present.sum()))
  assert np.array_equal(z['primary_gene_coverage_mask'],gate) and np.array_equal(z['TPM_gene_coverage_mask'],agate);A=csr_matrix((np.ones(len(ar),np.int64),(ar,ac)),shape=(shape[1],pcount));parts=[];start=0
  while start<n:
   end=min(start+500,n,int(np.searchsorted(rr,rr[start]+5000,side='left')));take=rr[start:end];lo=int(take[0]);hi=int(take[-1])+1;ip=np.asarray(x['indptr'][lo:hi+1],np.int64);aa,bb=int(ip[0]),int(ip[-1]);native=csr_matrix((x['data'][aa:bb],x['indices'][aa:bb],ip-aa),shape=(hi-lo,shape[1]))[take-lo].tocoo();assert allowed[native.col].all() and np.isfinite(native.data).all() and (native.data>=0).all() and np.array_equal(native.data,np.floor(native.data));parts.append(csr_matrix((native.data.astype(np.int64),(native.row,native.col)),shape=native.shape));start=end
  raw=vstack(parts,format='csr');del parts;tot=np.asarray(raw.sum(axis=1)).ravel();detall=np.asarray((raw>0).sum(axis=1)).ravel();assert np.array_equal(tot,z['full_RNA_total']) and np.array_equal(detall,z['full_detected_genes']);det=((raw>0).astype(np.int64)@A).toarray();frac=(raw@A).toarray()/tot[:,None];assert np.array_equal(det,z['detected_programme_genes']) and np.allclose(frac,z['programme_RNA_fraction'],atol=1e-7,rtol=1e-6);quantity_count+=2*n*pcount
  for cap in [1000,1500,3000]:
   y=z['scores_maxrank'+str(cap)];assert y.shape==(n,pcount) and np.array_equal(np.isfinite(y),np.broadcast_to(gate,(n,pcount))) and ((y[np.isfinite(y)]>=-1e-6)&(y[np.isfinite(y)]<=1+1e-6)).all()
  assert np.array_equal(np.isfinite(z['TPM_signature_scores_maxrank1500']),np.broadcast_to(agate,(n,pcount)))
  selected=np.unique(np.r_[0,n-1,rng.choice(n,min(plan['manual_cells_per_chunk']-2,n),replace=False)])
  for cell in selected:
   row=raw[cell,cols].toarray().ravel()
   for cap in [1000,1500,3000]:
    ranks=tied_ranks(row,cap)
    for j,ok in enumerate(gate):
     if ok:
      expected=score(ranks,primary[j],cap);err=abs(expected-float(z['scores_maxrank'+str(cap)][cell,j]));assert err<=1e-6;manual_count+=1;maxerr=max(maxerr,err)
   ranks=tied_ranks(row,1500)
   for j,ok in enumerate(agate):
    if ok:
     err=abs(score(ranks,alternate[j],1500)-float(z['TPM_signature_scores_maxrank1500'][cell,j]));assert err<=1e-6;manual_count+=1;maxerr=max(maxerr,err)
  # Two aggregation implementations on every chunk; no patient-independence assumption.
  local=np.searchsorted(a['unit_index'],uid);assert np.array_equal(a['unit_index'][local],uid);y=z['scores_maxrank1500'].astype(float);vals=np.nan_to_num(y);mask=np.isfinite(y);inc=csr_matrix((np.ones(n),(local,np.arange(n))),shape=(len(a['unit_index']),n));sums=inc@vals;direct=np.zeros_like(sums);np.add.at(direct,local,vals);assert np.allclose(sums,direct,atol=1e-10,rtol=1e-12);a['sums']+=sums;a['valid']+=np.asarray(inc@mask.astype(np.int64),np.int64);a['support']+=np.asarray(inc@((det>=5)&mask).astype(np.int64),np.int64);a['RNA_fraction_sum']+=inc@(frac*mask);a['cells']+=np.bincount(local,minlength=len(a['unit_index']));completed.append(dict(gse_id=sid,major=major,cells=n,manual_cells=len(selected)));(O/'progress.json').write_text(json.dumps(dict(chunks_checked=len(completed),cells_checked=int(seen.sum()),manual_scores_checked=manual_count,observed_UTC=datetime.now(timezone.utc).isoformat()),indent=2)+'\n');del raw,z,det,frac,A,inc,y;gc.collect()
 assert seen.all() and len(expected_rows)==4532933;aggregate_files=[]
 for i,major in enumerate(majors):
  a=aggs[major];assert np.array_equal(a['cells'],units.iloc[a['unit_index']].n_cells);means=np.divide(a['sums'],a['valid'],out=np.full(a['sums'].shape,np.nan),where=a['valid']>0);f=O/f'private_major_{i:02}_sample_node_programmes.npz';np.savez_compressed(f,unit_index=a['unit_index'],programme_ids=np.array([r['programme_id'] for r in a['programmes']]),cells=a['cells'],valid_score_cells=a['valid'],programme_mean=means,informative_gene_support_cells=a['support'],programme_RNA_fraction_sum=a['RNA_fraction_sum']);aggregate_files.append(dict(major=major,relative_file=f.name,sha256=sha(f),sample_node_units=len(a['unit_index']),programmes=len(a['programmes'])))
 out=dict(status='PASS_ALL4532933_IDENTITIES_SUPPORT_MASKS_AND_SAMPLED_MANUAL_UCELL_WITH_PRIVATE_UNIT_MEANS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),plan_sha256=sha(P),producer_inventory_sha256=sha(S/'inventory.json'),dictionary_sha256=sha(D),cells=len(expected_rows),chunks=len(completed),quantity_values_checked=quantity_count,manual_score_values_checked=manual_count,max_manual_score_error=maxerr,all_score_values_independently_rescored=False,all_sample_node_cells_match_frozen_counts=True,all_chunk_means_two_implementations_checked=True,private_aggregate_files=aggregate_files,formal_states=False,biological_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='private_aggregate_files'}),flush=True)
if __name__=='__main__':main()
