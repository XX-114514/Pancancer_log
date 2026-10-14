"""Full native UCell projection, original cell identities, bounded source-major chunks."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,shutil,gc,importlib.metadata
import numpy as np,pandas as pd,zarr,anndata as ad
from scipy.sparse import csr_matrix,vstack
import pyucell as uc

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);P=Path(os.environ['PANCANCER_UCELL_FULL_PROJECTION_PLAN']);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];assert importlib.metadata.version('pyucell')=='0.8.0'
 D=V/plan['dictionary_directory'];di=json.loads((D/'inventory.json').read_text());dp=D/'frozen_programme_dictionary.json';assert sha(dp)==di['dictionary_sha256'];d=json.loads(dp.read_text());assert d['schema']==plan['dictionary_schema'] and len(d['input_scopes'])==24
 # Independent dictionary acceptance must be available before any score output is created.
 audit=json.loads((V/plan['dictionary_acceptance_relative_file']).read_text());assert audit['status']=='ACCEPTED_ALL18_DICTIONARY_FOR_NATIVE_PROJECTION_NOT_FORMAL_STATES' and audit['dictionary_sha256']==sha(dp)
 recs=d['programmes'];assert len({r['programme_id'] for r in recs})==len(recs);F=R/'all_CM_state_native_input_preflight_20261014_v1';iv=json.loads((F/'inventory.json').read_text());assert sha(F/'inventory.json')==plan['scope_inventory_sha256'];units=pd.read_csv(F/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False);assert sha(F/'private_all_frozen_sample_node_units.tsv.gz')==iv['private_unit_table_sha256'];assert sha(F/'private_canonical_row_selection.npz')==iv['private_row_selection_sha256']
 with np.load(F/'private_canonical_row_selection.npz',allow_pickle=False) as z:rows=z['global_rows'];uid=z['unit_index']
 assert len(rows)==4532933 and np.all(np.diff(rows)>0) and np.array_equal(np.bincount(uid,minlength=len(units)),units.n_cells)
 U=AR/'20261010_uniform_expression_annotation_v1';assert sha(U/'input_contract.json')==plan['canonical_contract_sha256'];contract=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');x=g['X'];shape=contract['raw_shape'];key=g['obs'].attrs['_index'];unit_groups,group_names=pd.factorize(units.gse_id.astype(str)+'::'+units.major.astype(str),sort=True);row_groups=unit_groups[uid];order=np.argsort(row_groups,kind='stable');offset=np.r_[0,np.cumsum(np.bincount(row_groups))];group_lookup={str(n):i for i,n in enumerate(group_names)};srcinv=pd.read_csv(F/'all_source_native_scope.tsv',sep='\t').set_index('gse_id');majors=sorted(units.major.unique());bymajor={m:[r for r in recs if m in r['allowed_majors']] for m in majors};assert all(0<len(v)<=plan['programme_capacity_per_major'] for v in bymajor.values());expected=sum(int(units.loc[units.major.eq(m),'n_cells'].sum())*len(bymajor[m])*22 for m in majors)+plan['metadata_and_aggregate_budget_GiB']*2**30;assert expected<=plan['total_output_budget_GiB']*2**30 and shutil.disk_usage(V).free>=expected+plan['reserve_free_GiB']*2**30
 O=V/'all_source_native_UCell_projection_20261014_v1';assert not O.exists();O.mkdir(mode=0o700);files=[];coverage=[];seen=np.zeros(len(rows),bool);started=datetime.now(timezone.utc).isoformat()
 for sid in sorted(units.gse_id.unique()):
  gf=U/'by_gse'/sid/'genes.tsv';assert sha(gf)==srcinv.loc[sid,'native_gene_axis_sha256'];gt=pd.read_csv(gf,sep='\t');genes=pd.Index(gt.gene.astype(str));cols=gt.core_column.to_numpy(int);allowed=np.zeros(shape[1],bool);allowed[cols]=True
  for m in majors:
   group_id=group_lookup.get(sid+'::'+m);selected=order[offset[group_id]:offset[group_id+1]] if group_id is not None else np.array([],dtype=int)
   if not len(selected):continue
   records=bymajor[m];ids=np.array([r['programme_id'] for r in records]);sig={};alt={};ixs=[];gates=[];altgates=[]
   for r in records:
    names=r['primary_genes'];assert 0<len(names)<=50 and len(names)==len(set(names));ix=genes.get_indexer(names);present=ix>=0;ok=present.sum()>=plan['min_measured_genes'] and present.mean()>=plan['min_measured_fraction'];ixs.append(ix[present]);gates.append(ok)
    if ok:sig[r['programme_id']]=np.array(names)[present].tolist()
    an=r['TPM_sensitivity_genes'];ai=genes.get_indexer(an);ap=ai>=0;aok=bool(len(an) and ap.sum()>=plan['min_measured_genes'] and ap.mean()>=plan['min_measured_fraction']);altgates.append(aok)
    if aok:alt[r['programme_id']]=np.array(an)[ap].tolist()
    coverage.append(dict(gse_id=sid,major=m,programme_id=r['programme_id'],required_genes=len(names),measured_genes=int(present.sum()),primary_scorable=bool(ok),TPM_required_genes=len(an),TPM_measured_genes=int(ap.sum()),TPM_scorable=aok,missing_genes=';'.join(np.array(names)[~present])))
   folder=O/sid/f'major_{majors.index(m):02}';folder.mkdir(parents=True);pids=len(records)
   for chunk_start in range(0,len(selected),plan['output_cells_per_chunk']):
    assert shutil.disk_usage(V).free>=(plan['reserve_free_GiB']+1)*2**30
    selected_chunk=selected[chunk_start:chunk_start+plan['output_cells_per_chunk']];rr=rows[selected_chunk];parts=[];start=0
    while start<len(rr):
     end=min(start+plan['read_cells_per_chunk'],len(rr),int(np.searchsorted(rr,rr[start]+plan['max_core_row_span'],side='left')));assert end>start;take=rr[start:end];lo=int(take[0]);hi=int(take[-1])+1;ip=np.asarray(x['indptr'][lo:hi+1],np.int64);a,b=int(ip[0]),int(ip[-1]);raw=csr_matrix((x['data'][a:b],x['indices'][a:b],ip-a),shape=(hi-lo,shape[1]))[take-lo].tocsr();raw.sort_indices();assert raw.has_canonical_format and allowed[raw.indices].all();assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.array_equal(raw.data,np.floor(raw.data));parts.append(raw[:,cols].astype(np.float64));start=end
    counts=vstack(parts,format='csr');del parts;tot=np.asarray(counts.sum(axis=1)).ravel();detected=np.asarray((counts>0).sum(axis=1)).ravel();assert (tot>0).all();cellids=g['obs'][key].get_orthogonal_selection((rr,)).astype(str);assert len(set(cellids))==len(rr);a=ad.AnnData(X=counts,obs=pd.DataFrame(index=cellids),var=pd.DataFrame(index=genes));det=np.zeros((len(rr),pids),np.uint8);frac=np.zeros((len(rr),pids),np.float32)
    for j,ix in enumerate(ixs):det[:,j]=np.asarray((counts[:,ix]>0).sum(axis=1)).ravel();frac[:,j]=np.asarray(counts[:,ix].sum(axis=1)).ravel()/tot
    scores={}
    for rank in plan['max_rank_values']:
     if sig:uc.compute_ucell_scores(a,signatures=sig,max_rank=rank,ties_method='average',missing_genes='skip',chunk_size=plan['ucell_rank_batch_cells'],n_jobs=1,device='cpu',suffix='_UCell')
     out=np.full((len(rr),pids),np.nan,np.float32)
     for j,pid in enumerate(ids):
      if gates[j]:out[:,j]=a.obs[pid+'_UCell'].to_numpy()
     scores['scores_maxrank'+str(rank)]=out
    if alt:uc.compute_ucell_scores(a,signatures=alt,max_rank=1500,ties_method='average',missing_genes='skip',chunk_size=plan['ucell_rank_batch_cells'],n_jobs=1,device='cpu',suffix='_TPM_UCell')
    out=np.full((len(rr),pids),np.nan,np.float32)
    for j,pid in enumerate(ids):
     if altgates[j]:out[:,j]=a.obs[pid+'_TPM_UCell'].to_numpy()
    scores['TPM_signature_scores_maxrank1500']=out
    for v in scores.values():assert np.isfinite(v).sum()==np.sum(~np.isnan(v)) and ((v[np.isfinite(v)]>=-1e-6)&(v[np.isfinite(v)]<=1+1e-6)).all()
    path=folder/f'chunk_{chunk_start:08}.npz';np.savez_compressed(path,cell_ids=cellids,global_rows=rr,unit_index=uid[selected_chunk],programme_ids=ids,full_RNA_total=tot.astype(np.int64),full_detected_genes=detected.astype(np.int32),detected_programme_genes=det,programme_RNA_fraction=frac,primary_gene_coverage_mask=np.array(gates),TPM_gene_coverage_mask=np.array(altgates),**scores);assert not seen[selected_chunk].any();seen[selected_chunk]=True;files.append(dict(gse_id=sid,major=m,cells=len(rr),programmes=pids,relative_file=str(path.relative_to(O)),sha256=sha(path)));(O/'progress.json').write_text(json.dumps(dict(status='SCORING_FROZEN_ALL18_PROGRAMMES',observed_UTC=datetime.now(timezone.utc).isoformat(),cells_completed=int(seen.sum()),chunks_completed=len(files)),indent=2)+'\n');del a,counts,raw,scores,det,frac,out;gc.collect()
 assert seen.all();pd.DataFrame(coverage).to_csv(O/'all_source_major_programme_gene_coverage.tsv',sep='\t',index=False);out=dict(status='COMPLETE_ALL4532933_NATIVE_UCELL_PENDING_INDEPENDENT_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),dictionary_sha256=sha(dp),scope_inventory_sha256=sha(F/'inventory.json'),cells=len(rows),sources=52,majors=18,private_score_files=files,public_coverage_sha256=sha(O/'all_source_major_programme_gene_coverage.tsv'),original_identities_unchanged=True,formal_states=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='private_score_files'}),flush=True)
if __name__=='__main__':main()
