"""Fresh RNA IFNGR1/2 same-cell support in exact frozen Mac/T scope.
RNA joint positivity does not establish surface protein complex or signaling.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,zarr
from scipy.sparse import load_npz
from anndata._io.specs import read_elem

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';U=AR/'20261010_uniform_expression_annotation_v1';B=AR/'20261010_expression_annotation_downstream_v1';O=P/'mac_T_IFNG_receptor_single_cell_v1';planp=Path(os.environ['PANCANCER_SINGLE_CELL_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists();assert sha(U/'input_contract.json')==plan['input_contract_sha256'] and sha(U/'all_node_number_catalog.tsv')==plan['original_catalogue_sha256'];currentp=B/'downstream_review_v2/number_to_annotation_dictionary.tsv';assert sha(currentp)==plan['current_catalogue_sha256'];current=pd.read_csv(currentp,sep='\t').set_index('node_id');original=pd.read_csv(U/'all_node_number_catalog.tsv',sep='\t');context_plan=json.loads(Path(os.environ['PANCANCER_MAC_T_POOL_PLAN']).read_text());samples=set()
 for cid in context_plan['direction_contexts']['T_to_Macrophage']:samples.update(np.load(P/'major_receiver_TF_preflight_v3/private_context_arrays'/(cid+'.npz'))['samples'])
 contract=json.loads((U/'input_contract.json').read_text());assert contract['raw_layer']=='X';g=zarr.open_group(contract['core_zarr'],mode='r');X=g['X'];assert X.attrs['encoding-type']=='csr_matrix';sc=pd.Categorical(read_elem(g['obs']['sample_uid']));nc=pd.Categorical(read_elem(g['obs']['expression_v1_node_id']));slookup={str(v):i for i,v in enumerate(sc.categories)};nlookup={str(v):i for i,v in enumerate(nc.categories)};keys=sc.codes.astype(np.int64)*len(nc.categories)+nc.codes;allowed=np.zeros((len(sc.categories),len(nc.categories)),bool);lookup={};entries=[];columns={};sample_source={};started=datetime.now(timezone.utc).isoformat()
 for sid in plan['sources']:
  root=U/'by_gse'/sid;manifest=json.loads((root/'manifest.json').read_text());assert sha(root/'manifest.json')==plan['source_manifest_sha256'][sid]
  for name in ['units.tsv','genes.tsv','positive_counts.npz']:assert sha(root/name)==manifest['output_sha256'][name]
  units=pd.read_csv(root/'units.tsv',sep='\t');genes=pd.read_csv(root/'genes.tsv',sep='\t');positive=load_npz(root/'positive_counts.npz');columns[sid]={gene:int(genes[genes.gene==gene].iloc[0].core_column) for gene in plan['genes']};gi={v:i for i,v in enumerate(genes.gene)};sel=units[units.candidate_eligible&units.sample_uid.isin(samples)].copy();sel['node_id']=original.node_id.to_numpy()[sel.node_index];sel['major']=current.loc[sel.node_id,'majorCluster'].to_numpy();sel['subCluster']=current.loc[sel.node_id,'subCluster'].to_numpy();sel=sel[sel.major.isin(['T cell','Macrophage'])]
  for r in sel.itertuples():
   key=slookup[r.sample_uid]*len(nc.categories)+nlookup[r.node_id];assert key not in lookup;lookup[key]=len(entries);allowed[slookup[r.sample_uid],nlookup[r.node_id]]=True;sample_source[slookup[r.sample_uid]]=sid;entries.append(dict(unit_index=len(entries),unit_key=key,gse_id=sid,sample_uid=r.sample_uid,node_id=r.node_id,subCluster=r.subCluster,major=r.major,n_cells=int(r.n_cells),expected_IFNG=int(positive[r.unit_row,gi['IFNG']]),expected_IFNGR1=int(positive[r.unit_row,gi['IFNGR1']]),expected_IFNGR2=int(positive[r.unit_row,gi['IFNGR2']])))
 valid=(sc.codes>=0)&(nc.codes>=0);mask=np.zeros(len(keys),bool);mask[valid]=allowed.ravel()[keys[valid]];uk,counts=np.unique(keys[mask],return_counts=True);assert set(uk)==set(lookup) and all(int(n)==entries[lookup[int(k)]]['n_cells'] for k,n in zip(uk,counts));assert int(mask.sum())==plan['expected_cells'] and len(entries)==plan['expected_units'];bits=np.zeros(len(mask),np.uint8);source_code=np.array([plan['sources'].index(sample_source[i]) if i in sample_source else -1 for i in range(len(sc.categories))]);all_columns=sorted(set(v for d in columns.values() for v in d.values()));O.mkdir();ranges=0
 for start in range(0,len(mask),plan['row_chunk']):
  end=min(len(mask),start+plan['row_chunk'])
  if not mask[start:end].any():continue
  pointer=X['indptr'][start:end+1].astype(np.int64);lo=int(pointer[0]);hi=int(pointer[-1]);index=X['indices'][lo:hi];value=X['data'][lo:hi];assert (value>=0).all();pos=np.flatnonzero(np.isin(index,all_columns)&(value>0));rows=np.searchsorted(pointer-lo,pos,side='right')-1+start;take=mask[rows];rows=rows[take];hit_columns=index[pos[take]]
  for si,sid in enumerate(plan['sources']):
   same=source_code[sc.codes[rows]]==si
   for j,gene in enumerate(plan['genes']):rr=rows[same&(hit_columns==columns[sid][gene])];np.bitwise_or.at(bits,rr,np.uint8(1<<j))
  ranges+=1
  if ranges%10==0:print(json.dumps(dict(scanned_row_ranges=ranges,global_row=end,selected_cells_so_far=int(mask[:end].sum()))),flush=True)
 selected=np.flatnonzero(mask);unit=np.array([lookup[int(k)] for k in keys[selected]],int);b=bits[selected];E=pd.DataFrame(entries)
 for j,gene in enumerate(plan['genes']):E['positive_'+gene]=np.bincount(unit,weights=((b&(1<<j))>0).astype(int),minlength=len(E)).astype(int);assert np.array_equal(E['positive_'+gene],E['expected_'+gene])
 E['positive_IFNGR1_IFNGR2']=np.bincount(unit,weights=((b&6)==6).astype(int),minlength=len(E)).astype(int);E.to_csv(O/'private_all1832_sample_node_RNA_counts.tsv',sep='\t',index=False);np.savez_compressed(O/'private_selected_cell_gene_bits.npz',global_rows=selected,gene_bits=b,unit_index=unit,genes=np.array(plan['genes']),source_column_maps_json=json.dumps(columns));major=E.groupby(['gse_id','sample_uid','major'],as_index=False)[['n_cells','positive_IFNG','positive_IFNGR1','positive_IFNGR2','positive_IFNGR1_IFNGR2']].sum();major.to_csv(O/'private_all424_sample_major_RNA_counts.tsv',sep='\t',index=False);summary=major.groupby(['gse_id','major'],as_index=False).agg(cells=('n_cells','sum'),libraries=('sample_uid','nunique'),IFNG_positive=('positive_IFNG','sum'),IFNGR1_positive=('positive_IFNGR1','sum'),IFNGR2_positive=('positive_IFNGR2','sum'),IFNGR1_IFNGR2_positive=('positive_IFNGR1_IFNGR2','sum'));fine=E.groupby(['gse_id','major','subCluster'],as_index=False).agg(cells=('n_cells','sum'),libraries=('sample_uid','nunique'),IFNG_positive=('positive_IFNG','sum'),IFNGR1_positive=('positive_IFNGR1','sum'),IFNGR2_positive=('positive_IFNGR2','sum'),IFNGR1_IFNGR2_positive=('positive_IFNGR1_IFNGR2','sum'))
 for df in [summary,fine]:
  for gene in plan['genes']+['IFNGR1_IFNGR2']:df[gene+'_positive_fraction']=df[gene+'_positive']/df.cells
 summary.to_csv(O/'source_major_RNA_support.tsv',sep='\t',index=False);fine.to_csv(O/'source_frozen_subcluster_RNA_support.tsv',sep='\t',index=False);j=dict(status='FRESH_SAME_CELL_RNA_COUNTS_PRODUCED_PENDING_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),selected_cells=len(selected),sample_node_units=len(E),sample_major_units=len(major),fresh_all3_marginal_counts_equal_native_cache=True,producer_sha256=sha(__file__),plan_sha256=sha(planp),no_surface_complex_or_communication_claim=True,outputs={str(f.relative_to(O)):sha(f) for f in O.iterdir() if f.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}),flush=True)
if __name__=='__main__':main()
