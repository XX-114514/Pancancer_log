"""Descriptive native endpoint check; no new inference or identity relabeling."""
from pathlib import Path
from datetime import datetime, timezone
import os,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8388608),b''):h.update(block)
 return h.hexdigest()

def main():
 start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';A=P/'fullgene_receiver_target_association_v1';O=P/'ADAM12_native_endpoint_probe_v1';assert not O.exists()
 plan=json.loads(Path(os.environ['PANCANCER_ENDPOINT_PLAN']).read_text());assert sha(A/'inventory.json')==plan['input_inventory_sha256']
 rec=next(r for r in json.loads((P/'supported_LR_native_cell_panel_v1/inventory.json').read_text())['source_records'] if r['gse_id']=='GSE179994');assert sha(rec['native_H5'])==rec['native_H5_sha256'] and sha(rec['selected_cache'])==rec['cache_sha256']
 z=dict(np.load(rec['selected_cache'],allow_pickle=False));panel=plan['gene_panel']
 with h5py.File(rec['native_H5'],'r') as h:
  genes=pd.Index(h['var'][h['var'].attrs['_index']].asstr()[:]);x=csr_matrix((h['X/data'][:],h['X/indices'][:],h['X/indptr'][:]),shape=tuple(h['X'].attrs['shape']))[z['source_row']]
  obscols=list(h['obs'].keys())
 assert np.array_equal(x[:,genes.get_indexer(z['genes'])].toarray(),z['raw_counts']);assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(),z['full_RNA_total'])
 available=[g for g in panel if g in genes];raw=x[:,genes.get_indexer(available)].toarray();gi={g:i for i,g in enumerate(available)};assert np.all(raw>=0) and np.all(raw==np.floor(raw));fam=z['families'][z['family_code']];ss=pd.Index(z['samples']);sums=[];lineage=[]
 with np.load(P/'fullgene_receiver_target_readiness_v1/private_source_unit_arrays/GSE179994.npz',allow_pickle=False) as b:
  bn=list(b['nodes']);bg=pd.Index(b['genes']);bc=b['positive_cell_counts'].reshape(len(bn),len(ss),len(bg))
  for node in bn:
   for g in available:
    actual=np.bincount(z['sample_code'][(z['subCluster']==node)&(raw[:,gi[g]]>0)],minlength=len(ss));assert np.array_equal(actual,bc[bn.index(node),:,bg.get_loc(g)])
 for tid,node in [('NODE0001','T_cell_IL7R_28'),('NODE0002','T_cell_IL7R_36')]:
  d=dict(np.load(A/'private_designs'/(tid+'.npz'),allow_pickle=False));codes=ss.get_indexer(d['sample_keys']);assert (codes>=0).all();sm=(fam=='T_cell_CCL5') & np.isin(z['sample_code'],codes);rm=(z['subCluster']==node)&np.isin(z['sample_code'],codes);lp=sm&(raw[:,gi['ADAM12']]>0);rp=rm&(raw[:,gi['ITGB1']]>0)
  by=[]
  for c in codes:
   s=sm&(z['sample_code']==c);r=rm&(z['sample_code']==c);a=int((lp&(z['sample_code']==c)).sum());b1=int((rp&(z['sample_code']==c)).sum());joint=int((r&(raw[:,gi['ITGA9']]>0)&(raw[:,gi['ITGB1']]>0)).sum()) if 'ITGA9' in gi else None;by.append((a,b1,int(s.sum()),int(r.sum()),joint))
  counts=np.asarray([r[:4] for r in by]);summary=dict(test_id=tid,receiver_original_node=node,n_paired_libraries=len(codes),sender_cells=int(sm.sum()),ADAM12_positive_sender_cells=int(lp.sum()),ADAM12_sender_positive_fraction=float(lp.sum()/sm.sum()),receiver_cells=int(rm.sum()),ITGB1_positive_receiver_cells=int(rp.sum()),ITGA9_positive_receiver_cells=int((rm&(raw[:,gi['ITGA9']]>0)).sum()) if 'ITGA9' in gi else None,ITGA9_ITGB1_same_cell_positive=int(sum(v[4] for v in by)) if 'ITGA9' in gi else None,libraries_both_L_R_positive_ge={str(k):int(((counts[:,0]>=k)&(counts[:,1]>=k)).sum()) for k in [1,5,10,20]},ADAM12_positive_sender_per_library_quantiles=np.quantile(counts[:,0],[0,.25,.5,.75,1]).tolist(),ADAM12_sender_positive_fraction_library_quantiles=np.quantile(counts[:,0]/counts[:,2],[0,.25,.5,.75,1]).tolist());sums.append(summary)
  for subset,mask in [('all_sender',sm),('ADAM12_positive_sender',lp),('ADAM12_zero_sender',sm&~lp),('all_receiver',rm)]:
   for g in available:
    lineage.append(dict(test_id=tid,receiver_original_node=node,subset=subset,gene=g,n_cells=int(mask.sum()),positive_cells=int((mask&(raw[:,gi[g]]>0)).sum()),positive_fraction=float((mask&(raw[:,gi[g]]>0)).sum()/mask.sum()) if mask.sum() else None))
 O.mkdir();pd.DataFrame(lineage).to_csv(O/'source_endpoint_marker_detection.tsv',sep='\t',index=False,na_rep='NA')
 out=dict(schema='ADAM12_native_endpoint_descriptive_probe_v1',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),status='COMPLETE_NATIVE_PANEL_REPLAY_AND_FULLGENE_POSITIVE_COUNTS_MATCHED',script_sha256=sha(__file__),plan_sha256=sha(os.environ['PANCANCER_ENDPOINT_PLAN']),native_input_sha256=rec['native_H5_sha256'],selected_cache_sha256=rec['cache_sha256'],available_genes=available,missing_genes=[g for g in panel if g not in genes],native_obs_column_names=obscols,summary=sums,audit='Native selected panel and full RNA totals exact replay; every original-node/sample/gene positive count matches accepted independent fullgene arrays',no_new_statistical_tests=True,no_identity_changes=True,independent_biological_validation=False,mechanism_verified=False,outputs={'source_endpoint_marker_detection.tsv':sha(O/'source_endpoint_marker_detection.tsv')})
 (O/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
