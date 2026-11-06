import json
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';a=r/'stage_a_full_axis';d=r/'stage_d_frozen_marker_evidence';c=r/'stage_c_marker_identity_freeze'
assert json.loads((c/'independent_freeze_audit.json').read_text())['status']=='PASS';assert (d/'complete.json').exists()
t=pd.read_pickle(c/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');q=pd.read_pickle(r/'stage_b2_refined_components/all_cells_refined_identity_evidence_PRIVATE.pkl.gz');assert t.index.equals(q.index);m=pd.read_csv(d/'frozen_naming_marker_support.tsv',sep='\t');allm=pd.read_csv(d/'frozen_label_marker_expression.tsv',sep='\t');summ=pd.read_csv(d/'frozen_label_identity_state_summary.tsv',sep='\t');assert not allm.duplicated(['final_annotation','gene']).any()
raw=np.load(a/'marker_counts_PRIVATE.npy',mmap_mode='r');tot=np.load(a/'full_raw_library_totals_PRIVATE.npy');genes=np.load(a/'gene_symbols.npy').tolist();used=sorted(m.naming_marker.unique());ii=[genes.index(g) for g in used];av=np.load(a/'source_gene_measured.npy')[:,ii];src=pd.Index(pd.read_csv(a/'source_axis.tsv',sep='\t').source.astype(str));si=src.get_indexer(q.gse_id.astype(str));labels=pd.Index(summ.final_annotation);take=t.annotation_frozen.to_numpy();pos=np.flatnonzero(take);cc=labels.get_indexer(t.loc[take,'final_annotation'].astype(str));K=len(labels);G=len(used);n=np.bincount(cc,minlength=K);den=np.zeros((K,G),np.int64);nz=np.zeros_like(den);cp=np.zeros((K,G));sourcecounts=np.zeros((K,len(src)),np.int64)
for left in range(0,len(pos),50000):
 right=min(left+50000,len(pos));ids=pos[left:right];cl=cc[left:right];x=np.asarray(raw[ids])[:,ii];measured=av[si[ids]]
 for j in range(G):
  den[:,j]+=np.bincount(cl,weights=measured[:,j],minlength=K).astype(np.int64);nz[:,j]+=np.bincount(cl,weights=x[:,j]>0,minlength=K).astype(np.int64);cp[:,j]+=np.bincount(cl,weights=10000.0*x[:,j]/np.maximum(tot[ids],1),minlength=K)
 for i in range(K):sourcecounts[i]+=np.bincount(si[ids[cl==i]],minlength=len(src))
errmax=0.
for z in m.to_dict('records'):
 i=labels.get_loc(z['final_annotation']);j=used.index(z['naming_marker']);assert z['n_frozen_cells']==n[i] and z['n_measured']==den[i,j];assert abs(z['coverage_fraction']-den[i,j]/n[i])<1e-12;assert den[i,j]>0 and z['pct_expressing_measured']>0
 err=abs(z['mean_CP10K_measured']-cp[i,j]/den[i,j]);errmax=max(errmax,err);assert err<1e-7 and abs(z['pct_expressing_measured']-nz[i,j]/den[i,j])<1e-12
 full=allm[allm.final_annotation.eq(z['final_annotation'])&allm.gene.eq(z['naming_marker'])].iloc[0];assert abs(full.mean_cp10k_measured-z['mean_CP10K_measured'])<1e-10 and full.n_expressing==nz[i,j]
assert np.array_equal(n,summ.n_frozen_cells.to_numpy()) and n.sum()==take.sum()
for i,z in enumerate(summ.to_dict('records')):
 assert int((sourcecounts[i]>0).sum())==z['n_sources'] and abs(sourcecounts[i].max()/sourcecounts[i].sum()-z['max_source_fraction'])<1e-12
assert allm.groupby('final_annotation').n_frozen_cells.nunique().eq(1).all() and allm.n_expressing.le(allm.n_measured).all() and allm.n_measured.le(allm.n_frozen_cells).all()
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_frozen_cells':int(take.sum()),'n_labels':K,'n_all_marker_records':len(allm),'n_naming_marker_records':len(m),'n_distinct_naming_genes':G,'all_naming_gene_counts_denominators_CP10K_and_source_concentration_independently_verified':True,'max_absolute_naming_marker_CP10K_error':errmax,'no_fine_subtype_or_function_inferred_from_names':True};(d/'independent_naming_marker_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
