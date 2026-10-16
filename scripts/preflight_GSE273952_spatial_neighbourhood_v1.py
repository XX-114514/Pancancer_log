"""Spatial input geometry and historical partial reference audit; no outcome inference."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);O=V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def strs(x):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in x])
def field(obs,key):
 z=obs[key]
 if isinstance(z,h5py.Group):
  c=z['codes'][:];cats=strs(z['categories'][:]);return np.array([cats[k] if k>=0 else '' for k in c])
 return strs(z[:])
I=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';iv=json.loads((I/'inventory.json').read_text());old=AR.parent.parent/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';path=Path(json.loads((old/'outputs/spatial/GSE273952_spatial_contract_audit.json').read_text())['input_h5ad']);assert sha(path)==iv['input_file_sha256'];A=V/'existing_GSE273952_Visium_input_audit_20261016_v1';ai=json.loads((A/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL71398');Q=np.load(A/'private_all_spot_derived_QC_and_grid.npz',allow_pickle=False);codes=Q['library_codes'];rr=Q['array_row'].astype(int);cc=Q['array_col'].astype(int);n=len(rr);assert n==71398
with h5py.File(path,'r') as f:
 obs=f['obs'];rawids=strs(obs[obs.attrs['_index']][:]);tls=field(obs,'TLS_annotation_TLS_annotation');tls2=field(obs,'TLS_annotation_TLS_2_cat');tumor=field(obs,'TumoralZone_annotation_Tumor_zone');x=np.asarray(obs['pxl_col_in_fullres'][:],float);y=np.asarray(obs['pxl_row_in_fullres'][:],float)
assert len(set(rawids))==n;fold=np.array(['folded' in (a+' '+b+' '+c).lower() for a,b,c in zip(tls,tls2,tumor)]);keep=~fold
refdir=old/'outputs/spatial/cell2location/final';auditpath=refdir/'GSE273952_cell2location_merge_audit.json';mj=json.loads(auditpath.read_text());ref=refdir/'GSE273952_cell2location_lightweight.h5ad'
with h5py.File(ref,'r') as f:
 names=strs(f['obs'][f['obs'].attrs['_index']][:]);assert len(set(names))==n and set(names)==set(rawids);idx={k:i for i,k in enumerate(names)};order=np.array([idx[k] for k in rawids]);major=f['obsm/cell2location_q05_major_abundance'][:][order];fine=f['obsm/cell2location_q05_abundance'][:][order]
majors=mj['major_names'];factors=mj['factor_names'];assert len(majors)==6 and major.shape==(n,6) and fine.shape==(n,40);assert np.isfinite(major).all() and (major>=0).all();resum=np.column_stack([fine[:,[k for k,f in enumerate(factors) if f.startswith(m+'_')]].sum(1) for m in majors]);err=float(np.max(abs(resum-major)));assert np.allclose(resum,major,atol=1e-5);fractions=major/major.sum(1)[:,None]
# Hexagonal Visium array coordinates: six immediate neighbours, no cross-section edges.
index={(int(codes[i]),int(rr[i]),int(cc[i])):i for i in range(n)};assert len(index)==n;offsets=[(-1,-1),(-1,1),(0,-2),(0,2),(1,-1),(1,1)];a=[];b=[]
for i in range(n):
 if not keep[i]:continue
 for dr,dc in offsets:
  j=index.get((int(codes[i]),int(rr[i]+dr),int(cc[i]+dc)))
  if j is not None and keep[j]:a.append(i);b.append(j)
g=csr_matrix((np.ones(len(a),np.uint8),(a,b)),shape=(n,n));assert (g-g.T).nnz==0 and g.diagonal().sum()==0;deg=np.asarray(g.sum(1)).ravel();assert deg.max()<=6;assert np.all(codes[np.array(a)]==codes[np.array(b)]);dist=np.sqrt((x[np.array(a)]-x[np.array(b)])**2+(y[np.array(a)]-y[np.array(b)])**2);assert np.isfinite(dist).all() and (dist>0).all()
for name,val in [('TLS_primary',tls),('TLS_secondary',tls2),('tumor_zone',tumor)]:
 levels=sorted(set(val));enc=np.array([levels.index(v) for v in val],np.int16)
 if name=='TLS_primary':tls_codes=enc;tls_levels=levels
 elif name=='TLS_secondary':tls2_codes=enc;tls2_levels=levels
 else:tumor_codes=enc;tumor_levels=levels
np.savez_compressed(O/'private_spatial_graph_and_partial_composition.npz',graph_indptr=g.indptr,graph_indices=g.indices,library_codes=codes,array_row=rr,array_col=cc,pixel_x=x,pixel_y=y,folded_excluded=fold,degree=deg,TLS_primary=tls_codes,TLS_secondary=tls2_codes,tumor_zone=tumor_codes,legacy6_major_fractions=fractions)
(O/'private_region_and_partial_reference_schema.json').write_text(json.dumps(dict(TLS_primary=tls_levels,TLS_secondary=tls2_levels,tumor_zone=tumor_levels,legacy_majors=majors,missing_reference_majors=['Endothelial','Fibroblast','other_12_major_unrepresented'],not_full_TME_composition=True,not_exact_original172_identity=True),indent=2)+'\n')
counts=[]
for name,val in [('TLS_primary',tls),('TLS_secondary',tls2),('tumor_zone',tumor)]:
 for level in sorted(set(val)):
  mask=val==level;counts.append(dict(annotation=name,region=level or 'UNKNOWN_EMPTY',spots=int(mask.sum()),nonfolded_spots=int((mask&keep).sum()),libraries_with_region=int(len(np.unique(codes[mask])))))
pd.DataFrame(counts).to_csv(O/'public_author_region_coverage_aggregate.tsv',sep='\t',index=False)
lib=[]
for c in np.unique(codes):
 ix=np.flatnonzero((codes==c)&keep);sub=g[ix][:,ix];nc,lab=connected_components(sub,directed=False);lib.append(dict(spots=len(ix),edges=int(sub.nnz//2),isolated_spots=int((deg[ix]==0).sum()),components=nc,largest_component_spots=int(np.bincount(lab).max()) if len(ix) else 0))
# No per-library rows or identifiers in public payload.
ld=pd.DataFrame(lib);summary={col:dict(min=int(ld[col].min()),median=float(ld[col].median()),max=int(ld[col].max())) for col in ld}
out=dict(status='PASS_ALL71398_ARRAY_GEOMETRY_AUTHOR_REGIONS_AND_LEGACY6_MAJOR_ALIGNMENT_NOT_FULL_COMPOSITION_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),raw_input_sha256=iv['input_file_sha256'],input_audit_sha256=sha(A/'inventory.json'),historical_composition_file_sha256=sha(ref),historical_composition_audit_sha256=sha(auditpath),spots=n,libraries=19,folded_region_excluded_spots=int(fold.sum()),nonfolded_spots=int(keep.sum()),undirected_immediate_neighbour_edges=int(g.nnz//2),max_degree=int(deg.max()),isolated_nonfolded_spots=int(((deg==0)&keep).sum()),pixel_neighbour_distance_quantiles=[float(t) for t in np.quantile(dist,[.01,.5,.99])],library_geometry_summary=summary,legacy6_major_names=majors,legacy40_to6_max_sum_error=err,legacy_composition_full_TME=False,legacy_composition_missing_endothelium_and_fibroblast=True,patient_bridge_reviewed=False,clinical_outcomes_used=False,biological_spatial_validation=False,public_outputs={'public_author_region_coverage_aggregate.tsv':sha(O/'public_author_region_coverage_aggregate.tsv')},private_input_sha256=sha(O/'private_spatial_graph_and_partial_composition.npz'));(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
