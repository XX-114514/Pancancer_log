from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,h5py
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'GSE273952_spatial_neighbourhood_preflight_20261016_v1';j=json.loads((D/'inventory.json').read_text());assert j['status'].startswith('PASS_ALL71398');z=np.load(D/'private_spatial_graph_and_partial_composition.npz',allow_pickle=False);n=len(z['degree']);g=csr_matrix((np.ones(len(z['graph_indices'])),z['graph_indices'],z['graph_indptr']),shape=(n,n));a,b=g.nonzero();actual=set(zip(a.tolist(),b.tolist()));expected=set()
# Independent Euclidean hex embedding and tree radius, not producer offset lookup.
for lib in np.unique(z['library_codes']):
 ix=np.flatnonzero((z['library_codes']==lib)&~z['folded_excluded']);xy=np.column_stack([z['array_col'][ix]/2,z['array_row'][ix]*np.sqrt(3)/2]);tree=cKDTree(xy)
 for i,k in tree.query_pairs(1.001):expected.add((int(ix[i]),int(ix[k])));expected.add((int(ix[k]),int(ix[i])))
assert actual==expected;assert np.array_equal(np.asarray(g.sum(1)).ravel(),z['degree']);assert np.allclose(z['legacy6_major_fractions'].sum(1),1);assert np.isfinite(z['legacy6_major_fractions']).all() and (z['legacy6_major_fractions']>=0).all();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();O=V/'GSE273952_spatial_geometry_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);out=dict(status='PASS_ALL71398_HEX_GRAPH_INDEPENDENT_KDTREE_EXACT_EDGE_RECONSTRUCTION_AND_PARTIAL_FRACTIONS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),input_inventory_sha256=sha(D/'inventory.json'),private_geometry_sha256=sha(D/'private_spatial_graph_and_partial_composition.npz'),spots=n,undirected_edges=len(actual)//2,old_fractions_complete_TME=False,clinical_outcomes_used=False,biological_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
