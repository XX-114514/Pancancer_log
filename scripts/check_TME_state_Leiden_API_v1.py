"""Real Scanpy graph, Leiden and frozen projection on an orthogonal toy model."""
import argparse,datetime,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,scanpy as sc,anndata
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import adjusted_rand_score

def main(R):
 rng=np.random.default_rng(42);latent=np.repeat([0,1,2],120);centres=np.eye(3)*8;y=centres[latent]@rng.normal(size=(3,12))+rng.normal(0,.25,(360,12));source=np.tile(np.arange(3),120)
 train=np.r_[np.arange(0,80),np.arange(120,200),np.arange(240,320)];test=np.setdiff1d(np.arange(360),train)
 scale=StandardScaler().fit(y[train]);pca=PCA(n_components=8,svd_solver='full').fit(scale.transform(y[train]));pc=pca.transform(scale.transform(y[train]));a=anndata.AnnData(scale.transform(y[train]).astype(np.float32));a.obsm['X_state_PCA']=pc
 sc.pp.neighbors(a,n_neighbors=15,use_rep='X_state_PCA',random_state=42);sc.tl.leiden(a,resolution=.6,random_state=42,flavor='igraph',directed=False,n_iterations=2,key_added='state_pattern');label=a.obs.state_pattern.astype(str).to_numpy();ari=adjusted_rand_score(latent[train],label);assert ari>.8
 model=KNeighborsClassifier(n_neighbors=15,weights='distance').fit(pc,label);prediction=model.predict(pca.transform(scale.transform(y[test])));testari=adjusted_rand_score(latent[test],prediction);assert testari>.8
 assert np.unique(source[train]).size==3 and np.unique(source[test]).size==3
 p=R/'outputs/accepted_state_Leiden_API_synthetic_v1.json';p.write_text(json.dumps(dict(schema='TME_state_Leiden_projection_synthetic_v1',all_tests_passed=True,n_tests=3,training_ARI=ari,frozen_projection_ARI=testari,source_balanced_synthetic_design=True,official_scanpy_version=sc.__version__,not_human_state_validation=True),indent=2))
 (R/'provenance/terminal_v1.json').write_text(json.dumps(dict(status='OFFICIAL_SCANPY_STATE_GRAPH_LEIDEN_AND_PROJECTION_SYNTHETIC_PASS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  (R/'provenance/failure_v1.json').write_text(json.dumps(dict(traceback=traceback.format_exc()),indent=2));raise
