"""RNA marker-defined spatial pilot, independent of location-derived vendor labels.
Single donor. Candidate identities only; no LR, survival, or patient-level tests.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd,h5py
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
ROOT=Path(sys.argv[1]);src=ROOT/'inputs';out=ROOT/'spatial_pilot_v1';out.mkdir(exist_ok=False)
rec=json.loads((src/'complete_PRIVATE.json').read_text());files={}
for row in rec['members']:
 p=src/row['local_name'];assert p.stat().st_size==row['bytes'];assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'];files[row['member']]=p
with h5py.File(files['cell_feature_matrix.h5'],'r') as f:
 m=f['matrix'];x=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr();genes=[v.decode() for v in m['features']['name'][:]];ids=[v.decode() for v in m['barcodes'][:]];ft=[v.decode() for v in m['features']['feature_type'][:]]
c=pd.read_parquet(files['cells.parquet']);assert c.cell_id.is_unique;c=c.set_index('cell_id').reindex(ids);assert c[['x_centroid','y_centroid']].notna().all().all()
rna=np.array([v=='Gene Expression' for v in ft]);tot=np.asarray(x[:,rna].sum(1)).ravel();ng=x[:,rna].getnnz(1)
panels={'B':['CD19','MS4A1','CD79A'],'T':['CD3D','CD3E','TRAC'],'NK':['NKG7','GNLY','KLRD1','NCAM1','PRF1'],'Blood_EC':['PECAM1','VWF','EMCN','ACKR1'],'Lymphatic':['PROX1','LYVE1','PDPN'],'Mural':['PDGFRB','ACTA2','MYH11']}
needed=sorted(set(sum(panels.values(),[])));gix={g:genes.index(g) for g in needed if g in genes};a={g:x[:,j].toarray().ravel() for g,j in gix.items()}
def pos(g):return a.get(g,np.zeros(len(ids)))>=1
def count(p):return np.sum([pos(g) for g in panels[p]],axis=0)
b=count('B')>=2;t=count('T')>=2;ec=(count('Blood_EC')>=2)&(pos('PECAM1')|pos('VWF'));ly=count('Lymphatic')>=2
nk=pos('KLRD1')&(pos('NKG7')|pos('GNLY'))&(count('T')==0)&~b&~ec
m=(count('Mural')>=2)&~ec&~b&~t
masks={'B':b,'T':t,'NK':nk,'Blood_EC':ec&~ly,'Lymphatic':ec&ly,'Mural':m}
conf=np.sum(list(masks.values()),axis=0)>1;qc=(tot>=50)&(ng>=20)&~conf
identity=np.full(len(ids),'Other',dtype=object)
for k,v in masks.items():identity[v&~conf]=k
c['identity_candidate']=identity;c['QC_keep']=qc;c['rna_total']=tot;c['rna_genes']=ng;c.to_csv(out/'identity_candidates_PRIVATE.tsv.gz',sep='\t')
review=[]
for k in masks:
 sel=qc&(identity==k)
 for p,gg in panels.items():
  for g in gg:review.append({'identity_candidate':k,'n':int(sel.sum()),'gene':g,'marker_identity':p,'available':g in a,'fraction_positive':float(pos(g)[sel].mean()) if sel.any() and g in a else None})
pd.DataFrame(review).to_csv(out/'marker_identity_review.tsv',sep='\t',index=False)
coords=c[['x_centroid','y_centroid']].to_numpy();ecsel=qc&(identity=='Blood_EC');assert ecsel.sum()>20
valid=qc&(identity!='Blood_EC')&(identity!='Lymphatic');dist=np.full(len(ids),np.nan);dist[valid]=cKDTree(coords[ecsel]).query(coords[valid])[0]
rows=[]
for k in ['B','NK','T','Mural','Other']:
 z=dist[qc&(identity==k)];z=z[np.isfinite(z)]
 if not len(z):continue
 row={'identity_candidate':k,'n':len(z),'median_nearest_EC_centroid_um':float(np.median(z))}
 for radius in [25,50,100,200]:row[f'fraction_within_{radius}um']=float(np.mean(z<=radius))
 rows.append(row)
pd.DataFrame(rows).to_csv(out/'distance_descriptive.tsv',sep='\t',index=False)
# Fixed-coordinate null: shuffle immune identity within 250um tiles, preserving
# local immune abundance. Single-section conditional test only, not patient replication.
immune=qc&np.isin(identity,['B','NK','T']);ii=np.flatnonzero(immune);xy=coords[ii];lab=identity[ii];di=dist[ii]
tiles=np.floor(xy/250).astype(int);_,inv=np.unique(tiles,axis=0,return_inverse=True);groups=[np.flatnonzero(inv==j) for j in np.unique(inv)];rng=np.random.default_rng(20261108);nullrows=[]
for k in ['B','NK']:
 if (lab==k).sum()<20:continue
 observed=float(np.mean(di[lab==k]<=50));null=[]
 for repeat in range(999):
  shuffled=lab.copy()
  for ix in groups:shuffled[ix]=rng.permutation(lab[ix])
  null.append(float(np.mean(di[shuffled==k]<=50)))
 null=np.array(null);nullrows.append({'identity_candidate':k,'radius_um':50,'tile_um':250,'n_permutations':999,'observed_fraction':observed,'null_mean':float(null.mean()),'difference':float(observed-null.mean()),'p_one_sided_conditional':float((1+(null>=observed).sum())/1000),'interpretation':'single-section exploratory; fixed local immune counts; no patient-level inference'})
pd.DataFrame(nullrows).to_csv(out/'within_tile_immune_null.tsv',sep='\t',index=False)
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(8,7));sample=np.flatnonzero(qc)[::15];ax.scatter(coords[sample,0],coords[sample,1],s=.1,c='#dddddd',rasterized=True)
for k,color in [('Blood_EC','#d64a3b'),('B','#2676b8'),('NK','#33a36b')]:
 sel=qc&(identity==k);ax.scatter(coords[sel,0],coords[sel,1],s=.4,color=color,label=f'{k}: marker-supported candidate (n={sel.sum():,})',rasterized=True)
ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('Native x (micrometre)');ax.set_ylabel('Native y (micrometre)');ax.legend(markerscale=5,fontsize=7);ax.set_title('ccRCC single-donor spatial pilot; candidate identities');fig.tight_layout();fig.savefig(out/'spatial_candidates.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(12,6));types=list(masks);geneorder=[(g,p) for p,gg in panels.items() for g in gg]
for iy,k in enumerate(types):
 sel=qc&(identity==k)
 for ix,(g,p) in enumerate(geneorder):
  if g in a and sel.any():dots=ax.scatter(ix,iy,s=pos(g)[sel].mean()*100,c=[np.log1p(a[g][sel]).mean()],vmin=0,vmax=2,cmap='viridis')
fig.colorbar(dots,ax=ax,label='Mean log1p RNA UMI',fraction=.025,pad=.02);ax.set_xticks(range(len(geneorder)),[g+' | '+p for g,p in geneorder],rotation=90,fontsize=7);ax.set_yticks(range(len(types)),types);ax.set_title('Marker + expected identity; size = positive fraction, colour = mean log1p UMI');fig.tight_layout();fig.savefig(out/'marker_identity_panel.png',dpi=180);plt.close(fig)
receipt={'cells_input':len(ids),'QC_pass':int(qc.sum()),'candidate_counts':{k:int((qc&(identity==k)).sum()) for k in masks},'feature_types':{str(k):int(v) for k,v in pd.Series(ft).value_counts().items()},'distance':rows,'conditional_null':nullrows,'status':'pilot_complete_annotation_not_frozen','caveats':['one donor','RNA marker candidates require cluster/protein review','no vendor spatial labels used for inference','nearest EC centroid, not vessel boundary','TLS and histological confounding not excluded'],'communications_run':False}
(out/'complete.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
