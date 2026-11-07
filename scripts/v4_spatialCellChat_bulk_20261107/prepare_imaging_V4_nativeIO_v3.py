"""Native panel transfer diagnostics and conservative broad spatial-CellChat inputs."""
import sys,json,hashlib,tarfile,gzip,io,warnings
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,h5py,joblib
from scipy.io import mmwrite
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_fscore_support
r=PROJECT_ROOT/'runs/20261107_v4_spatialCellChat_bulk_v1';job=json.loads((r/'config/imaging_jobs_PRIVATE.json').read_text())[int(sys.argv[1])];out=r/'imaging'/job['unit_code'];out.mkdir(exist_ok=False)
refpath=PROJECT_ROOT/'runs/20261107_v4_Visium_author_validation_v2/reference/reference_and_heldout_PRIVATE.h5ad';ref=ad.read_h5ad(refpath)
if job['format']=='h5ad':
 a=ad.read_h5ad(job['path']);X=sp.csr_matrix(a.layers['counts']);genes=pd.Index(a.var_names);obs=a.obs.copy();xy=np.asarray(a.obsm['spatial'],float);fov=np.array(['one_tissue']*len(obs))
elif job['format']=='10x_h5':
 with h5py.File(job['path'],'r') as h:
  m=h['matrix'];decode=lambda x:[v.decode() if isinstance(v,bytes) else str(v) for v in x];allgenes=pd.Index(decode(m['features/name'][:]));typ=np.array(decode(m['features/feature_type'][:]));full=sp.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]));keep=typ=='Gene Expression';X=full[keep,:].T.tocsr();genes=allgenes[keep];bar=pd.Index(decode(m['barcodes'][:]))
 obs=pd.read_csv(job['coordinates_path']).set_index('cell_id');assert obs.index.is_unique and bar.is_unique and set(bar)<=set(obs.index);obs=obs.loc[bar];xy=obs[['x_centroid','y_centroid']].to_numpy(float);fov=np.array(['one_tissue']*len(obs))
else:
 ex=None;obs=None
 with tarfile.open(job['path'],'r|gz') as tf:
  for member in tf:
   name=member.name.lower()
   if name.endswith('exprmat_file.csv'):ex=pd.read_csv(io.BytesIO(tf.extractfile(member).read()))
   elif name.endswith('metadata_file.csv'):obs=pd.read_csv(io.BytesIO(tf.extractfile(member).read()))
   if ex is not None and obs is not None:break
 assert ex is not None and obs is not None,'Both expression and metadata required from original archive'
 cell='cell_ID' if 'cell_ID' in ex else 'cell_id';keys=['fov',cell];assert set(keys)<=set(obs) and set(keys)<=set(ex);ex=ex.loc[ex[cell].ne(0)].copy();obs=obs.loc[obs[cell].ne(0)].copy()
 ix=pd.MultiIndex.from_frame(ex[keys]);mi=pd.MultiIndex.from_frame(obs[keys]);assert ix.is_unique and mi.is_unique;obs.index=mi;obs=obs.reindex(ix);assert obs[keys].notna().all().all();good=ex[cell].ne(0).to_numpy();drop=set(keys+['Run_Tissue_name']);genes=pd.Index([c for c in ex if c not in drop and not any(z in c.lower() for z in ['negprb','negative','systemcontrol','falsecode'])]);X=sp.csr_matrix(ex.loc[good,genes].to_numpy(np.float32));obs=obs.loc[good].copy();obs.index=pd.Index(['c%08d'%i for i in range(len(obs))]);xy=obs[['CenterX_global_px','CenterY_global_px']].to_numpy(float)*job['pixel_um'];fov=obs.fov.astype(str).to_numpy()
assert genes.is_unique and np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.round(X.data));assert len(obs)==len(xy)==X.shape[0];n_original=len(obs);finite=np.isfinite(xy).all(1);native_total=np.asarray(X.sum(1)).ravel();common=genes.intersection(ref.var_names,sort=False);assert len(common)>=50
Q=sp.csr_matrix(X[:,genes.get_indexer(common)].astype(np.float64));R=sp.csr_matrix(ref[:,common].X.astype(np.float64));qgenes=np.asarray((Q>0).sum(1)).ravel();qtotal=np.asarray(Q.sum(1)).ravel();rtotal=np.asarray(R.sum(1)).ravel();qc=finite&(native_total>=50)&(qgenes>=20)&(qtotal>0)
def norm(z,total):
 z=sp.diags(1e4/np.maximum(total,1))@z;z=z.tocsr();z.data=np.log1p(z.data);return z
qn=norm(Q,qtotal);rn=norm(R,rtotal);train=ref.obs['split'].eq('reference').to_numpy();hold=~train;reports=[];pred={};prob={};goodclass={}
def fit_or_reuse(model,mode,y):
 # Reuse only this run's completed exact same panel/reference/parameters.
 for previous in sorted((r/'imaging').glob('I*')):
  if previous==out or not (previous/'prepare_complete.json').exists():continue
  gf=previous/'genes.tsv';mf=previous/(mode+'_classifier_PRIVATE.joblib')
  if not gf.exists() or not mf.exists():continue
  previousgenes=pd.Index(pd.read_csv(gf,header=None)[0]);pc=previousgenes.intersection(ref.var_names,sort=False)
  if pc.equals(common):
   cached=joblib.load(mf);assert np.array_equal(cached.classes_,np.unique(y[train]));print(json.dumps({'stage':'reuse_same_panel_reference_classifier','resolution':mode,'same_common_gene_order':True}),flush=True);return cached
 model.fit(rn[train],y[train]);return model
for name,column in [('fine','node'),('broad','major')]:
 y=ref.obs[column].astype(str).to_numpy();model=LogisticRegression(C=1,class_weight='balanced',solver='lbfgs',max_iter=1000,random_state=2026110703);model=fit_or_reuse(model,name,y);hp=model.predict(rn[hold]);pr,rc,f1,support=precision_recall_fscore_support(y[hold],hp,labels=model.classes_,zero_division=0);qual={}
 for lab,pv,rv,fv,n in zip(model.classes_,pr,rc,f1,support):
  passed=bool(pv>=.7 and fv>=.7 and n>=20);qual[lab]=passed;reports.append({'resolution':name,'label':lab,'precision':float(pv),'recall':float(rv),'F1':float(fv),'heldout_support':int(n),'reference_replay_gate':passed,'not_independent_spatial_truth':True})
 pp=[];vv=[]
 for start in range(0,len(obs),20000):
  z=model.predict_proba(qn[start:start+20000]);rank=np.argsort(z,axis=1);pp.extend(model.classes_[rank[:,-1]]);vv.extend(np.column_stack([np.max(z,axis=1),np.max(z,axis=1)-np.partition(z,-2,axis=1)[:,-2]]))
 pred[name]=np.array(pp);prob[name]=np.asarray(vv);goodclass[name]=qual
 joblib.dump(model,out/(name+'_classifier_PRIVATE.joblib'))
# Explicit marker support for CM-relevant broad identities; absent markers remain unmeasured.
lookup={g:i for i,g in enumerate(genes)}
def det(g):return np.asarray(X[:,lookup[g]].toarray()).ravel()>0 if g in lookup else np.zeros(len(obs),bool)
def nmarkers(gg):return np.column_stack([det(g) for g in gg]).sum(1)
identity_support={'Bcell':nmarkers(['MS4A1','CD79A','CD79B','CD37','CD74'])>=2,'Endothelial':nmarkers(['PECAM1','VWF','KDR','EMCN','CDH5','ESAM','PLVAP'])>=2,'NK':(nmarkers(['NKG7','KLRD1','GNLY','NCR1'])>=2)&(nmarkers(['CD3D','CD3E','TRAC'])<2),'Macrophage':nmarkers(['LST1','CSF1R','CD68','CTSS','C1QA','C1QB','FCER1G','TYROBP'])>=2,'Monocyte':nmarkers(['LST1','FCN1','LYZ','S100A8','S100A9','CTSS'])>=2}
broad_ok=qc&(prob['broad'][:,0]>=.7)&(prob['broad'][:,1]>=.2)&np.array([goodclass['broad'].get(x,False) for x in pred['broad']])
for label,evidence in identity_support.items():broad_ok&=(pred['broad']!=label)|evidence
fine_ok=qc&(prob['fine'][:,0]>=.7)&(prob['fine'][:,1]>=.2)&np.array([goodclass['fine'].get(x,False) for x in pred['fine']]);fine_major=ref.obs.groupby('node',observed=True).major.first().astype(str).to_dict();fine_ok&=np.array([fine_major[x] for x in pred['fine']])==pred['broad'];fine_ok&=broad_ok
# All naming-marker tokens in the reference axis must be measured and positive.
for node in np.unique(pred['fine']):
 marker=list(dict.fromkeys(t for t in node.split('_') if t in set(ref.var_names)));has=all(t in lookup for t in marker);m=pred['fine']==node
 if not has or not marker:fine_ok[m]=False
 else:fine_ok[m]&=np.column_stack([det(t) for t in marker])[m].all(1)
rep=pd.DataFrame(reports);rep.to_csv(out/'reference_panel_heldout_diagnostics.tsv',sep='\t',index=False);pd.DataFrame({'gene':genes}).to_csv(out/'native_panel_genes.tsv',sep='\t',index=False)
counts=pd.DataFrame({'cell_id':['c%08d'%i for i in range(len(obs))],'broad':pred['broad'],'fine':pred['fine'],'QC':qc,'broad_accepted':broad_ok,'fine_accepted':fine_ok,'probability_broad':prob['broad'][:,0],'probability_fine':prob['fine'][:,0],'source_FOV':fov,'x_um':xy[:,0],'y_um':xy[:,1]}).set_index('cell_id');counts.to_csv(out/'cell_transfer_PRIVATE.tsv.gz',sep='\t')
summary=counts.groupby(['broad','broad_accepted'],observed=True).size().rename('n_cells').reset_index();summary.to_csv(out/'broad_transfer_counts.tsv',sep='\t',index=False)
# Canonical native full-panel normalized expression: only accepted broad cells, no imputed genes.
selected=np.flatnonzero(broad_ok);x=X[selected].astype(np.float64);x=norm(x,native_total[selected]);mmwrite(out/'normalized_expression_PRIVATE.mtx',x.T);pd.Series(genes).to_csv(out/'genes.tsv',index=False,header=False,sep='\t');meta=counts.iloc[selected][['broad','source_FOV','x_um','y_um']].copy();meta.to_csv(out/'metadata_PRIVATE.tsv',sep='\t')
members=pd.read_csv(PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1/audit/network_members.tsv',sep='\t');cmrows=[]
for cm in ['CM%02d'%i for i in range(1,15)]:
 nodes=members.loc[members.cm.eq(cm),'subCluster'].tolist();measured={node:all(g in lookup for g in node.split('_') if g in set(ref.var_names)) for node in nodes};gate=bool(nodes and all(measured.values()) and all(goodclass['fine'].get(node,False) for node in nodes));accepted={node:int(((pred['fine']==node)&fine_ok).sum()) for node in nodes};eligible=gate and min(accepted.values())>=10 if accepted else False
 cmrows.append({'cm':cm,'n_original_members':len(nodes),'all_naming_markers_measured':bool(nodes and all(measured.values())),'all_original_identity_replay_gates_pass':gate,'all_members_min10_marker_supported_cells':eligible,'exact_CM_interpretation_eligible':eligible,'minimum_member_accepted_cells':min(accepted.values()) if accepted else 0,'status':'EXACT_IDENTITY_OPPORTUNITY_PENDING_SPATIAL_RECURRENCE' if eligible else 'NOT_EXACT_CM_MEASUREMENT_BROAD_CONTEXT_ONLY'})
pd.DataFrame(cmrows).to_csv(out/'all14_exact_CM_measurement_eligibility.tsv',sep='\t',index=False)
receipt={'status':'NATIVE_IMAGING_BROAD_INPUT_PREPARED_PENDING_CELLCHAT_AND_AUDIT','UTC':datetime.now(timezone.utc).isoformat(),'dataset':job['dataset'],'platform':job['platform'],'n_original_cells':n_original,'n_measured_genes':len(genes),'n_reference_common_genes':len(common),'n_QC_cells':int(qc.sum()),'n_broad_accepted_cells':len(selected),'n_broad_groups':int(meta.broad.nunique()),'n_source_FOV':int(meta.source_FOV.nunique()),'CM03_exact_eligible':bool(next(x for x in cmrows if x['cm']=='CM03')['exact_CM_interpretation_eligible']),'original_annotation_changed':False,'reference_replay_not_independent_truth':True,'coordinate_unit':'micrometre','cross_FOV_neighbors_forbidden':True};(out/'prepare_complete.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
