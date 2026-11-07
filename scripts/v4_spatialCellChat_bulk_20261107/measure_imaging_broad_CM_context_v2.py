"""All14 broad-parent context opportunities; within-FOV fixed-count spatial null."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,numpy as np,pandas as pd
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
r=PROJECT_ROOT/'runs/20261107_v4_spatialCellChat_bulk_v1';j=json.loads((r/'config/imaging_jobs_PRIVATE.json').read_text())[int(sys.argv[1])];out=r/'imaging'/j['unit_code'];meta=pd.read_csv(out/'metadata_PRIVATE.tsv',sep='\t',index_col=0);members=pd.read_csv(PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1/audit/network_members.tsv',sep='\t');ref=pd.read_csv(PROJECT_ROOT/'runs/20261107_v4_Visium_author_validation_v2/reference/node_reference_coverage.tsv',sep='\t').set_index('node');labels=sorted(meta.broad.unique());K=len(labels);codes=pd.Categorical(meta.broad,categories=labels).codes.astype(np.int64);groups=[];src=[];tgt=[]
for _,ii in meta.groupby('source_FOV',sort=True).indices.items():
 ii=np.asarray(ii);groups.append(ii)
 if len(ii)<2:continue
 xy=meta.iloc[ii][['x_um','y_um']].to_numpy();k=min(7,len(ii));d,nn=cKDTree(xy).query(xy,k=k,workers=1);aa=np.repeat(np.arange(len(ii)),k-1);bb=nn[:,1:].ravel();ok=(d[:,1:].ravel()<=100)&(aa!=bb);src.extend(ii[aa[ok]]);tgt.extend(ii[bb[ok]])
src=np.asarray(src,int);tgt=np.asarray(tgt,int);N=len(src);masks=[];records=[];counts=meta.broad.value_counts()
for cm in ['CM%02d'%i for i in range(1,15)]:
 nodes=members.loc[members.cm.eq(cm),'subCluster'];parents=sorted(set(ref.loc[nodes,'major'])) if len(nodes) else [];present=all(counts.get(x,0)>=10 for x in parents);mask=np.zeros((K,K),bool)
 for a in parents:
  for b in parents:
   if a!=b and a in labels and b in labels:mask[labels.index(a),labels.index(b)]=True
 eligible=len(parents)>=2 and present and mask.any() and N>0;masks.append(mask.ravel().astype(float)/mask.sum() if eligible else np.zeros(K*K));records.append({'cm':cm,'broad_parents':';'.join(parents),'n_broad_parents':len(parents),'n_spatial_edges':N,'eligible':bool(eligible),'status':'broad_context_only' if eligible else 'empty_single_parent_or_insufficient_supported_cells','exact_subtype_CM_replication':False})
M=np.vstack(masks);obs=M@np.bincount(codes[src]*K+codes[tgt],minlength=K*K)/max(N,1);rng=np.random.default_rng(2026110703);null=np.empty((999,14))
for b in range(999):
 y=codes.copy()
 for ii in groups:y[ii]=rng.permutation(y[ii])
 null[b]=M@np.bincount(y[src]*K+y[tgt],minlength=K*K)/max(N,1)
for i,z in enumerate(records):
 sd=float(null[:,i].std(ddof=1));z.update(observed_mean_directed_parent_pair_edge_fraction=float(obs[i]) if z['eligible'] else np.nan,null_mean=float(null[:,i].mean()) if z['eligible'] else np.nan,null_sd=sd if z['eligible'] else np.nan,p=float((1+np.sum(null[:,i]>=obs[i]))/1000) if z['eligible'] and sd>0 else np.nan,z_enrichment=float((obs[i]-null[:,i].mean())/sd) if z['eligible'] and sd>0 else np.nan)
t=pd.DataFrame(records);t['BH_all14']=multipletests(t.p.fillna(1),method='fdr_bh')[1];t.loc[t.p.isna(),'BH_all14']=np.nan;t.to_csv(out/'all14_broad_CM_context.tsv',sep='\t',index=False);np.savez_compressed(out/'broad_context_null_PRIVATE.npz',observed=obs,null=null,seed=np.array([2026110703]));proof={'status':'ALL14_BROAD_PARENT_SPATIAL_CONTEXT_MEASURED_PENDING_INDEPENDENT_REVIEW','UTC':datetime.now(timezone.utc).isoformat(),'n_permutations':999,'neighbors':6,'maximum_edge_distance_um':100,'n_directed_edges':N,'n_assessable_CM_contexts':int(t.p.notna().sum()),'fixed_FOV_label_counts':True,'segmented_cell_contact_or_exact_CM_or_patient_claim':False};(out/'context_complete.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps(proof),flush=True)
