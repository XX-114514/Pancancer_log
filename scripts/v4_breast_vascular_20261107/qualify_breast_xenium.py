"""Native-scale and author-label qualification plus descriptive vascular profiles.
No patient-level inference, exact CM03 recovery, HEV assignment or causal claim.
"""
from pathlib import Path
import json,sys,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
PANELS={'EC':['PECAM1','VWF','CDH5','KDR','EMCN'],'B':['MS4A1','CD79A','CD79B','CD37'],'T':['CD3E','TRBC1','TRBC2'],'NK_anchor':['NCR1','KLRD1'],'NK_effector':['GNLY','FCGR3A','PRF1'],'lymphatic':['PROX1','PDPN','LYVE1','FLT4'],'mural':['RGS5','PDGFRB','ACTA2','TAGLN','MYH11'],'fibroblast':['COL1A1','COL1A2','DCN','PDGFRA']}
AXES=['ACKR1','CCL19','CCL21','CCR7','CXCL12','CXCR4','IGHD','XCL1']
def run(run_root,unit,input_dir):
 r=Path(run_root);inp=Path(input_dir);out=r/'spatial'/unit;out.mkdir(parents=True,exist_ok=False)
 receipt=json.loads((inp/'complete_PRIVATE.json').read_text());a=ad.read_h5ad(inp/receipt['members'][0]['local_name'],backed='r')
 assert a.obs_names.is_unique and a.var_names.is_unique
 genes=a.var_names.astype(str)
 if 'gene_symbol' in a.var:genes=a.var['gene_symbol'].astype(str)
 assert genes.is_unique
 needed=sorted(set(AXES+sum(PANELS.values(),[])));ix={g:i for i,g in enumerate(genes)};available=[g for g in needed if g in ix];order=sorted(ix[g] for g in available);names=[str(genes[i]) for i in order]
 x=a[:,order].X;x=x.toarray() if sparse.issparse(x) else np.asarray(x);assert np.isfinite(x).all() and (x>=0).all()
 counts=np.expm1(x.astype(float));err=np.abs(counts-np.rint(counts));assert np.all(err<=1e-3+1e-5*np.abs(counts)), 'X not consistent with author native log1p counts'
 counts=np.rint(counts);cols={g:counts[:,j] for j,g in enumerate(names)}
 coords=a.obs[['x_centroid','y_centroid']].to_numpy(float);assert coords.shape==(a.n_obs,2) and np.isfinite(coords).all()
 if 'xenium_um_spatial' in a.obsm:assert np.allclose(np.asarray(a.obsm['xenium_um_spatial']),coords,atol=1e-3)
 maps=json.loads((r/'sources/mapping/complete_PRIVATE.json').read_text());sample=Path(receipt['members'][0]['member']).stem;labels={};mapping_qual={}
 for lev in ['major','minor','subset']:
  m=[z for z in maps['members'] if f'/{lev}_level/' in z['member'] and Path(z['member']).name.startswith(sample+'_')];assert len(m)==1
  t=pd.read_csv(r/'sources/mapping'/m[0]['local_name'],index_col=0,dtype={'pred_cell_type':str},float_precision='round_trip');assert t.index.is_unique
  assert set(a.obs_names)==set(t.index), 'Author mapping native barcode coverage differs'
  t=t.reindex(a.obs_names);assert t['pred_cell_type'].notna().all();labels[lev]=t['pred_cell_type'].astype(str)
  q=t['pred_cell_type_prob'].to_numpy(float);mapping_qual[lev]={'n_labels':int(t['pred_cell_type'].nunique()),'n_rows':len(t),'score_min':float(np.min(q)),'score_max':float(np.max(q)),'score_median':float(np.median(q)),'score_semantics':'cluster-to-location Tangram mass, not calibrated class probability; no 0.7 confidence threshold'}
 def positive(panel,min_genes):
  measured=[g for g in PANELS[panel] if g in cols]
  return np.sum([cols[g]>=1 for g in measured],axis=0)>=min_genes if len(measured)>=min_genes else np.zeros(a.n_obs,dtype=bool)
 major=labels['major'].str.casefold();minor=labels['minor'].str.casefold();sub=labels['subset'].str.casefold()
 endo=major.str.contains('endo').to_numpy()&positive('EC',2);lymph=endo&positive('lymphatic',2);blood=endo&~lymph
 b=major.str.fullmatch(r'b[\s_-]*cells?').to_numpy()&positive('B',2)
 t=major.str.contains(r'^t[\s_-]*cell',regex=True).to_numpy()&positive('T',2)
 nklabel=(minor.str.contains(r'\bnk\b',regex=True)|sub.str.contains(r'\bnk\b',regex=True)).to_numpy()
 nk=nklabel&positive('NK_anchor',1)&positive('NK_effector',1)&~positive('T',2)
 conflicts=(b.astype(int)+t.astype(int)+nk.astype(int)+blood.astype(int)+lymph.astype(int))>1
 assert not conflicts.any(),'Lineage masks overlap'
 identity={'B':b,'T':t,'NK':nk,'blood_EC':blood,'lymphatic_EC':lymph,'ACKR1_positive_blood_EC':blood&(cols.get('ACKR1',np.zeros(a.n_obs))>=1),'ACKR1_undetected_blood_EC':blood&(cols.get('ACKR1',np.zeros(a.n_obs))<1)}
 rows=[]
 for g in needed:rows.append({'gene':g,'measured':g in cols,'n_detected':int((cols[g]>=1).sum()) if g in cols else None})
 pd.DataFrame(rows).to_csv(out/'measurement_coverage.tsv',sep='\t',index=False)
 labrows=[]
 for level,v in labels.items():
  for k,n in v.value_counts().items():labrows.append({'level':level,'author_label':k,'n_cells':int(n)})
 pd.DataFrame(labrows).to_csv(out/'author_label_counts.tsv',sep='\t',index=False)
 expr=[]
 for label in sorted(major.unique()):
  mask=major.eq(label).to_numpy()
  for g in AXES:
   expr.append({'author_major':label,'gene':g,'measured':g in cols,'n_cells':int(mask.sum()),'n_detected':int((cols[g][mask]>=1).sum()) if g in cols else None,'mean_native_UMI':float(cols[g][mask].mean()) if g in cols else None})
 pd.DataFrame(expr).to_csv(out/'chemokine_source_descriptive.tsv',sep='\t',index=False)
 profiles=[];edges=np.array([0,25,50,100,200,400,np.inf]);nonendo=~endo
 for anchor in ['blood_EC','ACKR1_positive_blood_EC','ACKR1_undetected_blood_EC']:
  anchor_mask=identity[anchor]
  if not anchor_mask.any():continue
  d=cKDTree(coords[anchor_mask]).query(coords,k=1)[0];bins=np.searchsorted(edges,d,side='right')-1
  for j in range(len(edges)-1):
   den=nonendo&(bins==j)
   for typ in ['B','T','NK']:
    num=identity[typ]&den;profiles.append({'anchor':anchor,'lower_um':edges[j],'upper_um':edges[j+1],'cell_type':typ,'n_nonendothelial_cells':int(den.sum()),'n_lineage_cells':int(num.sum()),'fraction':float(num.sum()/den.sum()) if den.any() else None,'distance_definition':'native micron centroid to nearest qualifying endothelial centroid; not vessel wall'})
 pd.DataFrame(profiles).to_csv(out/'descriptive_distance_profiles.tsv',sep='\t',index=False)
 result={'UTC':datetime.now(timezone.utc).isoformat(),'unit_alias':unit,'status':'NATIVE_SCHEMA_AND_DESCRIPTIVE_PROFILES_COMPLETE_NOT_FORMAL_CM_OR_COMMUNICATION_VALIDATION','n_cells':a.n_obs,'n_genes':a.n_vars,'n_marker_genes_reviewed':len(available),'X_scale':'native log1p counts independently expm1/integer checked on reviewed genes','coordinates':'native x_centroid/y_centroid microns; image-space obsm spatial not used','mapping':mapping_qual,'identities':{k:int(v.sum()) for k,v in identity.items()},'missing_exact_CM03_genes':[g for g in ['IGHD','XCL1'] if g not in cols],'no_automatic_HEV_or_PVL_chemokine_assignment':True,'no_patient_level_inference':True,'no_state_prerequisite':True}
 (out/'complete.json').write_text(json.dumps(result,indent=2)+'\n');a.file.close();print(json.dumps(result),flush=True);return result
if __name__=='__main__':run(*sys.argv[1:])
