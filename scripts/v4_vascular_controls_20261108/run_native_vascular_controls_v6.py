"""Independent native-UMI identity/QC view and prospective vascular controls."""
from pathlib import Path
import sys,json,hashlib
from datetime import datetime,timezone
import anndata as ad,numpy as np,pandas as pd
from scipy import sparse
from native_coordinate_contract_v1 import native_micron_coordinates
from vascular_spatial_controls_v2 import neighborhood,regional_ES,edge_expression_null,AXES
PANELS={'B':['MS4A1','CD79A','CD79B','CD37'],'T':['CD3E','CD3G'],'NK_anchor':['NCR1','KLRD1'],'NK_effector':['GNLY','FCGR3A','PRF1'],'EC':['PECAM1','VWF','CDH5','KDR','EMCN'],'Lymphatic':['PROX1','PDPN','LYVE1','FLT4'],'Mural':['RGS5','PDGFRB','ACTA2','TAGLN','MYH11'],'Fibroblast':['COL1A1','COL1A2','DCN','PDGFRA'],'Venular_context':['ACKR1','SELE','SELP','CCL14'],'Arterial_context':['GJA5','EFNB2','SOX17'],'Capillary_context':['CA4','RGCC','EMCN']}
def run(w,unit,native_folder,mapping_folder,expected_member=None,panel_gene_names=None):
 w=Path(w);out=w/'spatial_controls_v1'/unit;out.mkdir(parents=True,exist_ok=False);native_folder=Path(native_folder);mapping_folder=Path(mapping_folder)
 z=json.loads((native_folder/'complete_PRIVATE.json').read_text());assert len(z['members'])==1
 if expected_member is not None:assert z['members'][0]['member']==expected_member, 'Wrong source object for expected donor/section'
 f=native_folder/z['members'][0]['local_name'];assert f.stat().st_size==z['members'][0]['bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==z['members'][0]['sha256']
 a=ad.read_h5ad(f,backed='r');assert a.obs_names.is_unique and a.var_names.is_unique
 coords,coordinate_receipt=native_micron_coordinates(a,author_pixel_size_um=0.2125)
 symbols=pd.Index(a.var['gene_symbol'].astype(str) if 'gene_symbol' in a.var else a.var_names.astype(str));assert symbols.is_unique
 needed=sorted(set(sum(PANELS.values(),[])+['ACKR1']+sum(([l,r] for l,r in AXES),[])));idx=sorted(symbols.get_loc(g) for g in needed if g in symbols);genes=list(symbols[idx]);x=a[:,idx].X;x=x.toarray() if sparse.issparse(x) else np.asarray(x);x=np.expm1(x.astype(float));assert np.isfinite(x).all() and np.all(np.abs(x-np.rint(x))<=1e-3+1e-5*np.abs(x));x=np.rint(x)
 totals=np.zeros(a.n_obs);detected=np.zeros(a.n_obs,int)
 for start in range(0,a.n_obs,10000):
  v=a.X[start:start+10000]
  if sparse.issparse(v):
   v=v.copy().tocsr();v.data=np.expm1(v.data.astype(float));assert np.all(np.abs(v.data-np.rint(v.data))<=1e-3+1e-5*np.abs(v.data));v.data=np.rint(v.data);totals[start:start+len(v.indptr)-1]=np.asarray(v.sum(axis=1)).ravel();detected[start:start+len(v.indptr)-1]=np.asarray((v>0).getnnz(axis=1))
  else:
   v=np.expm1(np.asarray(v,float));assert np.all(np.abs(v-np.rint(v))<=1e-3+1e-5*np.abs(v));v=np.rint(v);totals[start:start+len(v)]=v.sum(axis=1);detected[start:start+len(v)]=(v>0).sum(axis=1)
 mappings=json.loads((mapping_folder/'complete_PRIVATE.json').read_text());sample=Path(z['members'][0]['member']).stem;lab={};axis_hash=hashlib.sha256('\n'.join(a.obs_names.astype(str)).encode()).hexdigest()
 for level in ['major','minor','subset']:
  item=[v for v in mappings['members'] if '/'+level+'_level/' in v['member'] and Path(v['member']).name.startswith(sample+'_')];assert len(item)==1;v=item[0];path=mapping_folder/v['local_name'];assert hashlib.sha256(path.read_bytes()).hexdigest()==v['sha256'];t=pd.read_csv(path,index_col=0,dtype={'pred_cell_type':str});assert t.index.is_unique and set(t.index)==set(a.obs_names);lab[level]=t.reindex(a.obs_names).pred_cell_type.astype(str).to_numpy()
 c=pd.DataFrame({'x_um':coords[:,0],'y_um':coords[:,1],'fov':'one_native_section','native_total':totals,'native_genes':detected,**{'author_'+k:v for k,v in lab.items()}});major=pd.Series(lab['major']).str.casefold();minor=pd.Series(lab['minor']).str.casefold();sub=pd.Series(lab['subset']).str.casefold();geneix={g:i for i,g in enumerate(genes)}
 def pos(g):return x[:,geneix[g]]>=1 if g in geneix else np.zeros(a.n_obs,bool)
 def program(k,n):
  v=[pos(g) for g in PANELS[k] if g in geneix];return np.sum(v,axis=0)>=n if len(v)>=n else np.zeros(a.n_obs,bool)
 ec=program('EC',2);bp=program('B',2);tp=program('T',2);nkp=program('NK_anchor',2)&~tp;conflict=(ec&bp)|(ec&tp)|(bp&tp)|(ec&nkp)|(bp&nkp);qc=(totals>=50)&(detected>=20)&~conflict
 nklabel=(minor.str.contains(r'(?<![a-z])nk(?![a-z])',regex=True)|sub.str.contains(r'(?<![a-z])nk(?![a-z])',regex=True)).to_numpy()
 endothelial=major.str.contains('endo').to_numpy()&ec;lymphatic=endothelial&(program('Lymphatic',2)|minor.str.contains('lymphatic').to_numpy()|sub.str.contains('lymphatic').to_numpy());blood=endothelial&~lymphatic
 b=major.str.fullmatch(r'b[\s_-]*cells?').to_numpy()&bp;t=major.str.contains(r'^t[\s_-]*cell',regex=True).to_numpy()&tp&~nklabel;nk=nklabel&major.str.contains(r'^(?:t[\s_-]*cell|nk)',regex=True).to_numpy()&program('NK_anchor',2)&~tp
 mural=major.str.contains(r'pvl|mural|perivascular|smooth',regex=True).to_numpy()&program('Mural',2)&(pos('RGS5')|pos('PDGFRB')|pos('MYH11'));fibro=major.str.contains(r'caf|fibro',regex=True).to_numpy()&program('Fibroblast',2)
 masks={'B':b,'T':t,'NK':nk,'Blood_EC':blood,'Lymphatic_EC':lymphatic,'Mural':mural,'Fibroblast':fibro};overlap=np.sum(list(masks.values()),axis=0)>1;qc&=~overlap;conflict|=overlap
 for k in masks:masks[k]&=~overlap
 c['is_endothelial_reference']=major.str.contains('endo').to_numpy();c['cell_type']='Other'
 for k,v in masks.items():c.loc[v,'cell_type']=k
 c['QC_keep']=qc;c['cross_lineage_marker_review']=conflict;c.to_csv(out/'cell_identity_view_PRIVATE.tsv.gz',sep='\t',index=False);np.savez_compressed(out/'selected_native_counts_PRIVATE.npz',counts=x,genes=np.array(genes))
 plotrows=[];normalized=np.log1p(x/np.maximum(totals,1)[:,None]*10000)
 for level in ['major','minor']:
  for label in sorted(np.unique(lab[level])):
   mask=(lab[level]==label)&(totals>=50)&(detected>=20)
   for g in needed:
    j=geneix.get(g);plotrows.append({'level':level,'author_label':label,'gene':g,'expected_identity_or_axis':next((k for k,gs in PANELS.items() if g in gs),'vascular_target' if g=='ACKR1' else 'recruitment_axis'),'n_cells':int(mask.sum()),'matrix_gene_available':j is not None,'assay_status':('available_in_matrix' if j is not None else 'assayed_but_unavailable_in_processed_matrix' if panel_gene_names is not None and g in panel_gene_names else 'not_assayed_in_author_panel' if panel_gene_names is not None else 'matrix_unavailable_assay_status_unknown'),'fraction_positive':float((x[mask,j]>=1).mean()) if j is not None and mask.any() else np.nan,'mean_log1p_CP10K':float(normalized[mask,j].mean()) if j is not None and mask.any() else np.nan})
 pd.DataFrame(plotrows).to_csv(out/'author_marker_review.tsv',sep='\t',index=False)
 cc=c.loc[qc].reset_index(drop=True);xx=x[qc];neighborhood(cc,genes,xx).to_csv(out/'neighborhood.tsv',sep='\t',index=False);regional_ES(cc,genes,xx).to_csv(out/'regional_ES.tsv',sep='\t',index=False);edge_expression_null(cc,genes,xx).to_csv(out/'conditional_LR_edges.tsv',sep='\t',index=False)
 receipt={'UTC':datetime.now(timezone.utc).isoformat(),'unit_alias':unit,'status':'PROSPECTIVE_NATIVE_VASCULAR_DESCRIPTIVE_AND_CONDITIONAL_NULL_COMPLETE_NOT_CAUSAL_OR_COHORT_REPLICATION','n_input':a.n_obs,'n_QC':int(qc.sum()),'n_cross_lineage_review_excluded_from_primary_view':int(conflict.sum()),'n_by_broad_identity':cc.cell_type.value_counts().to_dict(),'native_UMI_and_coordinate_checks':True,'coordinate_contract':coordinate_receipt,'RNA_total_semantics':'sum over available native count-matrix genes after author gene filtering; not claimed as unfiltered whole-library totals','native_cell_axis_SHA256':axis_hash,'native_source_file_sha256':z['members'][0]['sha256'],'native_source_archive_member_PRIVATE':z['members'][0]['member'],'T_identity_uses_native_CD3E_CD3G_not_unmeasured_TRBC':True,'NK_identity_uses_NCR1_KLRD1_without_PRF1_FCGR3A_state_gate':True,'original_V4_and_author_labels_unchanged':True,'exact_CM03_IGHD_XCL1_measurable':all(g in symbols for g in ['IGHD','XCL1']),'patient_level_replication_pending':True,'scRNA_double_source_CellChat_support_pending':True,'true_vessel_boundary_contact_not_measured':True}
 (out/'complete.json').write_text(json.dumps(receipt,indent=2)+'\n');a.file.close();return receipt
if __name__=='__main__':print(json.dumps(run(*sys.argv[1:])),flush=True)
