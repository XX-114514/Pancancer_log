"""Synthetic full native-format test; no human records or biological conclusions."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from run_native_vascular_controls_v5 import run,PANELS
base=Path(sys.argv[1]);base.mkdir(parents=True,exist_ok=False);native=base/'native';maps=base/'mapping';native.mkdir();maps.mkdir()
types=['Blood_EC']*60+['B']*40+['T']*40+['NK']*40+['Mural']*30+['Fibroblast']*30+['Other']*150
symbols=sorted(set(sum(PANELS.values(),[])+['CCL19','CCL21','CCR7','CXCL12','CXCR4']+[f'SYNTH_BG{i}' for i in range(25)])-{'COL1A1','COL1A2','GNLY'});idx={g:i for i,g in enumerate(symbols)};x=np.zeros((len(types),len(symbols)),np.float32);xy=[];major=[];minor=[]
for i,t in enumerate(types):
 for g in [f'SYNTH_BG{k}' for k in range(25)]:x[i,idx[g]]=2
 gs={'Blood_EC':['PECAM1','VWF','CDH5'],'B':['MS4A1','CD79A','CD79B'],'T':['CD3E','CD3G'],'NK':['NCR1','KLRD1'],'Mural':['RGS5','PDGFRB'],'Fibroblast':['DCN','PDGFRA'],'Other':[]}[t]
 for g in gs:x[i,idx[g]]=2
 # Deliberately noPRF1/FCGR3A expression in NK: identity must survive without effector state.
 if t=='Blood_EC' and i<30:x[i,idx['ACKR1']]=1
 if t=='Mural' and i%2==0:x[i,idx['CCL19']]=1
 if t=='Fibroblast' and i%2==0:x[i,idx['CXCL12']]=1
 if t in ['B','T'] and i%2==0:x[i,idx['CCR7']]=1;x[i,idx['CXCR4']]=1
 side=(i//30)%2;xy.append([10+1000*side+i%7,20+i%13]);major.append({'Blood_EC':'Endothelial','B':'B-cells','T':'T-cells','NK':'T-cells','Mural':'PVL','Fibroblast':'CAFs','Other':'Other'}[t]);minor.append({'Blood_EC':'Endothelial ACKR1','B':'B cells Memory','T':'T cells CD4+','NK':'NK cells','Mural':'PVL Immature','Fibroblast':'CAFs MSC iCAF-like','Other':'Other'}[t])
ids=pd.Index([f'SYNTHETIC_CELL_{i:04d}' for i in range(len(types))]);a=ad.AnnData(sparse.csr_matrix(np.log1p(x)),obs=pd.DataFrame({'x_centroid':np.array(xy)[:,0],'y_centroid':np.array(xy)[:,1]},index=ids),var=pd.DataFrame(index=symbols));a.obsm['xenium_um_spatial']=np.array(xy);a.obsm['spatial']=np.array(xy)*2 # Different image pixels must not replace native microns.
f=native/'synthetic_PRIVATE.h5ad';a.write_h5ad(f);v={'member':'anndata_objects/xenium/SYNTH.h5ad','local_name':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()};(native/'complete_PRIVATE.json').write_text(json.dumps({'members':[v]}));rows=[]
for level,labels in [('major',major),('minor',minor),('subset',minor)]:
 f=maps/(level+'_PRIVATE.csv');pd.DataFrame({'pred_cell_type':labels,'pred_cell_type_prob':1e-5},index=ids).to_csv(f);rows.append({'member':f'02_cell_type_mapping/label_transfer/xenium/{level}_level/SYNTH_{level}_pred_cell_types.csv','local_name':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
(maps/'complete_PRIVATE.json').write_text(json.dumps({'members':rows}));z=run(base,'SYNTHETIC_ONLY',native,maps)
assert z['n_QC']==len(types) and z['n_by_broad_identity']['NK']==40 and z['n_by_broad_identity']['T']==40 and z['native_UMI_and_coordinate_checks']
out=base/'spatial_controls_v1/SYNTHETIC_ONLY';assert len(pd.read_csv(out/'neighborhood.tsv',sep='\t'))==12 and len(pd.read_csv(out/'regional_ES.tsv',sep='\t'))==27 and len(pd.read_csv(out/'conditional_LR_edges.tsv',sep='\t'))==27
print('PASS_SYNTHETIC_ONLY: raw log1p counts, native/image coordinate distinction, all3 barcode joins, T CD3E/G, NK without effectors, complete12+27+27 outputs')
