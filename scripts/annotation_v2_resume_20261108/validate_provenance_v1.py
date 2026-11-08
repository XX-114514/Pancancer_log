import pathlib,sys,json
import pandas as pd,numpy as np
h,r=map(pathlib.Path,sys.argv[1:3]);out=r/'working_annotation_v1'
L=pd.read_csv(h/'p0/cell_labels_r3.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('cell_id')
cols=list(L.columns);done=0
for d in pd.read_csv(out/'cell_labels_working_v1.tsv.gz',sep='\t',usecols=['cell_id']+cols,dtype=str,keep_default_na=False,chunksize=200000):
 d=d.set_index('cell_id');assert d[cols].equals(L.loc[d.index,cols]);done+=len(d)
assert done==len(L);del L
D=pd.read_csv(out/'cell_labels_working_v1.tsv.gz',sep='\t',usecols=['cell_id','transfer_label_raw','transfer_conf','annotation_status','major_label'],keep_default_na=False).set_index('cell_id')
seen=set();counts={}
for pool in ['T_NK','B_Plasma','Stromal','Microglia','Endothelial','Myeloid']:
 base=h/'p0/transfer' if pool in ['T_NK','B_Plasma','Stromal'] else r/'transfer_v1'
 t=pd.read_csv(base/f'{pool}.cells.tsv.gz',sep='\t',keep_default_na=False).set_index('cell_id');assert t.index.is_unique and not(seen&set(t.index));seen.update(t.index)
 x=D.loc[t.index];assert (x.transfer_label_raw==t.transfer_label).all();assert np.allclose(x.transfer_conf.astype(float),t.transfer_conf.astype(float),atol=1e-7)
 assert (x.annotation_status==np.where(t.transfer_status=='transferred','transferred','unresolved_low_depth')).all()
 assert (x.major_label=='r3_'+pool).all();counts[pool]=len(t)
assert seen==set(D.index[D.annotation_status.isin(['transferred','unresolved_low_depth'])])
a={'status':'passed','original_columns_preserved':cols,'rows_compared':done,'transfer_rows_by_pool':counts,'transfer_rows_total':len(seen),'lineage_pool_not_silently_changed':True}
(out/'provenance_acceptance.json').write_text(json.dumps(a,indent=2));print(json.dumps(a),flush=True)
