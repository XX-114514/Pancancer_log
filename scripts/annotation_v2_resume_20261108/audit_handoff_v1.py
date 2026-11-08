import pathlib,json,hashlib,sys,datetime
import pandas as pd,numpy as np
h=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
def emit(x):print(json.dumps(x),flush=True)
def sha(p):
 d=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):d.update(b)
 return d.hexdigest()
m=pd.read_csv(h/'MANIFEST.tsv',sep='\t');bad=[]
for r in m.itertuples(index=False):
 p=h/r.path
 if p.stat().st_size!=r.bytes or sha(p)!=r.sha256:bad.append(r.path)
assert not bad,bad
emit({'manifest_verified':len(m)})
L=pd.read_csv(h/'p0/cell_labels_r3.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('cell_id');assert L.index.is_unique
review=pd.read_csv(h/'p1/provenance/repo/inventories/annotation_v2/lineage_r2_reviewed.tsv',sep='\t',dtype=str,keep_default_na=False)
objects=[]; qcs=[]
for fp in sorted((h/'p1/obs').glob('*.obs.tsv.gz')):
 name=fp.name.replace('.obs.tsv.gz',''); rnd,pool=name.split('_',1)
 cols=['obs_name','row','cell_id','leiden_r1','n_genes','unit_id','source','sample_id']
 o=pd.read_csv(fp,sep='\t',usecols=cols,dtype={'obs_name':str,'cell_id':str,'leiden_r1':str,'source':str,'unit_id':str,'sample_id':str},keep_default_na=False)
 assert o.obs_name.is_unique and (o.obs_name==o.cell_id).all(); assert np.array_equal(o.row,np.arange(len(o)))
 with np.load(h/f'p1/emb/{name}.emb.npz',allow_pickle=False) as z:
  assert np.array_equal(z['cell_id'].astype(str),o.obs_name.values)
  assert z['X_emb'].shape==(len(o),30) and np.isfinite(z['X_emb']).all()
 oo=o.set_index('obs_name');ll=L.loc[oo.index]
 if rnd=='r2':
  rr=review[review.pool==pool].set_index('cluster'); assert set(oo.leiden_r1)==set(rr.index)
  assert (ll.origin_pool==pool).all();assert (ll.r2_cluster==oo.leiden_r1).all()
  assert (ll.r2_action==oo.leiden_r1.map(rr.action)).all();assert (ll.r2_label==oo.leiden_r1.map(rr.r2_label)).all()
  dest=ll.r2_action.map(lambda a:'r3_'+a[8:] if a.startswith('move_to_') else 'r3_'+pool if a.startswith('keep') else 'r2_'+a)
  assert (ll.major_label==dest).all()
 else: assert set(oo.index)==set(L.index[L.major_label=='r3_'+pool])
 seed=h/f'p0/consensus/{pool}.seeds.tsv.gz'
 if seed.exists() and (rnd=='r3' or pool not in ['Myeloid','Endothelial']):
  s=pd.read_csv(seed,sep='\t',dtype=str).set_index('cell_id');sc=[f'seed_{x}' for x in range(2026110800,2026110810)]
  assert s.index.is_unique and set(s.index)==set(oo.index) and all(c in s for c in sc)
  assert not s[sc].isna().any().any(); assert (s.loc[oo.index,'run_r1']==oo.leiden_r1).all()
  del s
 if rnd=='r3':
  for c,g in oo.groupby('leiden_r1'):
   # unit_id is the documented analysis/capture unit; also inspect source-qualified sample IDs.
   k=g.source.astype(str)+'::'+g.sample_id.astype(str)
   qcs.append({'pool':pool,'cluster':c,'cells':len(g),'max_unit_share':g.unit_id.value_counts(normalize=True).max(),'max_source_sample_share':k.value_counts(normalize=True).max(),'missing_sample_cells':int((g.sample_id=='').sum())})
 objects.append({'object':name,'cells':len(o),'alignment':'pass'});emit(objects[-1])
for pool in ['T_NK','B_Plasma','Stromal']:
 t=pd.read_csv(h/f'p0/transfer/{pool}.cells.tsv.gz',sep='\t',dtype={'cell_id':str});assert t.cell_id.is_unique
 assert (L.loc[t.cell_id,'major_label']=='r3_'+pool).all()
 assert t.transfer_conf.between(0,1).all()
 assert ((t.transfer_conf>=.6)==(t.transfer_status=='transferred')).all()
pd.DataFrame(qcs).to_csv(out/'r3_sample_concentration.tsv',sep='\t',index=False)
r={'status':'passed','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest_files':len(m),'label_rows':len(L),'objects':objects,'r3_label_sha256':sha(h/'p0/cell_labels_r3.tsv.gz'),'full_expression_values_checked':False}
(out/'acceptance.json').write_text(json.dumps(r,indent=2)+'\n');emit(r)
