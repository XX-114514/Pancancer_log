"""Run documented transfer against verified exported obs and X_emb, no expression reads."""
import sys,pathlib,json
import numpy as np,pandas as pd,pynndescent
h,r=map(pathlib.Path,sys.argv[1:3]);pool=sys.argv[3]
assert json.loads((r/'input_audit_v1/acceptance.json').read_text())['status']=='passed'
out=r/'transfer_v1';out.mkdir(exist_ok=True)
assert not (out/f'{pool}.cells.tsv.gz').exists()
rnd='r2' if pool=='Microglia' else 'r3';name=rnd+'_'+pool
rvpath=h/'p1/provenance/repo/inventories/annotation_v2/lineage_r2_reviewed.tsv' if rnd=='r2' else r/'review_v1/lineage_r3_reviewed.tsv'
R=pd.read_csv(rvpath,sep='\t',dtype=str,keep_default_na=False);R=R[R.pool==pool].set_index('cluster')
o=pd.read_csv(h/f'p1/obs/{name}.obs.tsv.gz',sep='\t',usecols=['obs_name','row','leiden_r1','n_genes'],dtype={'obs_name':str,'leiden_r1':str}).set_index('obs_name')
L=pd.read_csv(h/'p0/cell_labels_r3.tsv.gz',sep='\t',usecols=['cell_id','major_label','r1_action'],dtype=str,keep_default_na=False).set_index('cell_id').loc[o.index]
with np.load(h/f'p1/emb/{name}.emb.npz',allow_pickle=False) as z:
 assert np.array_equal(z['cell_id'].astype(str),o.index.values);Z=z['X_emb'].astype(np.float32)
act=o.leiden_r1.map(R.action).values;lab=o.leiden_r1.map(R[rnd+'_label']).values
assert pd.notna(lab).all()
if pool=='Microglia':lab=np.array(['Microglia' if x.startswith('Microglia') else 'Mac_glioma' if x.startswith('Mac_glioma') else x for x in lab])
low=(L.r1_action.values=='keep_flag_low_depth')|(act=='keep_flag_low_depth')
refacts=['keep','keep_tissue']+(['move_to_Myeloid'] if pool=='Microglia' else [])
ref=np.isin(act,refacts)&~low
q=low&(L.major_label.values=='r3_'+pool)&np.char.startswith(act.astype(str),'keep')
assert ref.sum()>31 and q.sum()>0 and not (ref&q).any()
seed=2026110704;K=30

def vote(zr,yr,zq):
 ix=pynndescent.NNDescent(zr,n_neighbors=K+1,random_state=seed,n_jobs=4);ix.prepare()
 i,d=ix.query(zq,k=K);w=1/(d+1e-6);w/=w.sum(1,keepdims=True)
 codes,names=pd.factorize(yr);s=np.zeros((len(zq),len(names)),np.float32)
 np.add.at(s,(np.arange(len(zq))[:,None].repeat(K,1),codes[i]),w)
 return names[s.argmax(1)],s.max(1)
g=o.n_genes.values
vd=pd.DataFrame({'lab':lab[ref],'g':g[ref],'i':np.flatnonzero(ref)});hold=vd.loc[vd.groupby('lab').g.rank(method='first',pct=True)<=.10,'i'].values
train=np.setdiff1d(np.flatnonzero(ref),hold)
print(json.dumps({'pool':pool,'reference':int(ref.sum()),'query':int(q.sum()),'holdout':len(hold),'step':'validation'}),flush=True)
p,cf=vote(Z[train],lab[train],Z[hold]);V=pd.DataFrame({'true':lab[hold],'pred':p,'conf':cf,'genes':g[hold]});V['correct']=V.true==V.pred
val=V.groupby('true').agg(cells=('pred','size'),accuracy=('correct','mean'),median_genes=('genes','median'));val['accuracy_if_confident']=V[V.conf>=.6].groupby('true').correct.mean()
val.to_csv(out/f'{pool}.validation.tsv',sep='\t');print(json.dumps({'pool':pool,'validation_macro_accuracy':float(val.accuracy.mean()),'step':'query'}),flush=True)
pq,cq=vote(Z[ref],lab[ref],Z[q]);T=pd.DataFrame({'cell_id':o.index[q],'pool':pool,'review_round':rnd,'cluster':o.leiden_r1.values[q],'cluster_label':lab[q],'cluster_action':act[q],'transfer_label':pq,'transfer_conf':cq.round(3)})
T['transfer_status']=np.where(T.transfer_conf>=.6,'transferred','unresolved');T.loc[T.transfer_status=='unresolved','transfer_label']=pool+'_unresolved_low_depth'
T.to_csv(out/f'{pool}.cells.tsv.gz',sep='\t',index=False,compression='gzip');T.groupby(['transfer_label','transfer_status']).size().rename('cells').reset_index().to_csv(out/f'{pool}.summary.tsv',sep='\t',index=False)
rpt={'pool':pool,'reference':int(ref.sum()),'query':int(q.sum()),'unresolved':int((T.transfer_status=='unresolved').sum()),'validation_macro_accuracy':float(val.accuracy.mean()),'seed':seed,'k':K,'min_conf':.6,'threads':4,'scope':'inferred sensitivity only; validation is a higher-depth upper bound','microglia_mac_predictions_do_not_change_pool':pool=='Microglia'}
(out/f'{pool}.DONE.json').write_text(json.dumps(rpt,indent=2));print(json.dumps(rpt),flush=True)
