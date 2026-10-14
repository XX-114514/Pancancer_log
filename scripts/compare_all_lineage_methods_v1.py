"""Compare every cNMF component with every permitted Spectra factor; no outcome selection."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,argparse
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def frame(p):
 with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def unit(x):
 x=np.asarray(x,float);d=np.linalg.norm(x,axis=1,keepdims=True)
 return np.divide(x,d,out=np.full(x.shape,np.nan),where=d>0)
def residual(x,g):
 y=np.asarray(x,float).copy()
 for v in np.unique(g):
  m=g==v;a=y[m];b=a-a.mean(axis=0);b[:,np.ptp(a,axis=0)==0]=0;y[m]=b
 return y
def ranked(x):
 y=rankdata(x,axis=0);y-=y.mean(axis=0);return unit(y.T)
def same(a,b):return bool(np.isclose(a,b,atol=2e-10,rtol=2e-10,equal_nan=True))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--major',required=True);ap.add_argument('--plan',required=True);args=ap.parse_args();P=Path(args.plan);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256']
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(A)==plan['accepted_input_inventory_sha256'];majors=sorted(r['major'] for r in json.loads(A.read_text())['files']);major=args.major;i=majors.index(major);block=next(b for b,mm in plan['blocks'].items() if major in mm)
 CD=(R/'all_lineage_cNMF_diagnostics_20261014_v1/major_00') if i==0 else V/'all_lineage_cNMF_diagnostics_20261014_v2'/f'major_{i:02}'
 CA=(R/'all_lineage_cNMF_diagnostic_audit_20261014_v2/major_00') if i==0 else V/'all_lineage_cNMF_diagnostic_audit_20261014_v3'/f'major_{i:02}'
 SD=(R/'all_lineage_Spectra_diagnostics_20261014_v1'/block) if block=='lymphoid_T' else V/'all_lineage_Spectra_diagnostics_20261014_v2'/block
 SA=(R/'all_lineage_Spectra_diagnostic_audit_20261014_v1'/block) if block=='lymphoid_T' else V/'all_lineage_Spectra_diagnostic_audit_20261014_v2'/block
 for d,a in [(CD,CA),(SD,SA)]:
  audit=json.loads((a/'independent_audit.json').read_text());assert audit['status'].startswith('PASS_ALL_') and audit['producer_inventory_sha256']==sha(d/'inventory.json')
  iv=json.loads((d/'inventory.json').read_text())
  for n,h in iv['output_sha256'].items():assert sha(d/n)==h
 ci=json.loads((CD/'inventory.json').read_text());si=json.loads((SD/'inventory.json').read_text());ct=pd.read_csv(CD/'all45_candidates.tsv',sep='\t');st=pd.read_csv(SD/'all_candidate_diagnostics.tsv',sep='\t',keep_default_na=False);st=st[st.scope.isin(['block_global',major])].copy();ids=st.anchor_factor.to_numpy(int);smroot=V/'all_lineage_Spectra_state_discovery_20261014_v2';models=[];metadata=[]
 for seed in plan['seeds']:
  p=smroot/block/f'seed{seed}';f=p/'factors_and_discovery_scores.npz';md=p/'private_discovery_metadata.tsv.gz'
  for x in [f,md]:assert sha(x)==si['input_seed_sha256'][str(x.relative_to(smroot))]
  with np.load(f,allow_pickle=False) as z:models.append({k:z[k].copy() for k in z.files})
  metadata.append(pd.read_csv(md,sep='\t',index_col=0,keep_default_na=False))
 m=metadata[0];mask=m.major.eq(major).to_numpy();meta=m.loc[mask];source=meta.gse_id.to_numpy(str);node=(meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy();maps=[ids,st.seed2_factor.to_numpy(int),st.seed3_factor.to_numpy(int)]
 for z,md in zip(models,metadata):assert md.equals(m) and np.array_equal(z['cell_ids'],m.index.to_numpy(str)) and np.array_equal(z['genes'],models[0]['genes'])
 sw=np.median(np.stack([unit(z['factors'])[ix] for z,ix in zip(models,maps)]),axis=0);sg=models[0]['genes'];ss=[z['cell_scores'][mask][:,ix] for z,ix in zip(models,maps)];stop=[set(str(v).split(';')) for v in st.primary50_median_L2_gene_weights];rows=[];checks=[];base=V/'all_lineage_cNMF_state_contrast_20261014_v2'/f'major_{i:02}'/'cnmf/cnmf_tmp'
 for k in plan['K']:
  cp=base/f'cnmf.gene_spectra_score.k_{k}.dt_0_5.df.npz';up=base/f'cnmf.usages.k_{k}.dt_0_5.consensus.df.npz'
  for p in [cp,up]:assert sha(p)==ci['input_output_sha256'][str(p.relative_to(base.parent.parent))]
  c=frame(cp);u=frame(up);assert u.index.equals(meta.index) and c.index.equals(u.columns);assert (u.values.sum(axis=1)>0).all();cu=u.values/u.values.sum(axis=1,keepdims=True);cw=np.maximum(c.values,0);cg=c.columns.to_numpy(str);loc=pd.Index(cg).get_indexer(sg);assert (loc>=0).all();common=cw[:,loc]
  cos=unit(common)@unit(sw).T;full=np.divide(common@sw.T,np.linalg.norm(cw,axis=1)[:,None]*np.linalg.norm(sw,axis=1)[None,:]);ret=np.divide(np.linalg.norm(common,axis=1),np.linalg.norm(cw,axis=1));cr={};sr={};rhos={}
  for label,g in [('pooled',None),('source',source),('source_node',node)]:
   cr[label]=cu if g is None else residual(cu,g);sr[label]=[x if g is None else residual(x,g) for x in ss];rhos[label]=np.stack([ranked(cr[label])@ranked(x).T for x in sr[label]])
  # Direct scalar dot products for every cosine, plus deterministic score checks against scipy.
  for a in range(k):
   for b in range(len(ids)):
    den=np.linalg.norm(common[a])*np.linalg.norm(sw[b]);v=float(np.dot(common[a],sw[b])/den) if den else np.nan;assert same(v,cos[a,b])
  for ix in np.unique(np.linspace(0,k*len(ids)-1,min(20,k*len(ids)),dtype=int)):
   a,b=divmod(int(ix),len(ids))
   for label in rhos:
    for seedidx in range(3):
     x=cr[label][:,a];y=sr[label][seedidx][:,b];v=float(spearmanr(x,y).statistic) if np.ptp(x)>0 and np.ptp(y)>0 else np.nan;assert same(v,rhos[label][seedidx,a,b]);checks.append(dict(K=k,component=int(c.index[a]),anchor_factor=int(ids[b]),seed=plan['seeds'][seedidx],score_context=label,scalar_rho=v,matrix_rho=float(rhos[label][seedidx,a,b])))
  lookup=ct[ct.K.eq(k)].set_index('component')
  for a,component in enumerate(c.index):
   cc=lookup.loc[component];ctop=set(str(cc.primary50_positive_coefficient_genes).split(';'))
   for b,sp in enumerate(st.itertuples()):
    rec=dict(major=major,cNMF_id=cc.candidate_id,Spectra_id=sp.candidate_id,K=k,component=int(component),anchor_factor=int(ids[b]),Spectra_scope=sp.scope,cNMF_repeat_crossK_stable=bool(cc.repeat_and_crossK_candidate),Spectra_context_stable=bool(sp.stable_context_candidate),common_gene_cosine=float(cos[a,b]),native_axis_padded_cosine=float(full[a,b]),cNMF_L2_fraction_on_Spectra_gene_axis=float(ret[a]),top50_jaccard=len(ctop&stop[b])/len(ctop|stop[b]),cells=len(meta),sources=len(set(source)),formal_state=False)
    for label in rhos:
     vals=rhos[label][:,a,b];rec[label+'_rho_seed1']=float(vals[0]);rec[label+'_rho_seed2']=float(vals[1]);rec[label+'_rho_seed3']=float(vals[2]);rec[label+'_rho_min']=float(np.min(vals));rec[label+'_rho_median']=float(np.median(vals))
    rows.append(rec)
 out=V/'all_lineage_cross_method_comparison_20261014_v1'/f'major_{i:02}';assert not out.exists();out.mkdir(parents=True);pd.DataFrame(rows).to_csv(out/'all_method_pairs.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(out/'scalar_score_checks.tsv',sep='\t',index=False)
 assert len(rows)==45*len(ids)
 inv=dict(status='COMPLETE_ALL_PERMITTED_METHOD_PAIRS_NOT_STATE_VALIDATION',major=major,block=block,observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),cNMF_inventory_sha256=sha(CD/'inventory.json'),Spectra_inventory_sha256=sha(SD/'inventory.json'),pairs=len(rows),permitted_Spectra_factors=len(ids),scalar_score_checks=len(checks),all_cosines_scalar_checked=True,all_factors_including_unstable=True,clinical_LR_p_used=False,formal_states=False,output_sha256={p.name:sha(p) for p in out.glob('*.tsv')});(out/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n');print(json.dumps(inv),flush=True)
if __name__=='__main__':main()
