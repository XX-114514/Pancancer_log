"""Independent frozen signature, patient-score/design binding and immune-statistic audit."""
from pathlib import Path
from datetime import datetime,timezone
import json,numpy as np,pandas as pd
from scipy.stats import rankdata,chi2
from scipy.special import stdtr
from scipy.io import mmread
r=PROJECT_ROOT/'runs/20261107_v4_spatialCellChat_bulk_v1';out=r/'bulk/analysis_v1';plan=json.loads((r/'config/plan_before_effects.json').read_text());p=plan['bulk'];programs=pd.read_csv(r/'bulk/frozen_signatures/all14_frozen_CMT_signature_genes.tsv',sep='\t');scope=pd.read_csv(r/'bulk/frozen_signatures/all14_signature_eligibility.tsv',sep='\t');allerr=[];scoreframes={};nscore=0
for f in sorted(out.glob('*_scores_clinical_PRIVATE.tsv.gz')):
 co=f.name.replace('_scores_clinical_PRIVATE.tsv.gz','');keycol=pd.read_csv(f,sep='\t',nrows=0).columns[0];ss=pd.read_csv(f,sep='\t',index_col=0,dtype={keycol:str});xx=pd.read_pickle(out/(co+'_selected_expression_PRIVATE.pkl.gz'));xx.columns=xx.columns.astype(str);assert ss.index.is_unique and xx.columns.is_unique and ss.index.equals(xx.columns)
 arr=xx.to_numpy(float);ranks=rankdata(arr,method='average',axis=1,nan_policy='omit')/np.isfinite(arr).sum(1)[:,None]
 for name,genes in list((cm,programs.loc[programs.cm.eq(cm),'gene'].tolist()) for cm in plan['all_CM'])+list(p['immune_panels'].items()):
  require=10 if name.startswith('CM') else len(genes);eligible=len(genes)==require and set(genes)<=set(xx.index) and np.isfinite(xx.loc[genes]).all().all()
  expected=ranks[xx.index.get_indexer(genes)].mean(0) if eligible else np.full(len(ss),np.nan);actual=ss[name].to_numpy(float);assert np.allclose(expected,actual,equal_nan=True,atol=1e-12), (co,name)
  nscore+=len(ss);allerr.extend(np.abs(expected[np.isfinite(expected)]-actual[np.isfinite(expected)]))
 scoreframes[co]=ss
rows=pd.read_csv(out/'all_frozen_model_slots_private.tsv',sep='\t');designs=json.loads((out/'fitted_designs_private.json').read_text())
for entry in designs:
 row=rows.iloc[entry['record_index']];co=('TCGA_'+row.cohort) if row.family=='TCGA_survival' else row.cohort;ss=scoreframes[co];d=pd.DataFrame(index=ss.index);d['joint']=ss[row.axis];d['pan_T']=ss.pan_T
 for name in ['age','sex','stage']:
  if name in ss:d[name]=ss[name]
 if entry['kind']=='logistic':d['y']=ss.response;required=['y']
 else:d['t']=ss[row.endpoint+'.time'] if row.family=='TCGA_survival' else ss[row.endpoint];d['e']=ss[row.endpoint] if row.family=='TCGA_survival' else ss[row.endpoint+'_event'];required=['t','e']
 used=d[required+entry['predictors']].replace([np.inf,-np.inf],np.nan).dropna();used=used[used.y.isin([0,1])] if entry['kind']=='logistic' else used[(used.t>0)&used.e.isin([0,1])];got=pd.DataFrame(entry['data']);assert np.allclose(used[got.columns],got,atol=1e-12),entry['record_index']
imm=pd.read_csv(out/'all_frozen_immune_slots_private.tsv',sep='\t');errs=[]
for row in imm[imm.status.eq('tested')].itertuples():
 ss=scoreframes['TCGA_'+row.cohort]
 if row.endpoint=='C1_C6_omnibus':
  d=ss[[row.cm,'immune_subtype']].dropna();d=d[d.immune_subtype.isin(['C1','C2','C3','C4','C5','C6'])];rk=rankdata(d[row.cm]);N=len(d);H=0.
  for _,ii in d.groupby('immune_subtype').indices.items():H+=rk[ii].sum()**2/len(ii)
  _,counts=np.unique(d[row.cm],return_counts=True);tie=1-np.sum(counts**3-counts)/(N**3-N);stat=(12*H/(N*(N+1))-3*(N+1))/tie;pv=chi2.sf(stat,d.immune_subtype.nunique()-1)
 else:
  d=ss[[row.cm,row.endpoint]].dropna();a=rankdata(d[row.cm]);b=rankdata(d[row.endpoint]);stat=np.corrcoef(a,b)[0,1];df=len(d)-2;pv=2*stdtr(df,-abs(stat)*np.sqrt(df/max(1-stat**2,1e-300)))
 assert abs(stat-row.statistic)<1e-8 and abs(pv-row.p)<1e-8;errs.append(max(abs(stat-row.statistic),abs(pv-row.p)))
from statsmodels.stats.multitest import multipletests
assert len(imm)==4620 and len(rows)==2100
for frame in [rows,imm]:
 for family,g in frame.groupby('family'):
  for column,method in [('BH_q','fdr_bh'),('BY_q','fdr_by')]:assert np.allclose(multipletests(g.p.fillna(1),method=method)[1],g[column],atol=1e-12)
# Independent selected-gene fold-change and top10 ordering against saved sample-level DE results.
x=mmread(r/'bulk/signature_inputs/mean_CP10K_PRIVATE.mtx').toarray();genes=pd.read_csv(r/'bulk/signature_inputs/genes.tsv',header=None)[0];meta=pd.read_csv(r/'bulk/signature_inputs/samples_PRIVATE.tsv',sep='\t',index_col=0);fce=[]
for z in scope.itertuples():
 if z.n_high_samples<3:continue
 de=pd.read_csv(r/'bulk/frozen_signatures'/(z.cm+'_CMT_DE.tsv'),sep='\t');fc='avg_log2FC' if 'avg_log2FC' in de else 'avg_logFC';mask=meta.CMT.eq(z.cm).to_numpy();ix=pd.Index(genes).get_indexer(de.gene);expected=np.log2(((x[ix][:,mask].sum(1)+1)/mask.sum())/((x[ix][:,~mask].sum(1)+1)/(~mask).sum()));assert np.allclose(expected,de[fc],atol=1e-8),z.cm;fce.extend(np.abs(expected-de[fc]));dd=de[(de.p_val_adj<.05)&(de[fc]>.25)].sort_values([fc,'gene'],ascending=[False,True]).head(10);assert dd.gene.tolist()==programs.loc[programs.cm.eq(z.cm),'gene'].tolist()
receipt={'status':'PASS_all_frozen_patient_scores_model_design_bindings_immune_statistics_full_FDR_and_signature_FC','UTC':datetime.now(timezone.utc).isoformat(),'n_patient_score_entries_checked':nscore,'n_model_designs_bound':len(designs),'n_immune_slots':len(imm),'n_immune_tested':len(errs),'maximum_score_error':float(max(allerr)) if len(allerr) else 0,'maximum_immune_statistic_p_error':float(max(errs)) if errs else 0,'maximum_signature_FC_error':float(max(fce)) if fce else 0,'independent_collection_or_clinical_validation_claim':False};(r/'audit/bulk_measurement_audit_v2.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
