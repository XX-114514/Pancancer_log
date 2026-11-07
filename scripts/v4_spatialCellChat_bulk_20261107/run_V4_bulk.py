"""Frozen same-V4 gene-program associations; no cutpoints or outcome-trained genes."""
import sys,json,gzip,hashlib,warnings
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import spearmanr,kruskal
from clinical_helpers_v4 import load_ici,load_tcga,helpers,collapse,fit_model,bh
r=PROJECT_ROOT/'runs/20261107_v4_spatialCellChat_bulk_v1';plan=json.loads((r/'config/plan_before_effects.json').read_text());p=plan['bulk'];out=r/'bulk/analysis_v1';out.mkdir(exist_ok=False)
sig=pd.read_csv(r/'bulk/frozen_signatures/all14_frozen_CMT_signature_genes.tsv',sep='\t');elig=pd.read_csv(r/'bulk/frozen_signatures/all14_signature_eligibility.tsv',sep='\t');program={cm:sig.loc[sig.cm.eq(cm),'gene'].tolist() for cm in plan['all_CM']};valid={cm:bool(elig.set_index('cm').loc[cm,'status']=='FROZEN10GENE_PROGRAM' and cm not in plan['empty_network_CM']) for cm in plan['all_CM']};panels=p['immune_panels'];wanted=set(sig.gene)|set(sum(panels.values(),[]));base=PROJECT_ROOT.parent/'Project_v5/bulk_cohorts_20260904';h=helpers(base);metadata=[];rows=[];designs=[];immune=[]
def scores(expr):
 ranks=expr.rank(axis=1,method='average',pct=True);sc=pd.DataFrame(index=expr.columns)
 for cm,gg in program.items():sc[cm]=ranks.loc[gg].mean() if len(gg)==10 and set(gg)<=set(ranks.index) and ranks.loc[gg].notna().all().all() else np.nan
 for name,gg in panels.items():sc[name]=ranks.loc[gg].mean() if set(gg)<=set(ranks.index) and ranks.loc[gg].notna().all().all() else np.nan
 return sc

def save_data(co,expr,c,sc,meta):
 expr.loc[expr.index.intersection(sorted(wanted))].to_pickle(out/(co+'_selected_expression_PRIVATE.pkl.gz'),compression='gzip');sc.join(c,rsuffix='_clinical').to_csv(out/(co+'_scores_clinical_PRIVATE.tsv.gz'),sep='\t');metadata.append(dict(meta,cohort=co));print(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'stage':'cohort_scored','cohort':co,'n_patients':len(c)}),flush=True)

def clinical(co,cm,sc,c,ep,model,tcga=False):
 family='TCGA_survival' if tcga else 'ICI_response' if ep=='response' else 'ICI_survival';d=pd.DataFrame(index=c.index);d['joint']=sc[cm];d['pan_T']=sc.pan_T
 if tcga:
  d['age']=c.age;d['sex']=c.sex;d['stage']=c.stage;d['t']=c[ep+'.time'];d['e']=c[ep];predictors=['joint','pan_T','age']+(['sex'] if c.sex.nunique()>1 else [])+(['stage'] if model=='stage_sensitivity' else [])
 else:
  d['y']=c.response;d['t']=c.get(ep,pd.Series(np.nan,index=c.index));d['e']=c.get(ep+'_event',pd.Series(np.nan,index=c.index));predictors=['joint','pan_T']
 kind='logistic' if ep=='response' else 'Cox';used=None
 if not valid[cm]:fit={'status':'unavailable','reason':'original_empty_network_or_insufficient_frozen_CMT_signature','n':len(c),'events':0}
 elif tcga and model=='stage_sensitivity' and (len(c)==0 or c.stage.notna().mean()<.7):fit={'status':'unavailable','reason':'overall_stage_availability_below70pct','n':len(c),'events':0}
 else:fit,used=fit_model(d,kind,predictors,tcga)
 rows.append(dict(family=family,cohort=co,axis=cm,endpoint=ep,model=model,**fit))
 if used is not None:designs.append(dict(record_index=len(rows)-1,kind=kind,predictors=predictors,data=used.to_dict(orient='list')))

for co in p['ICI_cohorts']:
 try:
  if co=='SU2C_MARK':
   source=r/'sources/SU2C_open_author_full';f=source/'SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct';assert f.stat().st_size==44557830 and hashlib.sha256(f.read_bytes()).hexdigest()=='892f8687ed42eb017894aa9c99d4d4a04f5af52037f88709ee7913df36166160'
   x=pd.read_csv(f,sep='\t',skiprows=2);gene=x.Description.astype(str);assert len(set(gene)&wanted)>=10;expr=x.iloc[:,2:].copy();expr.index=gene;expr=expr.loc[~expr.index.duplicated(keep=False)];expr=np.log2(expr.clip(lower=0)+1)
   cc=pd.read_excel(source/'Table_S1_Clinical_Annotations.xlsx',header=2);pid='Harmonized_SU2C_Participant_ID_v2';sid='Harmonized_SU2C_RNA_Tumor_Sample_ID_v2';cc=cc[cc.RNA_All.eq(1)&cc['Pre-treatment_RNA_Sample_QC'].eq('Keep')].set_index(sid);assert cc.index.is_unique;assert set(cc.index)<=set(expr.columns);expr=expr[cc.index];c=pd.DataFrame(index=cc.index);c['response']=cc.Harmonized_Confirmed_BOR.map({'CR':1,'PR':1,'SD':0,'PD':0})
   for ep in ['OS','PFS']:c[ep]=pd.to_numeric(cc['Harmonized_'+ep+'_Days'],errors='coerce');v=pd.to_numeric(cc['Harmonized_'+ep+'_Event'],errors='coerce');c[ep+'_event']=v.where(v.isin([0,1]))
   expr,c,meta=collapse(expr,c,cc[pid]);meta['time_status']='author_Keep_pretreatment_RNA'
  else:expr,c,meta=load_ici(base,co,h)
  expr=expr.loc[expr.index.intersection(sorted(wanted))];sc=scores(expr);save_data(co,expr,c,sc,meta)
 except (FileNotFoundError,AssertionError) as ex:
  # Keep the requested cohort and slots. No partial old panel or approximate clinical mapping.
  c=pd.DataFrame(columns=['response','OS','OS_event','PFS','PFS_event']);sc=pd.DataFrame(columns=plan['all_CM']+list(panels));metadata.append({'cohort':co,'input_status':'UNAVAILABLE_EXACT_FULL_INPUT_OR_KEY_REQUIREMENT','exception_type':type(ex).__name__})
 for cm in plan['all_CM']:
  for ep in p['ICI_endpoints']:clinical(co,cm,sc,c,ep,'pan_T_adjusted')

expr,c,meta=load_tcga(base,wanted);metadata.append(meta);sub=pd.read_csv(base/'tcga/TCGA_immune_subtype.txt.gz',sep='\t',dtype=str);sub['patient']=sub['sample'].str[:12];cats=sub.groupby('patient').Subtype_Immune_Model_Based.agg(lambda x:x.dropna().unique().tolist());cats=cats[cats.map(len).eq(1)].map(lambda x:x[0]);c['immune_subtype']=cats.reindex(c.index)
for co in p['TCGA_cancers']:
 cc=c[c.cancer.eq(co)].copy();xx=expr[cc.index];sc=scores(xx);save_data('TCGA_'+co,xx,cc,sc,{'n_patients':len(cc),'scope':'01_primary_solid_tumour','purity_adjustment':'not_available_in_this_input'})
 for cm in plan['all_CM']:
  for ep in p['TCGA_endpoints']:
   for model in p['TCGA_models']:clinical(co,cm,sc,cc,ep,model,True)
  for endpoint in ['C1_C6_omnibus']+[x for x in panels if x!='pan_T']:
   rec={'family':'TCGA_immune','cohort':co,'cm':cm,'endpoint':endpoint,'n':0,'status':'unavailable','reason':'','statistic':np.nan,'p':np.nan}
   if not valid[cm]:rec['reason']='original_empty_network_or_insufficient_frozen_signature'
   elif endpoint=='C1_C6_omnibus':
    d=pd.DataFrame({'score':sc[cm],'group':cc.immune_subtype}).dropna();d=d[d.group.isin(['C1','C2','C3','C4','C5','C6'])];gs=[t.score.to_numpy() for _,t in d.groupby('group')];rec['n']=len(d)
    if len(d)>=60 and len(gs)>=2 and min(map(len,gs))>=10 and d.score.nunique()>1:
     v=kruskal(*gs);rec.update(status='tested',statistic=float(v.statistic),p=float(v.pvalue),n_subtypes=len(gs))
    else:rec['reason']='insufficient_subtype_or_patient_counts'
   elif set(program[cm])&set(panels[endpoint]):rec['reason']='signature_endpoint_gene_overlap_not_an_independent_immune_assay'
   else:
    d=pd.DataFrame({'score':sc[cm],'immune':sc[endpoint]}).dropna();rec['n']=len(d)
    if len(d)>=60 and d.score.nunique()>1 and d.immune.nunique()>1:
     v=spearmanr(d.score,d.immune);rec.update(status='tested',statistic=float(v.statistic),p=float(v.pvalue))
    else:rec['reason']='missing_measurement_or_insufficient_variation_patients'
   immune.append(rec)
 (out/'progress.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'last_cancer':co,'clinical_slots_written':len(rows),'immune_slots_written':len(immune)}))
result=pd.DataFrame(rows);imm=pd.DataFrame(immune)
for frame in [result,imm]:
 if 'p' not in frame:frame['p']=np.nan
 for _,g in frame.groupby('family'):
  frame.loc[g.index,'BH_q']=bh(g.p.fillna(1));frame.loc[g.index,'BY_q']=bh(g.p.fillna(1),True)
assert len(result)==2100 and len(imm)==4620
result.to_csv(out/'all_frozen_model_slots_private.tsv',sep='\t',index=False);imm.to_csv(out/'all_frozen_immune_slots_private.tsv',sep='\t',index=False);(out/'fitted_designs_private.json').write_text(json.dumps(designs));(out/'cohort_input_audit.json').write_text(json.dumps(metadata,indent=2));receipt={'status':'ALL2100_CLINICAL_AND4620_IMMUNE_SLOTS_WRITTEN_PENDING_INDEPENDENT_AUDIT','UTC':datetime.now(timezone.utc).isoformat(),'n_slots':len(result),'n_immune_slots':len(imm),'n_fitted':len(designs),'program_not_exact_NMF_activity':True,'independent_patient_or_prediction_claim':False};(out/'terminal_public.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
