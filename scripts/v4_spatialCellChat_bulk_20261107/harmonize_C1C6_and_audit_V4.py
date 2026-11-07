"""Exact source-category harmonization; complete immune family recomputed, clinical models reused."""
from pathlib import Path
from datetime import datetime,timezone
import json,re,hashlib,numpy as np,pandas as pd
from scipy.stats import kruskal,rankdata,chi2
from statsmodels.stats.multitest import multipletests
r=PROJECT_ROOT/'runs/20261107_v4_spatialCellChat_bulk_v1';a=r/'bulk/analysis_v1';out=r/'bulk/immune_category_review_v2';out.mkdir(exist_ok=False)
plan=json.loads((r/'config/plan_before_effects.json').read_text());scope=pd.read_csv(r/'bulk/frozen_signatures/all14_signature_eligibility.tsv',sep='\t').set_index('cm');valid={cm:scope.loc[cm,'status']=='FROZEN10GENE_PROGRAM' and cm not in plan['empty_network_CM'] for cm in plan['all_CM']}
mapping={'Wound Healing (Immune C1)':'C1','IFN-gamma Dominant (Immune C2)':'C2','Inflammatory (Immune C3)':'C3','Lymphocyte Depleted (Immune C4)':'C4','Immunologically Quiet (Immune C5)':'C5','TGF-beta Dominant (Immune C6)':'C6'}
source=pd.read_csv(PROJECT_ROOT.parent/'Project_v5/bulk_cohorts_20260904/tcga/TCGA_immune_subtype.txt.gz',sep='\t',dtype=str);assert set(source['Subtype_Immune_Model_Based'].dropna())==set(mapping)
assert all(re.search(r'\(Immune (C[1-6])\)$',key).group(1)==value for key,value in mapping.items())
table=pd.read_csv(a/'all_frozen_immune_slots_private.tsv',sep='\t',float_precision='round_trip');assert len(table)==4620;original=table.copy();counts=[];descriptions=[];errors=[]
for co in plan['bulk']['TCGA_cancers']:
 f=a/('TCGA_'+co+'_scores_clinical_PRIVATE.tsv.gz');key=pd.read_csv(f,sep='\t',nrows=0).columns[0];ss=pd.read_csv(f,sep='\t',index_col=0,dtype={key:str},float_precision='round_trip');canonical=ss.immune_subtype.map(mapping);unknown=ss.immune_subtype.notna()&canonical.isna();assert not unknown.any();category_counts=canonical.value_counts();counts.append({'cohort':co,'n_source_patients':len(ss),'n_with_immune_subtype':int(canonical.notna().sum()),**{'n_'+x:int(category_counts.get(x,0)) for x in ['C1','C2','C3','C4','C5','C6']}})
 for cm in plan['all_CM']:
  d=pd.DataFrame({'score':ss[cm],'group':canonical}).dropna();groups=[g.score.to_numpy() for _,g in d.groupby('group')];which=table.index[table.cohort.eq(co)&table.cm.eq(cm)&table.endpoint.eq('C1_C6_omnibus')];assert len(which)==1;ix=which[0];table.loc[ix,['n','statistic','p']]=[len(d),np.nan,np.nan];table.loc[ix,'n_subtypes']=len(groups)
  if not valid[cm]:table.loc[ix,['status','reason']]=['unavailable','original_empty_network_or_insufficient_frozen_signature']
  elif len(d)>=60 and len(groups)>=2 and min(map(len,groups))>=10 and d.score.nunique()>1:
   v=kruskal(*groups);table.loc[ix,['status','reason','statistic','p']]=['tested','',float(v.statistic),float(v.pvalue)]
   N=len(d);rk=rankdata(d.score);term=0.
   for _,jj in d.groupby('group').indices.items():term+=rk[jj].sum()**2/len(jj)
   _,sizes=np.unique(d.score,return_counts=True);tie=1-np.sum(sizes**3-sizes)/(N**3-N);H=(12*term/(N*(N+1))-3*(N+1))/tie;pv=chi2.sf(H,len(groups)-1);assert abs(H-v.statistic)<1e-10 and abs(pv-v.pvalue)<1e-12;errors.append(max(abs(H-v.statistic),abs(pv-v.pvalue)))
  else:table.loc[ix,['status','reason']]=['unavailable','insufficient_subtype_or_patient_counts_minimum10_each_present_group']
  for group,g in d.groupby('group'):
   descriptions.append({'cohort':co,'cm':cm,'immune_subtype':group,'n_patients':len(g),'median_program':float(g.score.median()),'q25':float(g.score.quantile(.25)),'q75':float(g.score.quantile(.75)),'descriptive_not_extra_inferential_test':True})
# Every non-category statistic and P must be exactly unchanged.
mask=table.endpoint.ne('C1_C6_omnibus');assert np.allclose(table.loc[mask,'statistic'],original.loc[mask,'statistic'],equal_nan=True) and np.allclose(table.loc[mask,'p'],original.loc[mask,'p'],equal_nan=True)
for column,method in [('BH_q','fdr_bh'),('BY_q','fdr_by')]:table[column]=multipletests(table.p.fillna(1),method=method)[1]
# Independent complete-family adjustment by stable reverse cumulative minima.
p=table.p.fillna(1).to_numpy();order=np.argsort(p);m=len(p)
for column,factor in [('BH_q',1),('BY_q',np.sum(1/np.arange(1,m+1)))]:
 q=np.empty(m);q[order]=np.minimum(1,np.minimum.accumulate((p[order]*m*factor/np.arange(1,m+1))[::-1])[::-1]);assert np.allclose(q,table[column],atol=1e-12)
table.to_csv(out/'all4620_immune_slots_reviewed.tsv',sep='\t',index=False);pd.DataFrame(counts).to_csv(out/'all33_source_C1_C6_counts.tsv',sep='\t',index=False);pd.DataFrame(descriptions).to_csv(out/'all14_by33_C1_C6_descriptive_programs.tsv',sep='\t',index=False);(out/'exact_source_category_mapping.json').open('x').write(json.dumps(mapping,indent=2))
cat=table.endpoint.eq('C1_C6_omnibus');fit=table.status.eq('tested');proof={'status':'PASS_EXACT_SOURCE_C1_C6_MAPPING_AND_COMPLETE4620_FAMILY_HARMONIZATION_INDEPENDENT_H_KRUSKAL_BH_BY','UTC':datetime.now(timezone.utc).isoformat(),'original_immune_table_SHA256':hashlib.sha256((a/'all_frozen_immune_slots_private.tsv').read_bytes()).hexdigest(),'n_total_immune_slots':len(table),'n_immune_estimable':int(fit.sum()),'n_category_estimable':int((cat&fit).sum()),'n_BH05':int((fit&table.BH_q.lt(.05)).sum()),'n_BY05':int((fit&table.BY_q.lt(.05)).sum()),'n_category_BH05':int((cat&fit&table.BH_q.lt(.05)).sum()),'max_independent_H_or_p_error':float(max(errors)) if errors else 0.,'original_RNA_statistics_and_P_unchanged':True,'clinical_models_unchanged':True,'no_group_count_gene_or_P_threshold_relaxed':True,'descriptive_category_summaries_not_added_to_inferential_family':True};(out/'independent_audit.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps(proof),flush=True)
