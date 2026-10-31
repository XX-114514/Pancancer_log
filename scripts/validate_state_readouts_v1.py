"""Predeclared score-readout diagnostics; no state selection or module fitting.

Individual sidecars remain on the compute host. Split-half correlations are
within source and working subtype, with equal patient weighting across units.
Unverified patients never enter patient-level correlations. No clinical input.
"""
import argparse,datetime,hashlib,json,os,time,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def corr(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float)
 if len(a)<3:return np.full(a.shape[1],np.nan)
 a=rankdata(a,axis=0);b=rankdata(b,axis=0);a-=a.mean(0);b-=b.mean(0)
 den=np.sqrt((a*a).sum(0)*(b*b).sum(0))
 return np.divide((a*b).sum(0),den,out=np.full(a.shape[1],np.nan),where=den>0)
def seed(*items):return int.from_bytes(hashlib.sha256('|'.join(map(str,items)).encode()).digest()[:8],'little')
def main(R,P,A):
 t=time.monotonic();cfg=json.loads((R/'config.json').read_text())
 assert sha(P/'provenance/terminal_v1.json')==cfg['producer_terminal_SHA256']
 assert sha(A/'outputs/acceptance_v1.json')==cfg['acceptance_SHA256']
 axis=json.loads((P/'outputs/score_program_axis_v1.json').read_text());names=axis['program_ids'];meta=axis['reference_metadata'];idx=json.loads((P/'outputs/source_score_index_v1.json').read_text())
 assert len(names)==126 and len(idx)==37
 records=[];cohort=[];all_units=[];all_means=[];all_errors=[];total=0
 for src in idx:
  s=src['source'];p=P/'outputs'/('private_'+s+'_UCell_v1.npy');assert sha(p)==src['scores_SHA256']
  v=np.load(p,mmap_mode='r');o=pd.read_csv(P/'outputs'/('private_'+s+'_axis_v1.tsv.gz'),sep='\t',dtype={'sample':str,'subtype':str,'patient_id_source_scoped':str})
  assert len(v)==len(o) and v.shape[1]==len(names) and o.source.eq(s).all()
  assert o[['sample','subtype','role','patient_id_source_scoped']].notna().all().all()
  assert o.role.nunique()==1 and o.role.iloc[0] in ['discovery','description_only']
  unit=[];full=[];halfA=[];halfB=[];fullq=[];halfqA=[];halfqB=[];depth=[];ds={n:[] for n in cfg['downsample_sizes']}
  for (sample,subtype),positions in o.groupby(['sample','subtype'],sort=True).indices.items():
   g=o.iloc[positions];x=np.asarray(v[positions],float);n=len(x)
   assert g.patient_id_source_scoped.nunique()==1 and g.patient_verified.nunique()==1
   verified=bool(g.patient_verified.iloc[0]);patient=g.patient_id_source_scoped.iloc[0]
   assert not verified or patient not in ['__unknown__','nan','None','']
   u=dict(source=s,sample=sample,subtype=subtype,role=g.role.iloc[0],patient=patient,patient_verified=verified,n_cells=n)
   unit.append(u);full.append(x.mean(0));fullq.append(np.quantile(x,.95,axis=0));depth.append(float(g.detected_common_genes.mean()))
   ha=[];hb=[];qa=[];qb=[]
   for rep in range(cfg['split_repeats']):
    order=np.random.default_rng(seed(cfg['seed'],s,sample,subtype,'half',rep)).permutation(n);cut=n//2
    ha.append(x[order[:cut]].mean(0) if cut else np.full(len(names),np.nan));hb.append(x[order[cut:]].mean(0) if cut else np.full(len(names),np.nan))
    qa.append(np.quantile(x[order[:cut]],.95,axis=0) if cut else np.full(len(names),np.nan));qb.append(np.quantile(x[order[cut:]],.95,axis=0) if cut else np.full(len(names),np.nan))
   halfA.append(ha);halfB.append(hb);halfqA.append(qa);halfqB.append(qb)
   for size in ds:
    ds[size].append([x[np.random.default_rng(seed(cfg['seed'],s,sample,subtype,'down',size,rep)).choice(n,size,replace=False)].mean(0) for rep in range(cfg['downsample_repeats'])] if n>=size else np.full((cfg['downsample_repeats'],len(names)),np.nan))
  u=pd.DataFrame(unit);f=np.asarray(full);fq=np.asarray(fullq);ha=np.asarray(halfA);hb=np.asarray(halfB);qa=np.asarray(halfqA);qb=np.asarray(halfqB);d=np.asarray(depth)
  # Save raw measurement candidates and sampling uncertainty; all admission masks false.
  err=np.sqrt(np.mean((ha-hb)**2,axis=1))/2
  all_units.extend(unit);all_means.append(f);all_errors.append(err)
  u.to_csv(R/'outputs'/('private_'+s+'_units.tsv.gz'),sep='\t',index=False)
  np.savez_compressed(R/'outputs'/('private_'+s+'_readout_diagnostics.npz'),mean=f,q95=fq,half_mean_A=ha,half_mean_B=hb,half_q95_A=qa,half_q95_B=qb,sampling_se_diagnostic=err,**{'downsample_'+str(k):np.asarray(z) for k,z in ds.items()})
  for subtype,positions in u.groupby('subtype',sort=True).indices.items():
   for cutoff in cfg['cell_thresholds']:
    take=np.asarray([i for i in positions if u.n_cells.iloc[i]>=cutoff and u.patient_verified.iloc[i]],int)
    patient=u.patient.iloc[take].to_numpy();unique=np.unique(patient)
    cohort.append(dict(source=s,subtype=subtype,role=src['role'],min_cells=cutoff,n_units=len(take),n_patients=len(unique)))
    if len(unique)<cfg['min_patients']:continue
    groups=[take[patient==p] for p in unique]
    def pool(z):return np.asarray([z[g].mean(axis=0) for g in groups])
    pf=pool(f);pdpt=pool(d);pa=pool(ha);pb=pool(hb);pqa=pool(qa);pqb=pool(qb)
    meanrho=np.asarray([corr(pa[:,rep],pb[:,rep]) for rep in range(cfg['split_repeats'])]);qrho=np.asarray([corr(pqa[:,rep],pqb[:,rep]) for rep in range(cfg['split_repeats'])]);depthrho=corr(pf,np.repeat(pdpt[:,None],len(names),axis=1))
    # Patient bootstrap, units always stay inside their patient. First predeclared split.
    boots=[];rng=np.random.default_rng(seed(cfg['seed'],s,subtype,cutoff,'bootstrap'))
    for rep in range(cfg['bootstrap_repeats']):
     b=rng.choice(len(groups),len(groups),replace=True);boots.append(corr(pa[b,0],pb[b,0]))
    with np.errstate(all='ignore'):
     ci=np.nanquantile(np.asarray(boots),[.025,.975],axis=0)
    for j,name in enumerate(names):
     records.append(dict(source=s,subtype=subtype,role=src['role'],min_cells=cutoff,n_units=len(take),n_patients=len(unique),program_id=name,split_mean_rho_median=float(np.nanmedian(meanrho[:,j])) if np.isfinite(meanrho[:,j]).any() else np.nan,split_q95_rho_median=float(np.nanmedian(qrho[:,j])) if np.isfinite(qrho[:,j]).any() else np.nan,first_split_rho=float(meanrho[0,j]),first_split_patient_bootstrap_low=float(ci[0,j]),first_split_patient_bootstrap_high=float(ci[1,j]),depth_spearman=float(depthrho[j]),**{'downsample_'+str(size)+'_mean_abs_error':float(np.nanmean(np.abs(np.asarray(ds[size])[take,:,j]-f[take,j,None]))) if np.isfinite(np.asarray(ds[size])[take,:,j]).any() else np.nan for size in ds}))
  total+=len(o);save(R/'checkpoints'/(s+'.json'),dict(completed_UTC=now(),n_complete_cells=total,source=s,n_units=len(u)))
 assert total==776270
 tab=pd.DataFrame(records);tab.to_csv(R/'outputs/private_source_subtype_program_reliability.tsv.gz',sep='\t',index=False);pd.DataFrame(cohort).to_csv(R/'outputs/private_cohort_support.tsv.gz',sep='\t',index=False)
 units=pd.DataFrame(all_units);units.to_csv(R/'outputs/private_material_subtype_registry.tsv.gz',sep='\t',index=False)
 np.save(R/'outputs/private_S_candidate_means.npy',np.concatenate(all_means));np.save(R/'outputs/private_Q_sampling_SE_diagnostic.npy',np.concatenate(all_errors))
 # C is explicitly conditional on analyzed TNK scope, never whole-TME composition.
 c=units.pivot_table(index=['source','sample'],columns='subtype',values='n_cells',aggfunc='sum',fill_value=0);c.to_csv(R/'outputs/private_C_TNK_scope_counts.tsv.gz',sep='\t');c.div(c.sum(axis=1),axis=0).to_csv(R/'outputs/private_C_TNK_scope_proportions.tsv.gz',sep='\t')
 applicability=[dict(program_id=m['program_id'],reference_cell_type=m['reference_cell_type'],role=m['program_role'],module_admitted=False,reason='PENDING_PROVENANCE_PARENT_MATERIAL_TECHNICAL_AND_HELDOUT_REVIEW') for m in meta]
 save(R/'outputs/program_applicability_and_admission.json',applicability)
 agg=[]
 if len(tab):
  for (role,cutoff,name),g in tab.groupby(['role','min_cells','program_id']):
   agg.append(dict(role=role,min_cells=int(cutoff),program_id=name,n_source_subtype_strata=len(g),n_sources=int(g.source.nunique()),median_split_mean_rho=g.split_mean_rho_median.median(),median_split_q95_rho=g.split_q95_rho_median.median(),median_abs_depth_rho=g.depth_spearman.abs().median()))
 pd.DataFrame(agg).to_csv(R/'outputs/aggregate_program_reliability.tsv',sep='\t',index=False)
 save(R/'outputs/input_contract.json',dict(status='DIAGNOSTIC_CANDIDATE_NOT_MODULE_READY',C_scope='conditional_on_analyzed_TNK_cells_not_whole_TME',S='material_unit_by_working_subtype_program_means',M='all_module_admission_false_explicit_reasons_in_program_sidecar',Q='n_cells_and_split_sampling_SE_diagnostic',patient_independence='patient_equal_weight_within_source_subtype_diagnostics_only; material_compatibility_still_pending',split_half='cell_sampling_reliability_not_independent_cohort_validation',downsampling='nested_subsets_of_full_cells_not_independent_validation',no_module_fit=True,no_feature_selection=True,unknown_not_zero=True))
 summary=dict(completed_UTC=now(),n_cells=total,n_sources=len(idx),n_programs=len(names),n_material_subtype_units=len(units),n_verified_patient_keys=int(units.loc[units.patient_verified,['source','patient']].drop_duplicates().shape[0]),n_evaluable_stratum_program_rows=len(tab),n_discovery_rows=int(tab.role.eq('discovery').sum()) if len(tab) else 0,n_description_rows=int(tab.role.eq('description_only').sum()) if len(tab) else 0,all_program_module_admission=False,lockbox_expression_rows_read=0,wall_seconds=time.monotonic()-t,script_SHA256=sha(__file__),config_SHA256=sha(R/'config.json'))
 save(R/'outputs/summary.json',summary);save(R/'provenance/terminal.json',dict(status='READOUT_DIAGNOSTICS_COMPLETE_NOT_BIOLOGICAL_ACCEPTANCE',actual_exit_code=0,**summary))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('output');ap.add_argument('producer');ap.add_argument('audit');a=ap.parse_args();R=Path(a.output)
 try:main(R,Path(a.producer),Path(a.audit))
 except BaseException:save(R/'provenance/failure.json',dict(failed_UTC=now(),traceback=traceback.format_exc()));raise
