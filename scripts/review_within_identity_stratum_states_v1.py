"""Descriptive source-major state review; frozen sample-node units, no inference p."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rho(x,y):
 m=np.isfinite(x)&np.isfinite(y);x=x[m];y=y[m]
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 v=float(spearmanr(x,y).statistic);a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean();q=float(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)));assert np.isclose(v,q,atol=1e-10);return v
def demean(x,g):
 y=np.asarray(x,float).copy()
 for key in np.unique(g):
  m=g==key;v=y[m].copy();y[m]=v-v.mean()
  if np.ptp(v)==0:y[m]=0
 return y
def partial_node_rank(x,y,g):
 return rho_rank_residual(demean(rankdata(x),g),demean(rankdata(y),g))
def rho_rank_residual(x,y):
 if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return np.nan
 return float(np.dot(x,y)/np.sqrt(np.dot(x,x)*np.dot(y,y)))
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);A=V/'all_source_native_UCell_audit_aggregation_20261014_v2';S=V/'all_source_native_UCell_unit_sensitivity_20261014_v1';ai=json.loads((A/'inventory.json').read_text());si=json.loads((S/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL4532933') and si['status']=='COMPLETE_ALL4532933_UNIT_SCORE_SENSITIVITIES_QC_WITH_PRIMARY_AUDIT_AGREEMENT';assert si['audit_inventory_sha256']==sha(A/'inventory.json');D=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(D)==si['dictionary_sha256'];recs={r['programme_id']:r for r in json.loads(D.read_text())['programmes']};F=R/'all_CM_state_native_input_preflight_20261014_v1';fi=json.loads((F/'inventory.json').read_text());U=F/'private_all_frozen_sample_node_units.tsv.gz';assert sha(U)==fi['private_unit_table_sha256'];units=pd.read_csv(U,sep='\t',keep_default_na=False);assert len(units)==35567;rows=[];total=0
 for af in ai['private_aggregate_files']:
  major=af['major'];ap=A/af['relative_file'];sf=next(x for x in si['private_outputs'] if x['major']==major);sp=S/sf['relative_file'];assert sha(ap)==af['sha256'] and sha(sp)==sf['sha256'];a=dict(np.load(ap,allow_pickle=False));s=dict(np.load(sp,allow_pickle=False));assert np.array_equal(a['unit_index'],s['unit_index']) and np.array_equal(a['programme_ids'],s['programme_ids']);md=units.iloc[a['unit_index']].copy();assert np.array_equal(md.unit_index,a['unit_index']) and md.major.eq(major).all();assert np.array_equal(a['cells'],md.n_cells);assert np.allclose(a['programme_mean'],s['scores_maxrank1500_mean'],equal_nan=True,atol=1e-10);total+=int(a['cells'].sum())
  for (sid,nid,stratum),group in md.groupby(['gse_id','node_id','stratum'],sort=True):
   if group.effective_resampling_block.nunique()<10:continue
   allidx=np.flatnonzero((md.gse_id.eq(sid)&md.node_id.eq(nid)&md.stratum.eq(stratum)).to_numpy())
   for j,pid in enumerate(a['programme_ids']):
    r=recs[str(pid)];assert major in r['allowed_majors'];valid=a['valid_score_cells'][allidx,j]>=20;idx=allidx[valid];mm=md.iloc[idx];x=a['programme_mean'][idx,j];assert np.isfinite(x).all();node=mm.node_id.to_numpy(str);rna=np.log1p(s['mean_full_RNA'][idx]);det=np.log1p(s['mean_full_detected_genes'][idx]);support=np.divide(a['informative_gene_support_cells'][idx,j],a['valid_score_cells'][idx,j]);frac=np.divide(a['programme_RNA_fraction_sum'][idx,j],a['valid_score_cells'][idx,j]);record=dict(gse_id=sid,node_id=nid,stratum=stratum,major=major,score_IQR=float(np.percentile(x,75)-np.percentile(x,25)) if len(x) else np.nan,score_SD=float(np.std(x,ddof=1)) if len(x)>1 else np.nan,programme_id=str(pid),method=r['method'],algorithm_stable_candidate=r['algorithm_stable_candidate'],training_source_limited=r['source_limited'],sample_node_units_all=len(allidx),sample_node_units_atleast20_valid_cells=len(idx),distinct_libraries=int(mm.sample_uid.nunique()),effective_analysis_blocks=int(mm.effective_resampling_block.nunique()),reviewed_patient_blocks=int(mm.loc[mm.unit_kind.eq('reviewed_patient'),'effective_resampling_block'].nunique()),original_nodes=int(mm.node_id.nunique()),programme_score_median=float(np.median(x)) if len(x) else np.nan,support_fraction_median=float(np.median(support)) if len(x) else np.nan,programme_RNA_fraction_median=float(np.median(frac)) if len(x) else np.nan,raw_RNA_rho=rho(x,rna),raw_detected_rho=rho(x,det),after_original_node_RNA_partial_rank=partial_node_rank(x,rna,node) if len(x)>=3 else np.nan,after_original_node_detected_partial_rank=partial_node_rank(x,det,node) if len(x)>=3 else np.nan,descriptive_only_no_p=True,formal_states=False)
    resid=demean(x,node) if len(x) else np.array([]);den=float(np.dot(x-x.mean(),x-x.mean())) if len(x) else 0.;record['original_node_R2']=1-float(np.dot(resid,resid))/den if den>0 else np.nan
    for key in ['scores_maxrank1000','scores_maxrank3000','TPM_signature_scores_maxrank1500']:
     yy=s[key+'_mean'][idx,j];ok=np.isfinite(yy)&(s[key+'_valid_cells'][idx,j]>=20);record[key+'_comparable_units']=int(ok.sum());record[key+'_rho']=rho(x[ok],yy[ok]);record[key+'_after_node_partial_rank']=partial_node_rank(x[ok],yy[ok],node[ok]) if ok.sum()>=3 else np.nan
    rows.append(record)
 assert total==4532933;out=pd.DataFrame(rows);assert not out.duplicated(['gse_id','node_id','stratum','major','programme_id']).any();O=V/'within_identity_stratum_state_review_20261015_v1';assert not O.exists();O.mkdir();out.to_csv(O/'within_identity_stratum_programme_review.tsv',sep='\t',index=False);inv=dict(status='COMPLETE_WITHIN_IDENTITY_STRATUM_DESCRIPTIVE_REVIEW_NOT_FORMAL_STATES',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),audit_inventory_sha256=sha(A/'inventory.json'),sensitivity_inventory_sha256=sha(S/'inventory.json'),dictionary_sha256=sha(D),cells=total,sample_node_units=len(units),rows=len(out),no_inference_p_values=True,notes='Source-major records inherit repeated-library/patient units; descriptive correlations only. Node R2 and partial ranks measure association, not causal correction. Coverage includes all eligible measurements, no quantile state allocation.',output_sha256={'within_identity_stratum_programme_review.tsv':sha(O/'within_identity_stratum_programme_review.tsv')},formal_states=False);(O/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n');print(json.dumps(inv))
if __name__=='__main__':main()
