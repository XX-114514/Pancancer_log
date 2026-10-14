"""Frozen two-slot heldout patient test; no genes or model fitting on holdout.
Per-patient scores and identities remain private.
"""
import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,h5py
from scipy.stats import spearmanr,pearsonr,rankdata
from anndata._io.specs import read_elem
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'CRC_fine_pair_DIALOGUE_20261014_v1';plan=json.loads(Path(os.environ['PANCANCER_CRC_HOLDOUT_PLAN']).read_text());assert not (O/'heldout_patient_results_v1.tsv').exists();s=pd.read_csv(O/'private_frozen_projection_v1/all_patient_mean_scores_private.tsv',sep='\t');scope=pd.read_csv(O/'gene_supported_program_slots.tsv',sep='\t');p=pd.read_csv(O/'private_inputs/patient_library_partition.tsv',sep='\t');assert len(p)==38 and not p.analysis_block.duplicated().any();src=R/'state_singlecell_inputs_20261014_v1/private_native_sources/GSE178341.h5ad'
 with h5py.File(src) as f:obs=read_elem(f['obs'])
 majors=obs.groupby(['sample_uid','major'],observed=True).size();nodes=obs.groupby(['sample_uid','subCluster'],observed=True).size();labels=['CD8_T_cell_NKG7_9','Macrophage_IFI6_18'];rng=np.random.default_rng(20261014);results=[];leave=[]
 for program in ['MCP1','MCP2']:
  row=dict(program=program,status='no_two_gene_supported_endpoint_scores',n_holdout_patients=0,rho=np.nan,lower=np.nan,upper=np.nan,permutation_p=np.nan,Pearson=np.nan,abundance_partial_rho=np.nan,undefined_bootstrap_draws=0)
  for partition in ['discovery','holdout']:
   g=s[(s.program==program)&(s.partition==partition)];w=g.pivot(index='sample',columns='cell_type',values='score')
   if not {'Tnode9','Macnode18'}.issubset(w.columns):continue
   assert w.notna().all().all();assert set(w.index)==set(p.loc[p.partition==partition,'analysis_block']);x=w.Tnode9.to_numpy();y=w.Macnode18.to_numpy();rho=float(spearmanr(x,y).statistic)
   if partition=='discovery':row['discovery_rho']=rho;continue
   assert len(w)==12;row.update(status='frozen_holdout_estimated',n_holdout_patients=12,rho=rho,Pearson=float(pearsonr(x,y).statistic));rx=rankdata(x);ry=rankdata(y);rx=(rx-rx.mean())/np.linalg.norm(rx-rx.mean());ry=(ry-ry.mean())/np.linalg.norm(ry-ry.mean());extreme=0
   for _ in range(20000):extreme+=abs(rx@ry[rng.permutation(12)])>=abs(rho)-1e-12
   row['permutation_p']=(1+extreme)/20001;boots=[]
   for _ in range(2000):
    ix=rng.integers(0,12,12);r=spearmanr(x[ix],y[ix]).statistic
    if np.isfinite(r):boots.append(float(r))
    else:row['undefined_bootstrap_draws']+=1
   row['lower'],row['upper']=np.quantile(boots,[.025,.975]);lookup=p.set_index('analysis_block');fractions=[]
   for block in w.index:
    uid=lookup.loc[block,'sample_uid'];fractions.append([nodes.loc[uid,labels[0]]/majors.loc[uid,'T cell'],nodes.loc[uid,labels[1]]/majors.loc[uid,'Macrophage']])
   Z=np.column_stack((np.ones(12),np.apply_along_axis(rankdata,0,np.array(fractions))));sx=rankdata(x);sy=rankdata(y)
   if np.linalg.matrix_rank(Z)==3:
    ex=sx-Z@np.linalg.lstsq(Z,sx,rcond=None)[0];ey=sy-Z@np.linalg.lstsq(Z,sy,rcond=None)[0]
    if ex.std()>1e-12 and ey.std()>1e-12:row['abundance_partial_rho']=float(pearsonr(ex,ey).statistic)
   for i in range(12):mask=np.arange(12)!=i;leave.append(dict(program=program,left_out_ordinal=i+1,n_patients=11,rho=float(spearmanr(x[mask],y[mask]).statistic)))
  results.append(row)
 out=pd.DataFrame(results);pv=out.permutation_p.fillna(1).to_numpy();order=np.argsort(pv);qq=np.minimum(1,np.minimum.accumulate((2*pv[order]/np.arange(1,3))[::-1])[::-1]);q=np.empty(2);q[order]=qq;valid=out.status.eq('frozen_holdout_estimated');out['BH_two_slot_family']=np.where(valid,q,np.nan);out.to_csv(O/'heldout_patient_results_v1.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(leave).to_csv(O/'heldout_leave_one_patient_v1.tsv',sep='\t',index=False)
 j=dict(status='HELDOUT_PATIENT_EVALUATION_DONE_PENDING_INDEPENDENT_AUDIT',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(Path(os.environ['PANCANCER_CRC_HOLDOUT_PLAN'])),clinical_outcomes_used=False,project_wide_independent_validation=False,family_slots=2,results=json.loads(out.to_json(orient='records')));(O/'heldout_evaluation_inventory_v1.json').write_text(json.dumps(j,indent=2)+'\n');print(out.to_json(orient='records'))
if __name__=='__main__':main()
