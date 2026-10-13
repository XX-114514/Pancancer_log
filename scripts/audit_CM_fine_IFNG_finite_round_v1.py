"""All estimated coefficient full-X CR2 replay, full family and scalar score checks.
Same-data numerical verification; not independent biological validation.
"""
import os,json,hashlib,importlib.util
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.sparse import load_npz
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);C=Path(os.environ['PANCANCER_CODE_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);O=R/'immune_same_name_family_source_review_v1/CM_fine_IFNG_finite_round_v1';plan=json.loads((C/'CM_fine_IFNG_finite_round_plan_v1.json').read_text());assert sha(__file__)==plan['auditor_sha256'];inv=json.loads((O/'inventory.json').read_text())
 for name,h in inv['outputs'].items():assert sha(O/name)==h
 spec=importlib.util.spec_from_file_location('explicit',C/'audit_macrophage_T_multisource_response_v1.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a);out=pd.read_csv(O/'all336_fine_IFNG_models.tsv',sep='\t');assert len(out)==336 and not out.duplicated(['pair_id','score','control','program','test_id']).any();assert set(out.pair_id)=={'PAIR015','PAIR034','PAIR038'};assert set(out.test_id)==set(plan['contexts'])|{'POOLED'}
 for _,g in out.groupby(['score','control']):
  assert len(g)==84;order=np.argsort(g.p_exploratory.fillna(1).to_numpy());pv=g.p_exploratory.fillna(1).to_numpy()[order];s=np.minimum(1,np.minimum.accumulate((84*pv/np.arange(1,85))[::-1])[::-1]);q=np.empty(84);q[order]=s;valid=g.status.eq('inference_estimable').to_numpy();np.testing.assert_allclose(g.BH_exploratory.to_numpy()[valid],q[valid],atol=1e-12);np.testing.assert_allclose(g.BY_exploratory.to_numpy()[valid],np.minimum(1,q[valid]*sum(1/np.arange(1,85))),atol=1e-12);assert g.loc[~valid,['p_exploratory','BH_exploratory','BY_exploratory']].isna().all().all()
 checks=0;largest=0;seen=set()
 for f in sorted((O/'private_fixtures').glob('*.npz')):
  d=np.load(f);ix=int(d['row_index']);assert ix not in seen;seen.add(ix);row=out.iloc[ix];e=a.explicit(d['Z'],d['x'],d['Y'],d['blocks']);assert e is not None
  for metric in ['beta','se','df','lower','upper','p_exploratory']:
   value=float(e[metric]) if metric=='df' else float(e[metric][0]);actual=float(row[metric])
   if metric=='p_exploratory' and row.status!='inference_estimable':continue
   np.testing.assert_allclose(actual,value,rtol=2e-7,atol=2e-9,equal_nan=True);largest=max(largest,abs(actual-value));checks+=1
 assert set(out.index[out.beta.notna()])==seen
 dictionary=json.loads((C/'program_reference_coverage_v1.json').read_text());scalar=0
 for f in O.glob('private_score_example_*.npz'):
  sid=f.stem.replace('private_score_example_','');d=np.load(f);root=AR/'20261010_uniform_expression_annotation_v1/by_gse'/sid;genes=pd.read_csv(root/'genes.tsv',sep='\t').gene.tolist();gi={g:i for i,g in enumerate(genes)};values=load_npz(root/'sum_log1p_CP10K.npz')[int(d['unit_row'])].toarray().ravel()/float(d['n_cells']);pct=load_npz(root/'positive_counts.npz')[int(d['unit_row'])].toarray().ravel()/float(d['n_cells']);unique,counts=np.unique(values,return_counts=True);cum=np.cumsum(counts);rank_lookup={float(v):float((end-count+1+end)/2/len(genes)) for v,count,end in zip(unique,counts,cum)}
  for j,pg in enumerate(dictionary['programs']):
   idx=[gi[g] for g in sorted(set(pg['shared_primary_score_genes'])-{'IFNG','IFNGR1','IFNGR2'})];m=sum(float(values[i]) for i in idx)/len(idx);r=sum(rank_lookup[float(values[i])] for i in idx)/len(idx);det=sum(pct[i]>=.01 for i in idx);np.testing.assert_allclose([m,r,det],[d['meanlog'][j],d['rank'][j],d['detected'][j]],atol=1e-12);scalar+=3
 coverage=pd.read_csv(O/'all18_fine_context_coverage.tsv',sep='\t');assert len(coverage)==18 and (coverage.n_fine_both_ge20<=coverage.n_fine_both_captured).all();assert (coverage.n_fine_both_captured<=coverage.n_parent_libraries).all()
 result=dict(status='PASS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),all_model_keys=336,full_family_BH_BY_checked=True,full_X_numeric_values_checked=checks,all_estimated_coefficients_replayed=True,maximum_absolute_error=largest,independent_source_scalar_scores_checked=scalar,not_second_full_program_score_scan=True,not_independent_biological_validation=True,missing_retained_not_negative=True);(O/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
