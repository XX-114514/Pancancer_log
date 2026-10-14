"""Execute the six predeclared paired/coherence questions, with both slots.
No response labels used; neither causal treatment nor direct LR conclusions.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr


def seed(key):return 20261014+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    prior=root/'HNSCC_frozen_MCP_gene_axis_20261014_v1/prospective_paired_cohort_plan_before_scores.json'
    frozen=json.loads(prior.read_text());assert len(frozen['planned_tests'])==6 and frozen['family'].endswith('p1')
    out=root/'HNSCC_fixed_MCP_paired_results_20261014_v1';assert not out.exists();out.mkdir()
    (out/'execution_contract_before_paired_effects.json').write_text(json.dumps(dict(
      created_UTC=datetime.now(timezone.utc).isoformat(),original_plan_sha256=hashlib.sha256(prior.read_bytes()).hexdigest(),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      paired_shift_statistic='mean patient post-pre difference in full CRC-discovery patient z units',
      permutation='20000 patient sign flips for shifts;20000 patient-label permutations for correlations',
      family='six questions times two frozen slots=12; all unsupported slots internally p1',
      bootstrap='2000 paired-patient blocks;>=1800 valid for CI',min_paired_patients=10,
      no_response_labels=True,not_causal_treatment=True,no_new_gene_or_identity_selection=True),indent=2)+'\n')
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    pre=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    post=scores[scores.timepoint.eq('post')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    paired=pre.index.intersection(post.index).sort_values();pre=pre.loc[paired];post=post.loc[paired]
    assert len(paired)>=10 and pre.notna().all().all() and post.notna().all().all()
    questions=[('paired_shift_T',post.Tnode9.to_numpy()-pre.Tnode9.to_numpy(),None),
      ('paired_shift_Mac',post.Macnode18.to_numpy()-pre.Macnode18.to_numpy(),None),
      ('paired_shift_equal_mean',(post.mean(axis=1)-pre.mean(axis=1)).to_numpy(),None),
      ('pre_cross_type',pre.Tnode9.to_numpy(),pre.Macnode18.to_numpy()),
      ('post_cross_type',post.Tnode9.to_numpy(),post.Macnode18.to_numpy()),
      ('change_cross_type',post.Tnode9.to_numpy()-pre.Tnode9.to_numpy(),post.Macnode18.to_numpy()-pre.Macnode18.to_numpy())]
    rows=[]
    for program in ['MCP1','MCP2']:
        for name,x,y in questions:
            row=dict(program=program,question=name,status='unsupported_original_program_slot',n_paired_patients=len(paired),
                     effect=np.nan,p=np.nan,lower=np.nan,upper=np.nan)
            if program=='MCP1':
                rng=np.random.default_rng(seed(name));boots=[];extreme=0
                if y is None:
                    effect=float(x.mean())
                    for _ in range(20000):extreme+=abs(np.mean(x*rng.choice([-1,1],len(x))))>=abs(effect)-1e-12
                    for _ in range(2000):boots.append(float(x[rng.integers(0,len(x),len(x))].mean()))
                else:
                    effect=float(spearmanr(x,y).statistic)
                    rx=rankdata(x);ry=rankdata(y)
                    rx=(rx-rx.mean())/np.linalg.norm(rx-rx.mean());ry=(ry-ry.mean())/np.linalg.norm(ry-ry.mean())
                    for _ in range(20000):extreme+=abs(rx@ry[rng.permutation(len(x))])>=abs(effect)-1e-12
                    for _ in range(2000):
                        ix=rng.integers(0,len(x),len(x));r=spearmanr(x[ix],y[ix]).statistic
                        if np.isfinite(r):boots.append(float(r))
                low,high=map(float,np.quantile(boots,[.025,.975])) if len(boots)>=1800 else (np.nan,np.nan)
                row.update(status='estimated_partial_measurement_paired_cohort',effect=effect,p=(extreme+1)/20001,
                           lower=low,upper=high,valid_bootstrap=len(boots),invalid_bootstrap=2000-len(boots))
            rows.append(row)
    result=pd.DataFrame(rows);assert len(result)==12
    p=result.p.fillna(1).to_numpy();ix=np.argsort(p);q=np.empty(12)
    q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*12/np.arange(1,13))[::-1])[::-1])
    result['BH_twelve_slot_family']=np.where(result.p.notna(),q,np.nan)
    result.to_csv(out/'all12_private_paired_statistics.tsv',sep='\t',index=False,na_rep='NA')
    (out/'private_terminal.json').write_text(json.dumps(dict(status='ALL12_PAIRED_COHERENCE_QUESTIONS_COMPLETE_PRIVATE',
      finished_UTC=datetime.now(timezone.utc).isoformat(),n_paired_patients=len(paired),
      biopsy_vs_resection_confounding=True,not_causal_or_direct_interaction=True),indent=2)+'\n')
    print('ALL12_PAIRED_COHERENCE_QUESTIONS_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
