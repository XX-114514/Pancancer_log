"""All40 endpoint/all60 pair designs; freeze before clinical effects."""
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha, VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import COMMON, matched_quota
from HNSCC_author_state_clinical_kernel_v2 import bootstrap, seed


def main():
    data=load();root=data['root'];source=root/'NSCLC_native_MFA_clinical_coverage_20261016_v1'
    question=root/'NSCLC_native_MFA_clinical_identifiability_20261016_v2/PRIVATE_NEXT_ALL_NATIVE_FACTOR_CLINICAL_QUESTION.json'
    assert sha(question)=='0aad3b4277e60f10a2bd32ebba622971eb87e2bb57b8758d740eb908a4e79c77'
    assert sha(source/'private_terminal.json')=='76bf844a948caab15a05b2eb5a67a3a0e4e1ac04cf92e193ff981566e543f224'
    out=root/'NSCLC_native_MFA_clinical_inputs_20261016_v1';assert not out.exists();quota=matched_quota(root);out.mkdir(mode=0o700)
    kernel=Path(__file__).with_name('HNSCC_author_state_clinical_kernel_v2.py')
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),kernel_sha256=sha(kernel),
        question_sha256=sha(question),clinical_source_terminal_sha256=sha(source/'private_terminal.json'),
        family='All10 frozen native factors x4 own-view endpoints=40; all10 x6 pairs=60 primary conjunction slots.',
        primary='Original internal holdout, complete clinical author cases only; pCR/MPR1 vs non-MPR0. Same cases in both endpoints. Post-treatment association, not baseline prediction or independent patients.',
        outcome='Y=own-view frozen score minus discovery145 mean, divided by discovery145 SD(ddof1); X=author binary response. No clinical re-learning of factor or score scale.',
        controls=COMMON,quality='All original own-view technical variables, unchanged; original B-Plasma structural whole-parent fraction was never a varying control.',
        nuisance='Standardize continuous/numeric controls in each fixed analysis partition; all original columns and intercept retained as Zraw. SVD orthonormal basis for EXACT full nuisance span, relative cutoff1e-12; require retained relative singular>=1e-8, machine-rank agreement; no arbitrary column dropping. Separately identifiable response coefficient only; no individual aliased treatment coefficient.',
        gates='n>=50, response groups>=10, response outside nuisance span by norm>1e-8sqrt(n), df>=10, full reduced design rank, condition<1e8, full leverage<.999999; missing/failed slots retained.',
        primary_inference='19999 restricted-null HC3-residual Rademacher wild studentized HC3 endpoint tests; plus-one P. Each pairP=max(two-sidedPA,PB); whole60 BH/BY, simultaneous99% MonteCarlo-upper BH using40 underlying tests; endpoint40 separate auxiliary family. Require same direction separately, no joint-unique-increment requirement.',
        bootstrap='2000 case draws within ORIGINAL histology x agent blocks, full-case control and discovery score scales fixed; actual rank/condition failures retainedNA. CI only>=1900valid. All leave-one-case estimates descriptive.',
        sensitivities='All40 original internal-holdout mixture-residual scores, discovery140 and complete202 original scores: coefficients only, no additional P/family. Unadjusted rank-biserial only. All raw non-primary outcomes unexamined.',
        calibration='Before true effects: all4 real holdout designs, normal/t3 variance1/heteroskedastic normal errors;1024 independent null simulations each,1999wild each. Every12 rejection rate at alpha.05<=.075. Additional1024 strong coefficient1 normal power perdesign descriptive. If failed, do not compute real effects; no outcome-based rescue.',
        wild_draws=19999,bootstrap_draws=2000,min_valid_bootstrap=1900,calibration_seed=2026101701,
        bootstrap_design_check='All40 planned random streams with zero syntheticY before any real clinical fit; validates design estimability, not outcome.',
        finite_stop='One complete fixed family plus numerical repairs/audit. No post-effect endpoints, score, thresholds, precision or clinical definition variants.',
        quota=quota,new_expression_download=False,clinical_effects_already_computed=False)
    (out/'plan_before_any_native_factor_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    clinical_path=source/'private_all208_clinical_eligibility.tsv'
    source_terminal=json.loads((source/'private_terminal.json').read_text());assert sha(clinical_path)==source_terminal['files'][clinical_path.name]['sha256']
    c=pd.read_csv(clinical_path,sep='\t').set_index('sampleID').loc[data['cases'].sampleID]
    with np.load(root/'NSCLC_native_MFA_interpretation_20261016_v1/private_all_mixture_models_and_scores.npz',allow_pickle=False) as z:
        mixture={k:z[k] for k in z.files if k.startswith('mixture_residual_')}
    cases=data['cases'];tr=data['train'];records=[];view_record={};files=[]
    for partition in ['discovery','internal_holdout','all_frozen_four_view_cases']:
        keep=(np.ones(len(c),bool) if partition=='all_frozen_four_view_cases' else cases.partition.eq(partition).to_numpy())&c.clinical_eligible.to_numpy(bool)
        selected=c.loc[keep];x=selected.response.to_numpy(float);blocks=cases.holdout_block.to_numpy(str)[keep]
        for vi,view in enumerate(VIEWS):
            p=data['prepared'][vi];scores=data['model']['own_view_scores_'+str(vi)]
            mean=scores[tr].mean(axis=0);sd=scores[tr].std(axis=0,ddof=1);assert (sd>1e-12).all()
            y=(scores[keep]-mean)/sd;ym=mixture['mixture_residual_'+str(vi)][keep]/sd
            cov=np.column_stack([selected[COMMON].to_numpy(float),p['quality'][keep]])
            covmean=cov.mean(axis=0);covsd=cov.std(axis=0);raw=np.column_stack([np.ones(len(x)),(cov-covmean)/np.where(covsd>1e-12,covsd,1)])
            u,s,vt=np.linalg.svd(raw,full_matrices=False);rank=int((s>s[0]*1e-12).sum());machine=int((s>s[0]*np.finfo(float).eps*max(raw.shape)).sum())
            assert rank==machine and s[rank-1]/s[0]>=1e-8
            q=u[:,:rank];d=np.column_stack([q,x]);rx=x-q@(q.T@x);assert np.linalg.norm(rx)>1e-8*np.sqrt(len(x))
            assert np.linalg.matrix_rank(d)==d.shape[1] and np.linalg.cond(d)<1e8 and len(x)-len(d.T)>=10
            assert np.max(np.sum(d*np.linalg.pinv(d).T,axis=1))<.999999
            assert len(x)>=50 and min(int((x==0).sum()),int((x==1).sum()))>=10
            assert np.max(abs(raw-q@(q.T@raw)))<1e-10
            path=out/('private_'+partition+'_'+view+'.npz')
            np.savez_compressed(path,D=d,Zraw=raw,X=x,Y=y,Y_mixture=ym,cases=cases.sampleID.to_numpy(str)[keep],strata=blocks,
                source_score_mean=mean,source_score_sd=sd,raw_control_mean=covmean,raw_control_sd=covsd,
                raw_control_names=np.array(['intercept']+COMMON+p['covariate_names'].astype(str).tolist()),nuisance_singular=s)
            files.append(path)
            if partition=='internal_holdout':
                for k in range(10):
                    factor='NSCLC_MFA_20261016_F%02d'%(k+1);identifier='NATIVE_CLINICAL|'+factor+'|'+view
                    b=bootstrap(np.zeros(len(x)),d,blocks,2000,seed(identifier,'bootstrap'))
                    record=dict(record=identifier,factor=factor,view=view,view_index=vi,factor_index=k,n=len(x),response0=int((x==0).sum()),response1=int((x==1).sum()),
                        nuisance_original_columns=raw.shape[1],nuisance_rank=rank,design_columns=d.shape[1],residual_df=len(x)-d.shape[1],condition=float(np.linalg.cond(d)),
                        max_leverage=float(np.max(np.sum(d*np.linalg.pinv(d).T,axis=1))),bootstrap_design_valid=int(np.isfinite(b).sum()),status='READY_NO_CLINICAL_EFFECTS',design_sha256=sha(path))
                    records.append(record);view_record[(factor,view)]=identifier
    pairs=[]
    for k in range(10):
        factor='NSCLC_MFA_20261016_F%02d'%(k+1)
        for a,b in itertools.combinations(VIEWS,2):pairs.append(dict(slot=factor+'|'+a+'|'+b,factor=factor,view_A=a,view_B=b,record_A=view_record[(factor,a)],record_B=view_record[(factor,b)]))
    pd.DataFrame(records).to_csv(out/'private_all40_primary_models.tsv',sep='\t',index=False)
    pd.DataFrame(pairs).to_csv(out/'private_all60_primary_pairs.tsv',sep='\t',index=False)
    report=dict(status='COMPLETE_ALL40_NATIVE_FACTOR_AND60_PAIR_DESIGNS_EQUIVALENT_FULL_NUISANCE_SPANS_NO_CLINICAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},clinical_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
