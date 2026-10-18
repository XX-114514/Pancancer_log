"""All original native candidates: full/disjoint author-major RNA backgrounds.

Explicitly post-result interpretation following external source-specificity review.
No new tests or endpoint/label changes. Full native controls and gates retained.
"""
import json
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from describe_NSCLC_native_interface_mechanistic_followup_v1 import read_scope
from test_NSCLC_native_interface_association_v1 import inputs
from preflight_NSCLC_continuous_prior_native_association_v1 import make_case_partitions
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
import NSCLC_native_interface_association_math_v1 as kernel
from audit_NSCLC_external_programme_coordination_v3 import calculate


def checked(x,y,q,b):
    p=kernel.prepare(x,y,q,b);v=kernel.effect(p)
    other=calculate(x,y,q,b,minimum=50)
    assert np.array_equal(p['keep'],other['keep']) and abs(v-other['rho'])<1e-10
    return v,len(p['a']),p['design_rank'],p['residual_df']


def main():
    data=load();root=data['root'];budget=resource(root,8*1024**2)
    _,candidates=read_scope(data);_,_,_,z,clinical,_=inputs(root)
    parts=make_case_partitions(data,clinical)
    proof=root/'NSCLC_native_interface_association_results_20261018_v1/private_independent_native_interface_result_audit.json'
    assert sha(proof)=='cdad32bc8755db401dec97b3de3ab9884a6287d08ec0c4637e807cbd2bd78d91'
    out=root/'NSCLC_native_all_candidate_lineage_background_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    units=data['units'];caseids=data['cases'].sampleID.to_numpy(str);ci={v:i for i,v in enumerate(caseids)}
    assert len(set(units.major_cell_type))==3
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        native_result_audit_sha256=sha(proof),native_cache_terminal_sha256=sha(data['native']/'private_terminal.json'),
        scope='ALL four original internal candidate records, all three original author major labels, full-major and sender/receiver-disjoint scopes, all three original case partitions. Explicitly post-result source-specificity diagnostics, no new P/CI.',
        features='Cell-weighted mean per-cell logCP10K ligand from original complete unit cache; log1p background cell count and fraction among all native immune cells. Require20 measured background cells/case before existing non-singleton blocks and full native controls.',
        design='All original8 clinical covariates plus10 sender/receiver QC and full blocks retained, >=50 informative units/df10, original full-span stability. Same-case baseline and adjusted effect both independently OLS checked; no control dropping.',
        limitations='Author T/NK major remains T/NK, not purified T. No new secretion, cell-state, clinical or independent cohort claim. Missing uncaptured lineage not filled as zero ligand; no replacement of original sender.')
    (out/'plan_before_all_native_lineage_background_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    mean=np.load(data['native']/'mean_logCP10K.npy',mmap_mode='r')
    unit_case=units.sampleID.map(ci).fillna(-1).to_numpy(int);incase=unit_case>=0
    nall=np.bincount(unit_case[incase],weights=units.loc[incase,'n_cells'],minlength=len(caseids))
    rows=[];coverage=[]
    for item in candidates.itertuples():
        ligand_index=list(data['genes']).index(item.ligand)
        sv,rv=VIEWS.index(item.sender_view),VIEWS.index(item.receiver_view)
        gene=mean[:,ligand_index];assert np.isfinite(gene).all()
        excluded=set(data['labels'][item.sender_view])|set(data['labels'][item.receiver_view])
        for background_scope in ['full_author_major','sender_receiver_author_pool_disjoint']:
            for major in sorted(units.major_cell_type.unique()):
                take=incase&units.major_cell_type.eq(major).to_numpy()
                if background_scope=='sender_receiver_author_pool_disjoint':take &= ~units.sub_cell_type.isin(excluded).to_numpy()
                n=np.bincount(unit_case[take],weights=units.loc[take,'n_cells'],minlength=len(caseids))
                sums=np.bincount(unit_case[take],weights=units.loc[take,'n_cells'].to_numpy()*gene[take],minlength=len(caseids))
                values=np.divide(sums,n,out=np.full(len(caseids),np.nan),where=n>0)
                extra=np.column_stack([values,np.log1p(n),n/nall])
                # Retain all originally specified partitions; no outcome or expression-based case filtering.
                for partition in ['internal_holdout','discovery','all_frozen_cases']:
                    ix,q,b=parts[(item.sender_view,item.receiver_view,partition)]
                    keep=n[ix]>=20
                    use=ix[keep];qq=q[keep];bb=b[keep]
                    for j in ix:coverage.append(dict(record=item.record,scope=background_scope,major=major,partition=partition,case_index=int(j),background_cells=int(n[j]),eligible=bool(n[j]>=20)))
                    x=z['ligand_mean'][sv,use,item.ligand_axis_index];y=int(item.direction)*z['response'][item.slot_index,use]
                    baseline=adjusted=np.nan;nb=rank=df=na=arank=adf=0
                    try:
                        baseline,nb,rank,df=checked(x,y,qq,bb)
                        adjusted,na,arank,adf=checked(x,y,np.column_stack([qq,extra[use]]),bb)
                        assert nb==na;status='ESTIMABLE'
                    except ValueError as e:status=str(e)
                    rows.append(dict(record=item.record,scope=background_scope,major=major,partition=partition,
                        original_informative_cases=len(ix),background_eligible_cases=len(use),baseline_active_cases=nb,
                        baseline_same_case_rho=baseline,background_adjusted_rho=adjusted,baseline_rank=rank,baseline_df=df,
                        adjusted_rank=arank,adjusted_df=adf,status=status,
                        same_sign_as_baseline=bool(np.isfinite(adjusted) and np.sign(baseline)==np.sign(adjusted))))
    result=pd.DataFrame(rows);assert len(result)==len(candidates)*2*3*3
    result.to_csv(out/'private_all72_native_lineage_background_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(out/'private_all_native_case_background_coverage.tsv.gz',sep='\t',index=False,compression='gzip')
    files=[p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files)<8*1024**2
    terminal=dict(status='COMPLETE_ALL_ORIGINAL_NATIVE_CANDIDATE_LINEAGE_BACKGROUNDS_NO_NEW_P',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),planned_positions=len(result),all_accepted_fits_OLS_crosschecked=True,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},independent_raw_replay=False,main_or_clinical_changes=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
