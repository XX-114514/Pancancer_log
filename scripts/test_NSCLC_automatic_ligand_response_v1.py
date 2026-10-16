"""One complete, automatically enumerated ligand-response family.

Wait for full independent input audit and synthetic statistics calibration.
No candidate, control, mask or threshold changes after viewing effects.
"""
import concurrent.futures,hashlib,json,multiprocessing,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
import NSCLC_automatic_ligand_response_statistics_v1 as statistics

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def wait_terminal(root,relative,job_name,deadline):
    launch=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2'/job_name/'launch_receipt.json').read_text())
    path=root/relative
    while not path.exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();st=p.read_text().rsplit(')',1)[1].split()
        assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z' and time.monotonic()<deadline;time.sleep(15)
    return json.loads(path.read_text())

def worker(record,path,digest,settings):
    with threadpool_limits(limits=1):return statistics.model_statistics(record,path,digest,settings)

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_automatic_ligand_response_preflight_20261016_v1';deadline=time.monotonic()+43200
    audit=wait_terminal(R,Path(I.name)/'private_independent_audit.json','audit_NSCLC_automatic_ligand_response_inputs_v1',deadline)
    assert audit['status'].startswith('PASS_ALL45760_LOGICAL_SLOTS_ALL_UNIQUE_HYPOTHESES_')
    C=R/'NSCLC_ligand_response_statistics_calibration_20261016_v1'
    calibration=wait_terminal(R,Path(C.name)/'private_terminal.json','NSCLC_automatic_ligand_response_statistics_v1',deadline)
    assert calibration['status'].startswith('PASS_EXACT_FINITE_SEQUENTIAL_RULE_')
    assert calibration['statistics_source_sha256']==sha(statistics.__file__)
    scope_path=I/'private_all_unique_ligand_response_designs.tsv';assert sha(scope_path)==audit['unique_scope_sha256'];scope=pd.read_csv(scope_path,sep='\t')
    hashes_path=I/'private_all_ready_design_SHA256.tsv';assert sha(hashes_path)==audit['ready_design_hash_table_sha256']
    hashes=pd.read_csv(hashes_path,sep='\t').set_index('record').sha256
    ready=scope[scope.status.eq('READY_FOR_INDEPENDENT_DESIGN_AUDIT')].copy();assert set(ready.record)==set(hashes.index)
    settings=calibration['real_inference_settings']
    O=R/'NSCLC_automatic_ligand_response_results_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),statistics_code_sha256=sha(statistics.__file__),
        input_audit_sha256=sha(I/'private_independent_audit.json'),statistics_calibration_sha256=sha(C/'private_terminal.json'),unique_scope_sha256=sha(scope_path),
        complete_unique_hypothesis_family=len(scope),ready_designs=len(ready),all_original_logical_slots=45760,no_prior_effect_significance_filter=True,
        settings=settings,workers=4,BLAS_threads_per_worker=1,
        main='Rank sender literal ligand and frozen gene-disjoint receiver score; full18 nuisance plus intercept. Partial residual correlation and HC3 studentized rank-regression t. No pathology outcome control.',
        null='HC2-leverage-corrected null-model residuals, independent Rademacher signs by author case. Wild-model approximation, not exact treatment randomization or independent-patient proof.',
        sequential_rule='Besag-Clifford finite fixed-h rule: stop at FIRST100 absolute null statistics >= absolute observed, p=100/draws; otherwise at999999 draws p=(exceedances+1)/1000000. Ties conservative. No q-threshold or outcome-directed stopping.',
        references=['https://doi.org/10.1093/biomet/78.2.301','https://imsarchives.nus.edu.sg/oldwww/Programs/mcmc/files/besag_tl.pdf','https://academic.oup.com/bioinformatics/article/27/23/3235/236411'],
        family='All unique pair/orientation/ligand hypotheses, including input-unestimable p1 internally and scientific NA. Prior positive/negative slots mapped once. Empty candidate placeholders do not invent a ligand hypothesis.',
        multiplicity='Single complete BH primary screen; BY sensitivity from same P values for arbitrary dependence. No alternate family or denominator selected after results.',
        intervals='2000 case bootstrap draws, rerank X/Y each draw, full19 nuisance rank and20 full design retained. >=1900 valid for pointwise95% CI; not simultaneous or selection-adjusted.',
        diagnostics='All case leave-one-out correlations and required-full-rank failure counts; sender-programme and pathology additional-control points plus unmasked-receiver point descriptive only, no new P.',
        prioritization='For further validation: primary BH<.05, finite CI excluding0, >=95% of all cases have estimable LOO, and>=95% valid LOO signs match main. Still exploratory, source-restricted and not a receptor-specific or causal mechanism.',
        finite_stop='One complete automatic family and audit. No rank21 substitution, score/mask/control changes, relaxed gates or repeated threshold searches if no supported candidates.')
    (O/'plan_before_any_ligand_response_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    results={};loo=[];started=time.monotonic();stream=(O/'private_completed_models.jsonl').open('w',buffering=1)
    tasks=iter(ready.record.tolist());pending={}
    executor=concurrent.futures.ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn'))
    def submit():
        try:record=next(tasks)
        except StopIteration:return False
        path=I/'private_frozen_designs'/(record+'.npz');f=executor.submit(worker,record,str(path),str(hashes.loc[record]),settings);pending[f]=record;return True
    try:
        for _ in range(8):submit()
        while pending:
            done,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:
                record=pending.pop(f);result,ll=f.result();assert result['record']==record;results[record]=result;loo.extend(ll)
                stream.write(json.dumps(result)+'\n');submit()
            (O/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),ready_models=len(ready),completed=len(results),elapsed_seconds=time.monotonic()-started)))
    except BaseException as e:
        for f in pending:f.cancel()
        (O/'private_failure.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),exception_type=type(e).__name__,error=str(e),completed_models=len(results)),indent=2)+'\n')
        raise
    finally:stream.close();executor.shutdown(wait=True,cancel_futures=True)
    rows=[]
    for row in scope.to_dict('records'):
        if row['record'] in results:row.update(results[row['record']])
        rows.append(row)
    table=pd.DataFrame(rows)
    if 'sequential_wild_p' not in table:table['sequential_wild_p']=np.nan
    observed=table.sequential_wild_p.notna();pv=table.sequential_wild_p.fillna(1.).to_numpy(float)
    table['BH_q_complete_unique_family']=np.where(observed,multipletests(pv,method='fdr_bh')[1],np.nan)
    table['BY_q_dependence_sensitivity']=np.where(observed,multipletests(pv,method='fdr_by')[1],np.nan)
    assert table.record.tolist()==scope.record.tolist()
    table.to_csv(O/'private_all_unique_ligand_response_results.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(loo,columns=['record','removed_case_index','rho']).to_csv(O/'private_all_estimable_LOO.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
    report=dict(status='COMPLETE_AUTOMATIC_LIGAND_RESPONSE_FAMILY_INDEPENDENT_NUMERICAL_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        complete_unique_hypothesis_family=len(table),estimable_models=int(observed.sum()),status_counts=table.status.value_counts().to_dict(),
        result_sha256=sha(O/'private_all_unique_ligand_response_results.tsv'),LOO_sha256=sha(O/'private_all_estimable_LOO.tsv.gz'),
        source_selection_was_automatic=True,directed_receptor_specific_or_causal_interaction_validated=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
