"""Independently replay every matched-case two-lineage functional description by OLS."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from audit_NSCLC_RNA_composition_sensitivity_v1 import COMMON,QUALITY,QC,rho,full,close


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_matched_function_coordination_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_matched_function_coordination_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_matched_function_effect.json').read_text());assert plan['code_sha256']==receipt['source_sha256']
    for f,h in terminal['files'].items():assert sha(out/f)==h
    paths=dict(scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',QC=r/'NSCLC_author_RNA_composition_20261016_v1/private_all_native_unit_RNA_composition.tsv',readouts=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_all_native_functional_scores_and_QA.npz',scope=r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv',references=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/private_reference_axis.tsv')
    for k,p in paths.items():assert sha(p)==plan['input_sha256'][k]
    proofs=dict(QC=r/'NSCLC_author_RNA_composition_20261016_v1/private_independent_QC_audit.json',RNA_sensitivity=r/'NSCLC_RNA_composition_sensitivity_20261016_v1/private_independent_sensitivity_audit.json',measurement=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_independent_measurement_audit.json',original_coordination=r/'NSCLC_two_end_coordination_results_20261016_v1/private_independent_audit.json')
    for k,p in proofs.items():assert sha(p)==plan['input_audit_sha256'][k] and json.loads(p.read_text())['status'].startswith('PASS_')
    t=pd.read_csv(out/'private_all28600_matched_function_descriptions.tsv',sep='\t');scope=pd.read_csv(paths['scope'],sep='\t');refs=pd.read_csv(paths['references'],sep='\t').sort_values('reference_column');scores=pd.read_csv(paths['scores'],sep='\t').assign(score_row=lambda f:np.arange(len(f)));c=pd.read_csv(paths['clinical'],sep='\t').set_index('sampleID');qc=pd.read_csv(paths['QC'],sep='\t').set_index('unit_index');zfile=np.load(paths['readouts'],allow_pickle=False);arrays={k:zfile[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
    assert len(t)==28600 and len(scope)==572 and not t.duplicated(['record','reference']).any();assert t.record.drop_duplicates().tolist()==scope.record.tolist()
    groups={k:f for k,f in t.groupby('record',sort=False)};maximum=0.;positive=0;negative=0;checked=0;old_designs=0
    for pair in scope.itertuples():
        g=groups[pair.record];assert g.reference.tolist()==refs.reference.tolist()
        ends=[]
        for typ,name in [('typeA',pair.author_subtype_A),('typeB',pair.author_subtype_B)]:
            f=scores[scores.context.eq(pair.context)&scores.program.eq(pair.program)&scores.type.eq(typ)&scores.sub_cell_type.eq(name)].set_index('sampleID');assert f.index.is_unique;ends.append(f)
        a,b=ends;cases=sorted(set(a.index)&set(b.index));a=a.loc[cases];b=b.loc[cases];cc=c.loc[cases];take=a.n_cells.ge(20).to_numpy()&b.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);a=a.loc[take];b=b.loc[take];cc=cc.loc[take];n=len(a)
        assert n==pair.eligible_paired_cases and a.index.tolist()==b.index.tolist()==cc.index.tolist();assert np.array_equal(g.eligible_paired_cases,np.full(50,n))
        for col in ['context','program','author_subtype_A','author_subtype_B','family_A','family_B']:assert g[col].eq(getattr(pair,col)).all()
        assert g.original_pair_status.eq(pair.status).all();qa=qc.reindex(a.unit_index);qb=qc.reindex(b.unit_index);assert np.array_equal(qa.n_cells,a.n_cells) and np.array_equal(qb.n_cells,b.n_cells)
        x=a.raw_weighted_score.to_numpy(float);v=b.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float);ia=a.score_row.to_numpy(int);ib=b.score_row.to_numpy(int);counts=np.bincount(response.astype(int),minlength=2)
        assert g.response0.eq(counts[0]).all() and g.response1.eq(counts[1]).all()
        cov=np.column_stack([cc[COMMON].to_numpy(float),a[QUALITY].to_numpy(float),b[QUALITY].to_numpy(float)]);extended=np.column_stack([cov,qa[QC].to_numpy(float),qb[QC].to_numpy(float)])
        oldp=r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_pair_designs'/(pair.record+'.npz')
        if pair.status=='READY_FOR_FROZEN_TWO_END_TEST':
            old=np.load(oldp,allow_pickle=False);assert np.array_equal(old['cases'],a.index.to_numpy(str)) and np.array_equal(old['unitA'],a.unit_index) and np.array_equal(old['unitB'],b.unit_index)
            close(old['X'],x);close(old['Y'],v);close(np.delete(old['covariates'],8,axis=1),cov);old_designs+=1
        status='DESCRIBABLE_NO_NEW_P';Z=None;z0=None;zr=None;zp=None
        if n<50 or counts.min()<10:status='INSUFFICIENT_ORIGINAL_PAIRED_CASE_OR_RESPONSE_COVERAGE'
        elif not np.isfinite(extended).all() or np.any(extended.std(0)<=1e-12):status='CONSTANT_OR_NONFINITE_COMPLETE_PAIRED_CONTROL'
        else:
            Z=np.column_stack([np.ones(n),(extended-extended.mean(0))/extended.std(0)])
            if not full(Z,21):status='NONIDENTIFIABLE_FULL20_PAIRED_CONTROL_SPACE';Z=None
            else:
                z0=np.column_stack([np.ones(n),(cov-cov.mean(0))/cov.std(0)]);assert full(z0,17);candidate=np.column_stack([Z,response])
                if full(candidate,22):zr=candidate
                rx=rankdata(x);rv=rankdata(v)
                if rx.std()>1e-12 and rv.std()>1e-12:
                    candidate=np.column_stack([Z,(rx-rx.mean())/rx.std(),(rv-rv.mean())/rv.std()])
                    if full(candidate,23):zp=candidate
        for row in g.itertuples():
            j=int(refs.set_index('reference').loc[row.reference,'reference_column']);fa=arrays['ULM_score'][ia,j];fb=arrays['ULM_score'][ib,j];fga=arrays['target_mean'][ia,j];fgb=arrays['target_mean'][ib,j]
            close([row.median_expressed_targets_A,row.median_expressed_targets_B],[np.median(arrays['expressed_target_genes'][ia,j]) if n else np.nan,np.median(arrays['expressed_target_genes'][ib,j]) if n else np.nan])
            expected_status=status;main=np.nan;oldrho=np.nan;response_rho=np.nan;conditional=np.nan;owna=np.nan;ownb=np.nan;crossa=np.nan;crossb=np.nan;fgr=np.nan;bgr=np.nan
            if Z is not None:
                if not np.isfinite(fa).all() or not np.isfinite(fb).all():expected_status='NONFINITE_FUNCTION_ORIGINAL_CASES_NOT_REMOVED'
                elif not np.any(fga>0) or not np.any(fgb>0):expected_status='BACKGROUND_ONLY_ONE_OR_BOTH_FUNCTION_READOUTS'
                elif np.ptp(fa)<=1e-12 or np.ptp(fb)<=1e-12:expected_status='CONSTANT_ONE_OR_BOTH_FUNCTION_READOUTS'
                else:
                    ra=rankdata(fa);ra=(ra-ra.mean())/ra.std()
                    if not full(np.column_stack([Z,ra]),22):expected_status='NONIDENTIFIABLE_COMPLETE_FUNCTION_PAIR_DESIGN'
                    else:
                        main=rho(fa,fb,Z)
                        if not np.isfinite(main):expected_status='NO_FUNCTION_VARIATION_AFTER_COMPLETE_PAIRED_CONTROLS'
                        else:
                            oldrho=rho(fa,fb,z0);response_rho=rho(fa,fb,zr) if zr is not None and full(np.column_stack([zr,ra]),23) else np.nan;conditional=rho(fa,fb,zp) if zp is not None else np.nan;owna=rho(x,fa,Z);ownb=rho(v,fb,Z);crossa=rho(x,fb,Z);crossb=rho(v,fa,Z);fgr=rho(fga,fgb,Z);bgr=rho(arrays['nontarget_mean'][ia,j],arrays['nontarget_mean'][ib,j],Z)
            assert row.status==expected_status
            close([row.function_AB_rho_original16,row.function_AB_rho_QC20,row.function_AB_rho_QC20_plus_response,row.function_AB_rho_QC20_plus_both_programmes,row.Aprogramme_Afunction_rho,row.Bprogramme_Bfunction_rho,row.Aprogramme_Bfunction_rho,row.Bprogramme_Afunction_rho,row.target_mean_AB_rho_QC20,row.nontarget_mean_AB_rho_QC20],[oldrho,main,response_rho,conditional,owna,ownb,crossa,crossb,fgr,bgr])
            if np.isfinite(main):maximum=max(maximum,abs(row.function_AB_rho_QC20-main));positive+=int(main>0);negative+=int(main<0)
            checked+=1
        (out/'private_independent_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),rows_checked=checked))+'\n')
    assert checked==28600 and old_designs==int(scope.status.eq('READY_FOR_FROZEN_TWO_END_TEST').sum()) and positive>0 and negative>0 and t.function_AB_rho_QC20.isna().any()
    proof=dict(status='PASS_ALL28600_MATCHED_FUNCTION_CASE_CONTROL_KEYS_QA_STATUSES_ALL_OLS_PRIMARY_AND_SENSITIVITY_POINTS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),independent_OLS_helper_sha256=sha(Path(__file__).with_name('audit_NSCLC_RNA_composition_sensitivity_v1.py')),producer_terminal_sha256=sha(out/'private_terminal.json'),rows_checked=checked,original_complete_designs_crosschecked=old_designs,maximum_main_rho_difference=maximum,positive_negative_and_unavailable_cases_checked=True,new_P_BH=False)
    p=out/'private_independent_matched_function_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
