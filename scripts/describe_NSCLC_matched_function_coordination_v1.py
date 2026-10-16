"""Complete matched-case two-lineage heldout-function descriptions, no new P."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from describe_NSCLC_disjoint_hallmark_v1 import COMMON,QUALITY,correlation,identifiable

QC=['mean_author_pct_mito','mean_author_pct_ribo']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';qcroot=r/'NSCLC_author_RNA_composition_20261016_v1';measurement=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1'
    proofs=dict(QC=qcroot/'private_independent_QC_audit.json',RNA_sensitivity=r/'NSCLC_RNA_composition_sensitivity_20261016_v1/private_independent_sensitivity_audit.json',measurement=measurement/'private_independent_measurement_audit.json',original_coordination=r/'NSCLC_two_end_coordination_results_20261016_v1/private_independent_audit.json')
    for p in proofs.values():assert json.loads(p.read_text())['status'].startswith('PASS_')
    paths=dict(scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',QC=qcroot/'private_all_native_unit_RNA_composition.tsv',readouts=measurement/'private_all_native_functional_scores_and_QA.npz',scope=r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv',references=pre/'private_reference_axis.tsv')
    hashes={k:sha(p) for k,p in paths.items()};original=json.loads((pre/'plan_before_functional_RNA_measurements.json').read_text())
    assert hashes['scores']==original['input_sha256']['scores'] and hashes['clinical']==original['input_sha256']['clinical']
    for key,root in [('QC',qcroot),('readouts',measurement)]:assert hashes[key]==json.loads((root/'private_terminal.json').read_text())['files'][paths[key].name]
    out=r/'NSCLC_matched_function_coordination_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256=hashes,input_audit_sha256={k:sha(p) for k,p in proofs.items()},
              scope='All original572 programme/subtype pairs x all50 original Hallmarks =28600 slots; same exact original matched author cases, no old q or new functional effect selection.',
              main='Rank heldout function in A and B separately and residualize under8 clinical +4 original composition/quality per end + author mito/ribo per end (20 controls+intercept). Main does not condition on pathology outcome; extra-response sensitivity separately.',
              gates='Both endpoints>=20 cells, original complete clinical PD1, >=50 cases and10 per pathology group; no case deletion, all20 controls nonconstant, fullrank21 nuisance and fullrank22 function-A design, df>=10, condition<1e8.',
              descriptors='Function A-B rho under original16 controls, all20 controls, extra pathology control, and extra both ranked original programmes. Programme-own-function and cross-function rhos under all20 controls. Target/non-target RNA and support checks.',
              limits='Same post-treatment source and shared environment; not patient-independent replication, directed communication, protein activity or causal mediation. Gene masks exclude union of both original learned genes, all original weights fixed.',
              no_new_P_BH_or_old_family_changes=True)
    (out/'plan_before_any_matched_function_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    scores=pd.read_csv(paths['scores'],sep='\t');c=pd.read_csv(paths['clinical'],sep='\t').set_index('sampleID');qc=pd.read_csv(paths['QC'],sep='\t').set_index('unit_index');scope=pd.read_csv(paths['scope'],sep='\t');refs=pd.read_csv(paths['references'],sep='\t').sort_values('reference_column')
    assert len(scope)==572 and len(refs)==50;scores=scores.assign(score_row=np.arange(len(scores)))
    zz=np.load(paths['readouts'],allow_pickle=False);arrays={k:zz[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']};rows=[]
    for pair in scope.itertuples():
        endpoints=[]
        for typ,name in [('typeA',pair.author_subtype_A),('typeB',pair.author_subtype_B)]:
            f=scores[scores.context.eq(pair.context)&scores.program.eq(pair.program)&scores.type.eq(typ)&scores.sub_cell_type.eq(name)].set_index('sampleID');assert f.index.is_unique;endpoints.append(f)
        a,b=endpoints;cases=sorted(set(a.index)&set(b.index));a=a.loc[cases];b=b.loc[cases];cc=c.loc[cases];use=a.n_cells.ge(20).to_numpy()&b.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);a=a.loc[use];b=b.loc[use];cc=cc.loc[use]
        n=len(a);assert n==pair.eligible_paired_cases and a.index.tolist()==b.index.tolist()==cc.index.tolist();qa=qc.reindex(a.unit_index);qb=qc.reindex(b.unit_index)
        assert np.array_equal(qa.n_cells.to_numpy(int),a.n_cells.to_numpy(int)) and np.array_equal(qb.n_cells.to_numpy(int),b.n_cells.to_numpy(int))
        x=a.raw_weighted_score.to_numpy(float);v=b.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float);ia=a.score_row.to_numpy(int);ib=b.score_row.to_numpy(int)
        cov=np.column_stack([cc[COMMON].to_numpy(float),a[QUALITY].to_numpy(float),b[QUALITY].to_numpy(float)]);extended=np.column_stack([cov,qa[QC].to_numpy(float),qb[QC].to_numpy(float)])
        counts=np.bincount(response.astype(int),minlength=2);base='DESCRIBABLE_NO_NEW_P';Z=None;q=None;q0=None;zr=None;qp=None
        if n<50 or counts.min()<10:base='INSUFFICIENT_ORIGINAL_PAIRED_CASE_OR_RESPONSE_COVERAGE'
        elif not np.isfinite(extended).all() or np.any(extended.std(0)<=1e-12):base='CONSTANT_OR_NONFINITE_COMPLETE_PAIRED_CONTROL'
        else:
            Z=np.column_stack([np.ones(n),(extended-extended.mean(0))/extended.std(0)])
            if not identifiable(Z,21):base='NONIDENTIFIABLE_FULL20_PAIRED_CONTROL_SPACE'
            else:
                q=np.linalg.qr(Z,mode='reduced')[0];z0=np.column_stack([np.ones(n),(cov-cov.mean(0))/cov.std(0)]);assert identifiable(z0,17);q0=np.linalg.qr(z0,mode='reduced')[0]
                candidate=np.column_stack([Z,response])
                if identifiable(candidate,22):zr=candidate
                rx=rankdata(x);rv=rankdata(v)
                if np.std(rx)>1e-12 and np.std(rv)>1e-12:
                    zp=np.column_stack([Z,(rx-rx.mean())/rx.std(),(rv-rv.mean())/rv.std()])
                    if identifiable(zp,23):qp=np.linalg.qr(zp,mode='reduced')[0]
        for ref in refs.itertuples():
            j=int(ref.reference_column);fa=arrays['ULM_score'][ia,j];fb=arrays['ULM_score'][ib,j];fga=arrays['target_mean'][ia,j];fgb=arrays['target_mean'][ib,j]
            row=dict(record=pair.record,context=pair.context,program=pair.program,author_subtype_A=pair.author_subtype_A,author_subtype_B=pair.author_subtype_B,family_A=pair.family_A,family_B=pair.family_B,reference=ref.reference,eligible_paired_cases=n,response0=int(counts[0]),response1=int(counts[1]),original_pair_status=pair.status,status=base,
                     function_AB_rho_original16=np.nan,function_AB_rho_QC20=np.nan,function_AB_rho_QC20_plus_response=np.nan,function_AB_rho_QC20_plus_both_programmes=np.nan,Aprogramme_Afunction_rho=np.nan,Bprogramme_Bfunction_rho=np.nan,Aprogramme_Bfunction_rho=np.nan,Bprogramme_Afunction_rho=np.nan,target_mean_AB_rho_QC20=np.nan,nontarget_mean_AB_rho_QC20=np.nan,median_expressed_targets_A=float(np.median(arrays['expressed_target_genes'][ia,j])) if n else np.nan,median_expressed_targets_B=float(np.median(arrays['expressed_target_genes'][ib,j])) if n else np.nan)
            if q is not None:
                if not np.isfinite(fa).all() or not np.isfinite(fb).all():row['status']='NONFINITE_FUNCTION_ORIGINAL_CASES_NOT_REMOVED'
                elif not np.any(fga>0) or not np.any(fgb>0):row['status']='BACKGROUND_ONLY_ONE_OR_BOTH_FUNCTION_READOUTS'
                elif np.ptp(fa)<=1e-12 or np.ptp(fb)<=1e-12:row['status']='CONSTANT_ONE_OR_BOTH_FUNCTION_READOUTS'
                else:
                    ra=rankdata(fa);ra=(ra-ra.mean())/ra.std()
                    if not identifiable(np.column_stack([Z,ra]),22):row['status']='NONIDENTIFIABLE_COMPLETE_FUNCTION_PAIR_DESIGN'
                    else:
                        row['function_AB_rho_QC20']=correlation(fa,fb,q)
                        if not np.isfinite(row['function_AB_rho_QC20']):row['status']='NO_FUNCTION_VARIATION_AFTER_COMPLETE_PAIRED_CONTROLS'
                        else:
                            row['function_AB_rho_original16']=correlation(fa,fb,q0);row['function_AB_rho_QC20_plus_response']=correlation(fa,fb,np.linalg.qr(zr,mode='reduced')[0]) if zr is not None and identifiable(np.column_stack([zr,ra]),23) else np.nan
                            row['function_AB_rho_QC20_plus_both_programmes']=correlation(fa,fb,qp) if qp is not None else np.nan
                            row['Aprogramme_Afunction_rho']=correlation(x,fa,q);row['Bprogramme_Bfunction_rho']=correlation(v,fb,q);row['Aprogramme_Bfunction_rho']=correlation(x,fb,q);row['Bprogramme_Afunction_rho']=correlation(v,fa,q)
                            row['target_mean_AB_rho_QC20']=correlation(fga,fgb,q);row['nontarget_mean_AB_rho_QC20']=correlation(arrays['nontarget_mean'][ia,j],arrays['nontarget_mean'][ib,j],q)
            rows.append(row)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),function_rows=len(rows)))+'\n')
    table=pd.DataFrame(rows);assert len(table)==28600 and not table.duplicated(['record','reference']).any();table.to_csv(out/'private_all28600_matched_function_descriptions.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='COMPLETE_ALL572_PAIRS_ALL50_MATCHED_FUNCTION_DESCRIPTIONS_NO_NEW_P_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),rows=len(table),status_counts=table.status.value_counts().to_dict(),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
