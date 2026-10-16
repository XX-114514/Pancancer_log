"""All original functional slots with author RNA-composition sensitivity, no new P."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from describe_NSCLC_disjoint_hallmark_v1 import COMMON,QUALITY,correlation,identifiable

QC=['mean_author_pct_mito','mean_author_pct_ribo']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def design(cov,x):
    if not np.isfinite(cov).all() or np.any(cov.std(0)<=1e-12):return None,None,'CONSTANT_OR_NONFINITE_COMPLETE_CONTROL'
    z=np.column_stack([np.ones(len(x)),(cov-cov.mean(0))/cov.std(0)]);rx=rankdata(x);rx=(rx-rx.mean())/rx.std();d=np.column_stack([z,rx])
    if not identifiable(d,d.shape[1]):return None,None,'NONIDENTIFIABLE_COMPLETE_DESIGN'
    return z,np.linalg.qr(z,mode='reduced')[0],'DESCRIBABLE_NO_NEW_P'


def clinical_beta(y,response,z,q):
    if z is None or not np.isfinite(y).all() or np.ptp(y)<=1e-12 or not identifiable(np.column_stack([z,response]),z.shape[1]+1):return np.nan
    rr=response-q@(q.T@response);ys=(y-y.mean())/y.std(ddof=1)
    return float(rr@ys/(rr@rr)) if rr@rr>1e-12 else np.nan


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);qcroot=r/'NSCLC_author_RNA_composition_20261016_v1';qcproof=qcroot/'private_independent_QC_audit.json'
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_author_RNA_composition_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not qcproof.exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(qcproof.read_text())['status'].startswith('PASS_ALL_SOURCE_CELL_ORDER')
    assert sha(qcroot/'private_terminal.json')==json.loads(qcproof.read_text())['producer_terminal_sha256']
    oldroot=r/'NSCLC_disjoint_hallmark_functional_descriptions_20261016_v1';oldproof=oldroot/'private_independent_description_audit.json'
    assert json.loads(oldproof.read_text())['status'].startswith('PASS_ALL10600') and sha(oldroot/'private_terminal.json')==json.loads(oldproof.read_text())['producer_terminal_sha256']
    paths=dict(qc=qcroot/'private_all_native_unit_RNA_composition.tsv',old=oldroot/'private_all10600_functional_descriptions.tsv',scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',readouts=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_all_native_functional_scores_and_QA.npz',references=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/private_reference_axis.tsv')
    hashes={k:sha(p) for k,p in paths.items()}
    assert hashes['qc']==json.loads((qcroot/'private_terminal.json').read_text())['files'][paths['qc'].name]
    assert hashes['old']==json.loads((oldroot/'private_terminal.json').read_text())['files'][paths['old'].name]
    pre=json.loads((r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/plan_before_functional_RNA_measurements.json').read_text())
    assert hashes['scores']==pre['input_sha256']['scores'] and hashes['clinical']==pre['input_sha256']['clinical']
    mt=json.loads((r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_terminal.json').read_text());assert hashes['readouts']==mt['files'][paths['readouts'].name]
    out=r/'NSCLC_RNA_composition_sensitivity_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256=hashes,qc_audit_sha256=sha(qcproof),original_description_audit_sha256=sha(oldproof),
              scope='All212 original signature/author-subtype slots and all10600 fixed function slots. Exactly original clinical cases and >=20 cells; no new endpoint or function selection.',
              design='Original8 clinical + own4 composition/quality + author mean per-cell mitochondrial and ribosomal RNA percentages, all standardized. Intercept plus all14 controls; full16-column programme design, condition<1e8, df>=10. No controls dropped.',
              missing='Any missing QC for an original eligible case makes that original slot unestimable; no case deletion. Preserve old estimates and explicit reasons.',
              outputs='No new P/BH/CI. Programme-versus-each-QC rank correlation under original12 controls; all functional rank correlations with both QC controls, target/background sensitivities, additional response control, descriptive standardized pathology coefficients.',
              interpretation='Composition-sensitive association is not proof of technical artifact: mitochondrial and ribosomal fractions may reflect real biology or mediators. Unchanged correlation is not independent validation.',
              no_formal_cell_state_or_old_family_changes=True)
    (out/'plan_before_any_QC_association.json').write_text(json.dumps(plan,indent=2)+'\n')
    old=pd.read_csv(paths['old'],sep='\t');scores=pd.read_csv(paths['scores'],sep='\t');c=pd.read_csv(paths['clinical'],sep='\t').set_index('sampleID');qc=pd.read_csv(paths['qc'],sep='\t').set_index('unit_index');refs=pd.read_csv(paths['references'],sep='\t').set_index('reference')
    npz=np.load(paths['readouts'],allow_pickle=False);arrays={k:npz[k] for k in ['ULM_score','target_mean','nontarget_mean']};rows=[];slots=[]
    for rec,g in old.groupby('original_clinical_record',sort=True):
        first=g.iloc[0];f=scores[scores.context.eq(first.context)&scores.type.eq(first.type)&scores.program.eq(first.program)&scores.sub_cell_type.eq(first.sub_cell_type)]
        cc=c.reindex(f.sampleID);use=f.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);f=f.loc[use];cc=cc.loc[use];n=len(f);assert n==int(first.eligible_cases) and f.sampleID.is_unique
        extra=qc.reindex(f.unit_index);assert np.array_equal(extra.n_cells.to_numpy(int),f.n_cells.to_numpy(int));x=f.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float);ids=f.index.to_numpy(int)
        cov=np.column_stack([cc[COMMON].to_numpy(float),f[QUALITY].to_numpy(float)]);extended=np.column_stack([cov,extra[QC].to_numpy(float)]);counts=np.bincount(response.astype(int),minlength=2)
        status='DESCRIBABLE_NO_NEW_P';z=None;q=None;z0=None;q0=None;zr=None;qr=None;rank=np.nan;cond=np.nan
        if n<50 or counts.min()<10:status='INSUFFICIENT_ORIGINAL_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif np.ptp(x)<=1e-12:status='CONSTANT_ORIGINAL_PROGRAMME'
        else:
            z0,q0,base=design(cov,x)
            if base!='DESCRIBABLE_NO_NEW_P':status='ORIGINAL_'+base
            elif not np.isfinite(extra[QC].to_numpy(float)).all():status='INCOMPLETE_QC_ORIGINAL_CASES_NOT_REMOVED'
            else:
                z,q,status=design(extended,x)
                if z is not None:
                    rx=rankdata(x);rx=(rx-rx.mean())/rx.std();d=np.column_stack([z,rx]);rank=int(np.linalg.matrix_rank(d));cond=float(np.linalg.cond(d));candidate=np.column_stack([z,response])
                    if identifiable(np.column_stack([candidate,rx]),17):zr=candidate;qr=np.linalg.qr(candidate,mode='reduced')[0]
        slot=dict(original_clinical_record=rec,context=first.context,type=first.type,program=first.program,sub_cell_type=first.sub_cell_type,signature=first.signature,eligible_cases=n,response0=int(counts[0]),response1=int(counts[1]),QC_design_status=status,design_rank=rank,condition=cond,programme_mito_rho=np.nan,programme_ribo_rho=np.nan,programme_MPR_SD_beta_original=np.nan,programme_MPR_SD_beta_QC=np.nan)
        if q0 is not None:
            for label,column in [('mito',QC[0]),('ribo',QC[1])]:slot[f'programme_{label}_rho']=correlation(x,extra[column].to_numpy(float),q0)
            slot['programme_MPR_SD_beta_original']=clinical_beta(x,response,z0,q0)
        if q is not None:slot['programme_MPR_SD_beta_QC']=clinical_beta(x,response,z,q)
        slots.append(slot)
        for row in g.itertuples():
            j=int(refs.loc[row.reference,'reference_column']);y=arrays['ULM_score'][ids,j];fg=arrays['target_mean'][ids,j];bg=arrays['nontarget_mean'][ids,j]
            a=dict(description_id=row.description_id,original_clinical_record=rec,context=row.context,type=row.type,program=row.program,sub_cell_type=row.sub_cell_type,signature=row.signature,reference=row.reference,eligible_cases=n,original_status=row.status,QC_status=status,original_partial_rank_rho=row.partial_rank_rho,QC_partial_rank_rho=np.nan,QC_response_controlled_rho=np.nan,QC_target_mean_rho=np.nan,QC_nontarget_mean_rho=np.nan,original_functional_MPR_SD_beta=row.functional_MPR_minus_nonMPR_SD_beta,QC_functional_MPR_SD_beta=np.nan,QC_minus_original_rho=np.nan)
            if row.status!='DESCRIBABLE_NO_NEW_P':a['QC_status']='ORIGINAL_'+row.status
            elif q is not None:
                a['QC_partial_rank_rho']=correlation(x,y,q)
                if not np.isfinite(a['QC_partial_rank_rho']):a['QC_status']='NO_READOUT_VARIATION_AFTER_COMPLETE_QC_CONTROLS'
                else:
                    a['QC_response_controlled_rho']=correlation(x,y,qr) if qr is not None else np.nan;a['QC_target_mean_rho']=correlation(x,fg,q);a['QC_nontarget_mean_rho']=correlation(x,bg,q);a['QC_functional_MPR_SD_beta']=clinical_beta(y,response,z,q);a['QC_minus_original_rho']=a['QC_partial_rank_rho']-row.partial_rank_rho
            rows.append(a)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),slots=len(slots),function_rows=len(rows)))+'\n')
    table=pd.DataFrame(rows);slot_table=pd.DataFrame(slots);assert len(table)==10600 and len(slot_table)==212
    assert np.array_equal(table.description_id,old.description_id) and not slot_table.original_clinical_record.duplicated().any()
    table.to_csv(out/'private_all10600_functional_QC_sensitivities.tsv',sep='\t',index=False,na_rep='NA');slot_table.to_csv(out/'private_all212_programme_composition_descriptions.tsv',sep='\t',index=False,na_rep='NA')
    summary=[]
    for (signature,reference),f in table.groupby(['signature','reference'],sort=True):
        valid=f.QC_partial_rank_rho.notna();v=f.loc[valid];summary.append(dict(signature=signature,reference=reference,planned_author_subtypes=len(f),QC_describable_author_subtypes=len(v),old_median_same_subtypes=v.original_partial_rank_rho.median(),QC_median_same_subtypes=v.QC_partial_rank_rho.median(),median_change_same_subtypes=v.QC_minus_original_rho.median(),same_case_and_shared_subtypes_not_independent_replications=True))
    pd.DataFrame(summary).to_csv(out/'private_all1200_functional_QC_summary.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='COMPLETE_ALL212_PROGRAMME_AND_ALL10600_FUNCTIONAL_QC_SENSITIVITIES_NO_NEW_P_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),programme_slots=len(slot_table),functional_slots=len(table),status_counts=table.QC_status.value_counts().to_dict(),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
