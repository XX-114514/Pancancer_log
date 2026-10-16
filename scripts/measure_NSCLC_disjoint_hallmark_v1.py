"""All native disjoint functional scores/QA, streamed from existing full-gene cache."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from disjoint_hallmark_ULM_kernel_v1 import score


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()


def wait(r,path,job):
    p=r/path;rc=r/'HNSCC_receiver_followup_job_control_20261016_v2'/job;j=json.loads((rc/'launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not p.exists():
        q=Path('/proc')/str(j['pid'])/'stat';assert q.exists();s=q.read_text().rsplit(')',1)[1].split();assert s[19]==str(j['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(p.read_text())['status'].startswith('PASS_');return p


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1'
    audit=wait(r,pre.name+'/private_independent_scope_audit.json','audit_NSCLC_disjoint_hallmark_scope_v1')
    calibration=wait(r,'NSCLC_disjoint_hallmark_ULM_calibration_20261016_v1/private_terminal.json','calibrate_disjoint_hallmark_ULM_v1')
    plan=json.loads((pre/'plan_before_functional_RNA_measurements.json').read_text());kernel=Path(__file__).with_name('disjoint_hallmark_ULM_kernel_v1.py')
    assert sha(kernel)==plan['kernel_sha256']==json.loads(calibration.read_text())['kernel_sha256']
    terminal=json.loads((pre/'private_terminal.json').read_text());assert sha(pre/'private_terminal.json')==json.loads(audit.read_text())['producer_terminal_sha256']
    for n,h in terminal['files'].items():assert sha(pre/n)==h
    out=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    record=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),scope_audit_sha256=sha(audit),calibration_sha256=sha(calibration),
                original_plan_sha256=sha(pre/'plan_before_functional_RNA_measurements.json'),kernel_sha256=sha(kernel),measurement='All original source units and eligible masked reference readouts; no clinical association or significance P; raw and zero-expression QA retained.',missingness='Unmeasured numerical scores/means are NaN; unmeasured integer support=-1. Real observed zero support=0, never imputed for a missing readout.')
    (out/'plan_before_real_functional_measurements.json').write_text(json.dumps(record,indent=2)+'\n')
    axis=pd.read_csv(pre/'private_score_row_axis.tsv.gz',sep='\t');scope=pd.read_csv(pre/'private_all1200_functional_scope.tsv',sep='\t');references=pd.read_csv(pre/'private_reference_axis.tsv',sep='\t').reference.tolist();masks=json.loads((pre/'private_physical_column_masks.json').read_text())
    native=r/'NSCLC_full_native_gene_units_20261016_v1';m=np.load(native/'mean_logCP10K.npy',mmap_mode='r');N=len(axis)
    arrays={k:np.full((N,50),np.nan,dtype=np.float64) for k in ['ULM_score','target_mean','nontarget_mean']};arrays['expressed_target_genes']=np.full((N,50),-1,np.int16)
    processed=0
    for (context,program),f in axis.groupby(['context','program'],sort=True):
        assert f.unit_index.is_unique;mask=masks[f'C{context}_{program}'];keep=np.ones(m.shape[1],bool);keep[mask['forbidden_native_columns']]=False;bg=np.flatnonzero(keep)
        mapping=np.full(m.shape[1],-1,int);mapping[bg]=np.arange(len(bg))
        rows=scope[scope.context.eq(context)&scope.program.eq(program)].drop_duplicates('reference').set_index('reference')
        eligible=[name for name in references if rows.loc[name,'status']=='MEASURABLE_INPUT_NO_RNA_SCORE'];columns=[references.index(name) for name in eligible]
        targets=[mapping[np.asarray(mask['targets'][name],int)] for name in eligible];assert all((ix>=0).all() for ix in targets)
        if eligible:
            for start in range(0,len(f),128):
                take=f.iloc[start:start+128];unit=take.unit_index.to_numpy(int);dest=take.score_row.to_numpy(int)
                result=score(np.asarray(m[np.ix_(unit,bg)]),targets)
                for key,a in result.items():arrays[key][np.ix_(dest,columns)]=a
                processed+=len(take)
                (out/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),source_unit_score_rows_processed=processed,total_rows=N,last_context=int(context),last_program=program))+'\n')
        else:processed+=len(f)
    assert processed==N;path=out/'private_all_native_functional_scores_and_QA.npz';np.savez_compressed(path,**arrays)
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())+8*1024**2<plan['output_cap_bytes']
    summaries=[]
    for row in scope.itertuples():
        ids=axis.loc[axis.context.eq(row.context)&axis.type.eq(row.type)&axis.program.eq(row.program),'score_row'].to_numpy(int);j=row.reference_column;y=arrays['ULM_score'][ids,j];v=y[np.isfinite(y)]
        summaries.append(dict(readout_id=row.readout_id,signature=row.signature,reference=row.reference,native_units=len(ids),coverage_status=row.status,
                              finite_scores=len(v),mean_ULM=float(v.mean()) if len(v) else np.nan,SD_ULM=float(v.std()) if len(v) else np.nan,
                              positive_foreground_units=int((arrays['target_mean'][ids,j]>0).sum()),median_expressed_target_genes=float(np.median(arrays['expressed_target_genes'][ids,j])) if len(v) else np.nan))
    pd.DataFrame(summaries).to_csv(out/'private_all1200_measurement_summary.tsv',sep='\t',index=False,na_rep='NA')
    result=dict(status='COMPLETE_ALL24_NATIVE_DISJOINT_FUNCTIONAL_MEASUREMENTS_INDEPENDENT_AUDIT_PENDING_NO_ASSOCIATION_P',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
                source_score_rows=N,functional_slots=1200,finite_score_values=int(np.isfinite(arrays['ULM_score']).sum()),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_clinical_tests=False,new_download=False)
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
