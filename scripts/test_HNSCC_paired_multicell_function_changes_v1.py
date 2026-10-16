"""Complete prespecified patient-paired functional change families; no LR claim."""
import hashlib,itertools,json,math,os,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def exact_tests(d):
    nonzero=d[d!=0];n=len(nonzero)
    if n==0:return 1.,1.,0,0,0.
    positive=int((nonzero>0).sum());tail=min(positive,n-positive);sign=min(1.,2*sum(math.comb(n,k) for k in range(tail+1))/2**n)
    ranks=np.rint(2*rankdata(abs(nonzero),method='average')).astype(int);total=int(ranks.sum());dist=np.zeros(total+1,np.int64);dist[0]=1;end=0
    for rank in ranks:
        prior=dist[:end+1].copy();dist[rank:rank+end+1]+=prior;end+=int(rank)
    assert int(dist.sum())==2**n;observed=int(ranks[nonzero>0].sum());signed=min(1.,2*min(int(dist[:observed+1].sum()),int(dist[observed:].sum()))/2**n)
    return sign,signed,n,positive,observed/2


def calibrate():
    assert exact_tests(np.zeros(4))[:2]==(1.,1.)
    fixture=np.array([1.,1.,-2.,0.,3.,-3.]);d=fixture[fixture!=0];absolute=abs(d)
    ranks=np.array([int((absolute<x).sum())+((absolute==x).sum()+1)/2 for x in absolute]);obs=float(ranks[d>0].sum());all_values=np.array([sum(r for r,s in zip(ranks,bits) if s) for bits in itertools.product([0,1],repeat=len(d))]);expected=min(1.,2*min(float((all_values<=obs).mean()),float((all_values>=obs).mean())))
    assert abs(exact_tests(fixture)[1]-expected)<1e-14 and exact_tests(fixture)[:2]==exact_tests(-fixture)[:2]
    assert exact_tests(np.arange(1.,7.))[:2]==(2/64,2/64)


def main():
    calibrate()
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=r/'HNSCC_native_multicell_functions_20261016_v1';audit=source/'private_independent_multicell_function_audit.json';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_HNSCC_native_multicell_functions_v1';rc=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not audit.exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(audit.read_text())['status'].startswith('PASS_ALL_SELECTED_RAW_CELLS')
    paths={name:source/name for name in ['private_all228_pool_library_axis_and_quality.tsv','private_all16_mask_axis.tsv','private_all50_reference_axis.tsv','private_all_multicell_function_scores_and_QA.npz']};terminal=json.loads((source/'private_terminal.json').read_text())
    for name,p in paths.items():assert sha(p)==terminal['files'][name]
    out=r/'HNSCC_paired_multicell_function_changes_20261016_v1';assert not out.exists()
    fs=subprocess.check_output(['df','-P',str(r.resolve())],text=True).splitlines()[-1].split()[0];quota=subprocess.check_output(['quota','-g','-v','-w','USER002'],text=True);qr=[s.split() for s in quota.splitlines() if s.split() and s.split()[0]==fs];assert len(qr)==1;free=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024;assert free>500_000_000_000+1_073_741_824+268_435_456;out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),measurement_audit_sha256=sha(audit),input_sha256={name:sha(p) for name,p in paths.items()},question='Within reviewed HNSCC patients, which functional RNA patterns change between pre and post in each of the4 frozen pools? This is longitudinal cohort response, not direct communication, causal drug effect or baseline prediction.',scope='All16 original masks x50 functions x4 pools =3200 slots in each of two separate families. No effect-dependent selection. All missing/capture failures retained with internal P=1, displayed unavailable.',pairs='Same reviewed patient, exactly one pre and one ordinary post per pool, both>=20cells, same recorded tissue and treatment arm, finite function and positive foreground RNA support at both visits; recurrence/post2 excluded.>=20 patient pairs and>=10 nonzero changes.',primary='Exact two-sided paired sign test, independent patients. Separate companion exact conditional Wilcoxon signed-rank sign randomization uses doubled average ranks and integer subset-sum distribution, ties handled exactly, zero differences omitted. Wilcoxon additionally assumes sign exchangeability/symmetry; companion is not a rescue test.',numerics='Round post-pre ULM differences to8 decimal places before ranks/signs, fixed for all rows before any effect. Original raw differences also preserved. This numerical rule is not a biological threshold.',multiplicity='Primary3200-slot BH with BY dependence sensitivity; companion3200-slot BH/BY separate. No reuse of old clinical or LR q, no union of favourable families.',interval='2000 whole-patient bootstrap resamples of paired changes, percentile median interval, seed SHA256 of fixed slot; descriptive interval, no simultaneous coverage claim.',quality='Report pre/post target/nontarget RNA, foreground detected genes, pool fractions/counts and all native QC changes. Pairing does not remove time-dependent tissue/processing effects. No claim of abundance/quality-adjusted coordination or outcome prediction.',future='After independent full audit, inspect complete landscape, same-patient multi-pool concordance and QC. Significant one-pool shifts or parallel shifts alone do not prove cell-cell communication.',new_download=False)
    (out/'plan_before_any_paired_function_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    axis=pd.read_csv(paths['private_all228_pool_library_axis_and_quality.tsv'],sep='\t',dtype={'library_id':str,'patient_block':str});masks=pd.read_csv(paths['private_all16_mask_axis.tsv'],sep='\t');refs=pd.read_csv(paths['private_all50_reference_axis.tsv'],sep='\t');arrays=np.load(paths['private_all_multicell_function_scores_and_QA.npz'],allow_pickle=False)
    pools=sorted(axis.pool.unique());patients=sorted(axis.patient_block.unique());pairs={};exclusions=[]
    for pool in pools:
        indices=[]
        for patient in patients:
            g=axis[axis.pool.eq(pool)&axis.patient_block.eq(patient)];pre=g[g.timepoint.eq('pre')];post=g[g.timepoint.eq('post')];reason='READY'
            if len(pre)!=1 or len(post)!=1:reason='NO_UNIQUE_REVIEWED_PRE_POST'
            elif min(int(pre.n_cells.iloc[0]),int(post.n_cells.iloc[0]))<20:reason='CAPTURE_LT20'
            elif pre.tissue.iloc[0]!=post.tissue.iloc[0] or pre.treatment_arm.iloc[0]!=post.treatment_arm.iloc[0]:reason='INCOMPARABLE_RECORDED_TISSUE_OR_ARM'
            if reason=='READY':indices.append((patient,int(pre.pool_index.iloc[0]),int(post.pool_index.iloc[0])))
            exclusions.append(dict(pool=pool,patient_block=patient,status=reason,other_phase_libraries=int((~g.timepoint.isin(['pre','post'])).sum())))
        pairs[pool]=indices
    rows=[];case_rows=[]
    for mask_index,mask in enumerate(masks.itertuples()):
        for reference in refs.itertuples():
            j=int(reference.reference_index)
            for pool in pools:
                slot=f'{mask.mask_id}|{reference.reference}|{pool}';kept=[];candidates=pairs[pool]
                for patient,pre,post in candidates:
                    values=arrays['ULM_score'][mask_index,[pre,post],j];support=arrays['expressed_target_genes'][mask_index,[pre,post],j]
                    if not np.isfinite(values).all() or (support<=0).any():continue
                    raw=float(values[1]-values[0]);rounded=float(np.round(raw,8));kept.append((patient,pre,post,raw,rounded))
                    case_rows.append(dict(mask_id=mask.mask_id,context=mask.context,program=mask.program,reference=reference.reference,pool=pool,patient_block=patient,pre_pool_index=pre,post_pool_index=post,raw_difference=raw,rounded_difference=rounded))
                delta=np.array([x[4] for x in kept]);row=dict(mask_id=mask.mask_id,context=mask.context,program=mask.program,reference=reference.reference,pool=pool,paired_capture_patients=len(candidates),patients=len(kept),RNA_or_nonfinite_exclusions=len(candidates)-len(kept),nonzero_patients=int(np.count_nonzero(delta)),status='INSUFFICIENT_PAIRED_OR_NONZERO_PATIENTS',primary_sign_P=1.,companion_Wilcoxon_P=1.,median_change=np.nan,median_change_low=np.nan,median_change_high=np.nan,positive_fraction=np.nan,positive_nonzero_fraction=np.nan,Wilcoxon_positive_rank=np.nan,new_LR_or_clinical_prediction_claim=False)
                if len(delta)>=20 and np.count_nonzero(delta)>=10:
                    sp,wp,n,pos,rank=exact_tests(delta);seed=int.from_bytes(hashlib.sha256(slot.encode()).digest()[:8],'little');rng=np.random.default_rng(seed);sample=rng.integers(0,len(delta),size=(2000,len(delta)));boot=np.median(delta[sample],axis=1);lo,hi=np.quantile(boot,[.025,.975]);row.update(status='ESTIMATED_LONGITUDINAL_FUNCTIONAL_SHIFT_NOT_COMMUNICATION',primary_sign_P=sp,companion_Wilcoxon_P=wp,median_change=float(np.median(delta)),median_change_low=float(lo),median_change_high=float(hi),positive_fraction=float((delta>0).mean()),positive_nonzero_fraction=pos/n,Wilcoxon_positive_rank=rank)
                if kept:
                    before=np.array([x[1] for x in kept]);after=np.array([x[2] for x in kept])
                    for key in ['target_mean','nontarget_mean','expressed_target_genes']:row['median_change_'+key]=float(np.median(arrays[key][mask_index,after,j]-arrays[key][mask_index,before,j]))
                    for key in ['n_cells','fraction_of_all_source_library_cells','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo']:row['median_change_'+key]=float(np.median(axis.iloc[after][key].to_numpy()-axis.iloc[before][key].to_numpy()))
                rows.append(row)
    table=pd.DataFrame(rows);assert len(table)==3200 and not table.duplicated(['mask_id','reference','pool']).any()
    for prefix,column in [('primary_sign','primary_sign_P'),('companion_Wilcoxon','companion_Wilcoxon_P')]:
        for method in ['fdr_bh','fdr_by']:table[prefix+'_'+method]=multipletests(table[column].to_numpy(float),method=method)[1]
    table.to_csv(out/'private_all3200_paired_functional_change_results.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(case_rows).to_csv(out/'private_all_slot_exact_patient_pair_differences.tsv.gz',sep='\t',index=False);pd.DataFrame(exclusions).to_csv(out/'private_all_patient_pool_capture_exclusions.tsv',sep='\t',index=False);assert sum(p.stat().st_size for p in out.iterdir())<268_435_456
    terminal=dict(status='COMPLETE_ALL3200_PAIRED_FUNCTIONAL_CHANGE_SIGN_AND_COMPANION_WILCOXON_FAMILIES_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),slots=len(table),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},original_clinical_LR_families_unchanged=True,new_direct_communication_or_causal_claim=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
