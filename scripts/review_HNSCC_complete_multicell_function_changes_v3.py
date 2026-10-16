"""Complete longitudinal functional landscape and same-patient coordination, private."""
import hashlib,itertools,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rho(x,y):
    return float(spearmanr(x,y).statistic) if len(x)>=10 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan


def partial(x,y,controls):
    if len(x)<20:return np.nan,'INSUFFICIENT_SAME_PATIENTS'
    if controls.shape[1]+2>len(x)-10:return np.nan,'INSUFFICIENT_DF_WITH_COMPLETE_CONTROLS'
    if not np.isfinite(controls).all() or any(np.ptp(c)==0 for c in controls.T):return np.nan,'MISSING_OR_CONSTANT_REQUIRED_CONTROL'
    C=np.column_stack([rankdata(c) for c in controls.T]);C=(C-C.mean(axis=0))/C.std(axis=0);X=np.column_stack([np.ones(len(x)),C]);rx=rankdata(x);ry=rankdata(y);full=np.column_stack([X,rx])
    if np.ptp(x)==0 or np.ptp(y)==0 or np.linalg.matrix_rank(full)!=full.shape[1] or np.linalg.cond(full)>1e8:return np.nan,'CONSTANT_OR_COLLINEAR_COMPLETE_DESIGN'
    ex=rx-X@np.linalg.lstsq(X,rx,rcond=None)[0];ey=ry-X@np.linalg.lstsq(X,ry,rcond=None)[0]
    if np.linalg.norm(ex)<1e-10*np.linalg.norm(rx-rx.mean()) or np.linalg.norm(ey)<1e-10*np.linalg.norm(ry-ry.mean()):return np.nan,'FUNCTION_FULLY_EXPLAINED_BY_REQUIRED_CONTROLS'
    return float(np.corrcoef(ex,ey)[0,1]),'DESCRIPTIVE_COMPLETE_CONTROL_RANK_RESIDUAL_NO_P'


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=r/'HNSCC_paired_multicell_function_changes_20261016_v2';p=source/'private_independent_paired_change_audit.json';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_HNSCC_paired_multicell_function_changes_v2/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not p.exists():
        stat=Path('/proc')/str(rc['pid'])/'stat';assert stat.exists();a=stat.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(p.read_text())['status'].startswith('PASS_ALL3200');t=json.loads((source/'private_terminal.json').read_text())
    for name,h in t['files'].items():assert sha(source/name)==h
    results=pd.read_csv(source/'private_all3200_paired_functional_change_results.tsv',sep='\t');cases=pd.read_csv(source/'private_all_slot_exact_patient_pair_differences.tsv.gz',sep='\t',dtype={'patient_block':str});native=r/'HNSCC_native_multicell_functions_20261016_v1';axis=pd.read_csv(native/'private_all228_pool_library_axis_and_quality.tsv',sep='\t').set_index('pool_index');archive=np.load(native/'private_all_multicell_function_scores_and_QA.npz',allow_pickle=False);npz={k:archive[k] for k in archive.files};masks=pd.read_csv(native/'private_all16_mask_axis.tsv',sep='\t');refs=pd.read_csv(native/'private_all50_reference_axis.tsv',sep='\t');mi={x:i for i,x in enumerate(masks.mask_id)};ri={x:i for i,x in enumerate(refs.reference)}
    out=r/'HNSCC_complete_multicell_function_review_20261016_v3';assert not out.exists();out.mkdir(mode=0o700);rows=[]
    for (reference,pool),g in results.groupby(['reference','pool'],sort=True):
        assert len(g)==16;signed=2*g.positive_fraction-g.nonzero_patients/g.patients
        rows.append(dict(reference=reference,pool=pool,mask_count=16,estimable_masks=int(g.status.str.startswith('ESTIMATED_').sum()),min_paired_patients=int(g.patients.min()),max_paired_patients=int(g.patients.max()),primary_BH_masks=int((g.primary_sign_fdr_bh<.05).sum()),primary_BY_masks=int((g.primary_sign_fdr_by<.05).sum()),companion_BH_masks=int((g.companion_Wilcoxon_fdr_bh<.05).sum()),companion_BY_masks=int((g.companion_Wilcoxon_fdr_by<.05).sum()),all_mask_change_directions_equal=g.median_change.notna().all() and len(set(np.sign(g.median_change)))==1,median_patient_signed_fraction=float(signed.median()),median_raw_ULM_change=float(g.median_change.median()),min_primary_BH=float(g.primary_sign_fdr_bh.min()),max_primary_BH=float(g.primary_sign_fdr_bh.max()),median_foreground_RNA_change=float(g.median_change_target_mean.median()),median_nontarget_RNA_change=float(g.median_change_nontarget_mean.median()),independent_mask_replications=False))
    summary=pd.DataFrame(rows);assert len(summary)==200
    groups={k:g.set_index('patient_block') for k,g in cases.groupby(['mask_id','reference','pool'],sort=False)};pools=sorted(axis.pool.unique());all_pairs=[];quality=['n_cells','fraction_of_all_source_library_cells','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo'];third='extension::B_plus_Plasma'
    def deltas(frame,field):
        pre=axis.loc[frame.pre_pool_index,field].to_numpy(float);post=axis.loc[frame.post_pool_index,field].to_numpy(float)
        return np.log1p(post)-np.log1p(pre) if field=='n_cells' else post-pre
    for mask_id in masks.mask_id:
        for reference in refs.reference:
            for A,B in itertools.combinations(pools,2):
                ga=groups.get((mask_id,reference,A));gb=groups.get((mask_id,reference,B));assert ga is not None and gb is not None;people=sorted(set(ga.index)&set(gb.index));a=ga.loc[people];b=gb.loc[people];x=a.rounded_difference.to_numpy();y=b.rounded_difference.to_numpy();controls=np.column_stack([deltas(f,k) for f in [a,b] for k in quality]);adjusted,status=partial(x,y,controls)
                row=dict(mask_id=mask_id,reference=reference,pool_A=A,pool_B=B,same_patients=len(people),raw_change_rho=rho(x,y),quality12_adjusted_change_rho=adjusted,quality12_status=status,quality12_required_controls=12,third_B_same_patients=0,raw_rho_on_B_subset=np.nan,B_function_only_adjusted_rho=np.nan,B_only_status='B_ALREADY_ENDPOINT',quality12_plus_B_rho=np.nan,quality12_plus_B_status='B_ALREADY_ENDPOINT',new_P_or_clinical_claim=False)
                i=mi[mask_id];j=ri[reference]
                row['baseline_same_patient_rho']=rho(npz['ULM_score'][i,a.pre_pool_index.to_numpy(int),j],npz['ULM_score'][i,b.pre_pool_index.to_numpy(int),j]);row['post_same_patient_rho']=rho(npz['ULM_score'][i,a.post_pool_index.to_numpy(int),j],npz['ULM_score'][i,b.post_pool_index.to_numpy(int),j])
                if third not in [A,B]:
                    bg=groups.get((mask_id,reference,third));subset=sorted(set(people)&set(bg.index));aa=ga.loc[subset];bb=gb.loc[subset];background=bg.loc[subset].rounded_difference.to_numpy();xx=aa.rounded_difference.to_numpy();yy=bb.rounded_difference.to_numpy();qc=np.column_stack([deltas(f,k) for f in [aa,bb] for k in quality]);br,bs=partial(xx,yy,background[:,None]);qr,qs=partial(xx,yy,np.column_stack([qc,background]));row.update(third_B_same_patients=len(subset),raw_rho_on_B_subset=rho(xx,yy),B_function_only_adjusted_rho=br,B_only_status=bs,quality12_plus_B_rho=qr,quality12_plus_B_status=qs)
                all_pairs.append(row)
    pair_table=pd.DataFrame(all_pairs);assert len(pair_table)==4800 and not pair_table.duplicated(['mask_id','reference','pool_A','pool_B']).any()
    summary.to_csv(out/'private_all200_function_pool_summaries.tsv',sep='\t',index=False);pair_table.to_csv(out/'private_all4800_same_patient_change_coordination.tsv',sep='\t',index=False,na_rep='NA')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels={p:p.replace('candidate::','').replace('extension::','').replace('_plus_',' + ') for p in pools};order=refs.reference.tolist();names=[x.replace('HALLMARK_','').replace('_',' ') for x in order];v=summary.pivot(index='pool',columns='reference',values='median_patient_signed_fraction').reindex(index=pools,columns=order);support=summary.pivot(index='pool',columns='reference',values='primary_BH_masks').reindex(index=pools,columns=order)
    fig,axes=plt.subplots(2,1,figsize=(19,8),gridspec_kw={'height_ratios':[1,1]},constrained_layout=True)
    for ax,matrix,cmap,lo,hi,title in [(axes[0],v,'coolwarm',-1,1,'Paired functional direction: median across16 correlated exclusion masks'),(axes[1],support,'Blues',0,16,'Sensitivity support: number of16 masks with primary sign-test BH<0.05')]:
        im=ax.imshow(matrix.to_numpy(float),aspect='auto',cmap=cmap,vmin=lo,vmax=hi);ax.set_yticks(range(4),[labels[p] for p in pools]);ax.set_xticks(range(50),names,rotation=70,ha='right',fontsize=7);ax.set_title(title,fontsize=11);fig.colorbar(im,ax=ax,shrink=.6)
    fig.suptitle('HNSCC reviewed patient pre/post functions: complete50-function landscape',fontsize=14);fig.savefig(out/'private_all50_function_paired_changes.png',dpi=150);plt.close(fig)
    pair_labels=[f'{labels[a]} / {labels[b]}' for a,b in itertools.combinations(pools,2)];display=pair_table.copy();display['pair']=[f'{labels[a]} / {labels[b]}' for a,b in zip(display.pool_A,display.pool_B)];matrix=display.groupby(['pair','reference']).raw_change_rho.median().unstack('reference').reindex(index=pair_labels,columns=order)
    fig,ax=plt.subplots(figsize=(19,5),constrained_layout=True);im=ax.imshow(matrix.to_numpy(float),cmap='coolwarm',vmin=-1,vmax=1,aspect='auto');ax.set_yticks(range(6),pair_labels);ax.set_xticks(range(50),names,rotation=70,ha='right',fontsize=7);ax.set_title('Same-patient pre/post CHANGE correlation: all6 pool pairs, no new P; quality limitations separate');fig.colorbar(im,ax=ax,shrink=.7);fig.savefig(out/'private_all50_function_change_coordination.png',dpi=150);plt.close(fig)
    report='''# Complete HNSCC longitudinal functional review\n\nAll3200 primary and companion paired change slots independently accepted. All200 function/pool summaries and4800 same-patient change correlations retained, not only favourable functions. The16 masks are correlated sensitivity analyses, not16 replications. Primary sign BH/BY and companion Wilcoxon BH/BY are separate families.\n\nSignificance of change within a pool, parallel changes in several pools, and across-patient correlation of the changes answer different questions. The last uses exactly the same patients in both pools; B comparisons use the identical B-captured subset. Complete12 required QC changes are retained; insufficient df, constants or collinearity are explicit, with no automatic covariate deletion. B-only residual correlation is descriptive and does not replace full quality control. Treatment arm, tissue processing, within-pool subtype changes and clinical effect modification remain limitations.\n\nRaw functional patterns are not protein activity, cell-state labels, drug causality, baseline prediction or direct LR communication. Original HNSCC programme clinical/LR families unchanged. Source-selected after other exploration, not an untouched confirmatory validation cohort. Complete tables and figures server private; no scientific export.\n'''
    (out/'PRIVATE_COMPLETE_LONGITUDINAL_FUNCTION_REVIEW.md').write_text(report);terminal=dict(status='COMPLETE_ALL200_FUNCTION_POOL_SUMMARIES_ALL4800_SAME_PATIENT_COORDINATION_AND_COMPLETE_MAPS_REVIEW_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_audit_sha256=sha(p),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
