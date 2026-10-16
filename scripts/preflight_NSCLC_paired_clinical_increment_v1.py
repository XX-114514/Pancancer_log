"""Freeze all 572 matched-case designs; no clinical association effects."""
import hashlib,json,os,subprocess
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);old=r/'NSCLC_two_end_coordination_preflight_20261016_v1'
    for sub in ['NSCLC_frozen_signature_clinical_results_20261016_v1','NSCLC_two_end_coordination_results_20261016_v1']:
        assert json.loads((r/sub/'private_independent_audit.json').read_text())['status'].startswith('PASS_')
    scope_path=old/'private_all572_planned_pair_designs.tsv';oldplan=json.loads((old/'plan_before_any_two_end_effect.json').read_text())
    paths={'scores':r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',
           'case_eligibility':r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',
           'parameters':r/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'}
    for k,p in paths.items():assert sha(p)==oldplan['input_sha256'][k]
    kernel=Path(__file__).with_name('paired_clinical_increment_kernel_v1.py')
    out=r/'NSCLC_paired_programme_clinical_increment_preflight_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    df=subprocess.run(['df','-P',str(r.resolve())],capture_output=True,text=True,timeout=20)
    quota=subprocess.run(['quota','-g','-v','-w','USER002'],capture_output=True,text=True,timeout=20)
    assert df.returncode==quota.returncode==0;fs=df.stdout.splitlines()[-1].split()[0]
    qr=[x.split() for x in quota.stdout.splitlines() if x.split() and x.split()[0]==fs];assert len(qr)==1
    remaining=(int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024
    assert remaining>500_000_000_000+1_073_741_824+268_435_456,'Preserve original reserve and output cap'
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),
              original_scope_sha256=sha(scope_path),input_sha256={k:sha(p) for k,p in paths.items()},
              family=572,clinical_controls=COMMON,unit_controls=UNIT,columns=['intercept']+COMMON+['A_'+k for k in UNIT]+['B_'+k for k in UNIT]+['programme_A','programme_B'],
              question='Do BOTH original programmes add conditional clinical information on the SAME author cases beyond the other programme and all composition/technical/available clinical controls?',
              outcome='Previously frozen named PD1 and complete clinical author cases; author pCR/MPR=1 vs non-MPR=0. Post-treatment RNA, not baseline prediction.',
              primary='Full-rank linear conditional association. Separate restricted-null HC2 Rademacher wild HC3 tests of last-two coefficients; pairP=max(PA,PB). Complete572 BH; BY and simultaneous99% MCupper sensitivity.',
              gates='n>=50, each response group>=10, both subtype cells>=20, all19 columns fullrank/condition<1e8/df>=10/leverage<.999999; no covariate deletion.',
              bootstrap='2000 paired-case row draws on original full-case standardization, require1900 valid; all leave-one-case coefficients descriptive.',
              descriptive_only='Matched-single-end coefficients with same16 controls; partial incremental R2; shared-gene-removed conditional effects if both retained80% L1/95% measurable/20 genes.',
              clinical_independent_patient_identity_accepted=False,communication_or_causal_synergy_test=False,
              wild_draws=199999,bootstrap_draws=2000,min_valid_bootstrap=1900,seed=2026101611,
              quota_free_bytes=remaining,quota_reserve_decimal_bytes=500_000_000_000,output_cap_bytes=268_435_456,
              finite_stop='No post-effect endpoint, outcome, score, threshold or scope replacement. Full independent input/calibration before true effects.')
    (out/'plan_before_any_paired_clinical_increment_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    s=pd.read_csv(paths['scores'],sep='\t');c=pd.read_csv(paths['case_eligibility'],sep='\t').set_index('sampleID')
    assert c.index.is_unique and pd.api.types.is_bool_dtype(c.clinical_eligible)
    oldscope=pd.read_csv(scope_path,sep='\t');assert len(oldscope)==572 and oldscope.record.is_unique
    # Original nonshared readouts are reconstructed from frozen native genes, not clinical outcomes.
    params=pd.read_csv(paths['parameters'],sep='\t');nroot=r/'NSCLC_full_native_gene_units_20261016_v1'
    genes=pd.Index(pd.read_csv(nroot/'private_native_gene_axis.tsv',sep='\t').gene);m=np.load(nroot/'mean_logCP10K.npy',mmap_mode='r');masked={}
    for context,program in oldscope[['context','program']].drop_duplicates().itertuples(index=False):
        pp={t:params[params.context.eq(context)&params.program.eq(program)&params.type.eq(t)] for t in ['typeA','typeB']}
        shared=set(pp['typeA'].gene)&set(pp['typeB'].gene)
        for typ,p in pp.items():
            keep=~p.gene.isin(shared);available=genes.get_indexer(p.gene)>=0;w=np.abs(p.weight.to_numpy())
            if not keep.any() or w[keep].sum()/w.sum()<.8 or w[keep&available].sum()/w[keep].sum()<.95 or (keep&available).sum()<20:continue
            p=p.loc[keep&available];f=s[s.context.eq(context)&s.program.eq(program)&s.type.eq(typ)]
            values=np.asarray(m[np.ix_(f.unit_index.to_numpy(int),genes.get_indexer(p.gene))])@(p.weight/p.gene_sd).to_numpy()-float((p.weight*p.gene_mean/p.gene_sd).sum())
            masked[(context,program,typ)]=pd.Series(values,index=f.unit_index)
    designs=out/'private_designs';designs.mkdir(mode=0o700);records=[];hashes=[]
    for row in oldscope.itertuples():
        fa=s[s.context.eq(row.context)&s.program.eq(row.program)&s.type.eq('typeA')&s.sub_cell_type.eq(row.author_subtype_A)].set_index('sampleID')
        fb=s[s.context.eq(row.context)&s.program.eq(row.program)&s.type.eq('typeB')&s.sub_cell_type.eq(row.author_subtype_B)].set_index('sampleID')
        common=sorted(set(fa.index)&set(fb.index));a=fa.loc[common];b=fb.loc[common];cc=c.loc[common]
        take=a.n_cells.ge(20).to_numpy()&b.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool)
        a=a.loc[take];b=b.loc[take];cc=cc.loc[take];assert a.index.tolist()==b.index.tolist()==cc.index.tolist() and a.index.is_unique
        x=a.raw_weighted_score.to_numpy(float);v=b.raw_weighted_score.to_numpy(float);y=cc.response.to_numpy(float);n=len(y)
        cov=np.column_stack([cc[COMMON].to_numpy(float),a[UNIT].to_numpy(float),b[UNIT].to_numpy(float)]);raw=np.column_stack([cov,x,v])
        sd=raw.std(0) if n else np.zeros(18);means=raw.mean(0) if n else np.zeros(18)
        d=np.column_stack([np.ones(n),(raw-means)/np.where(sd>1e-12,sd,1)])
        sv=np.linalg.svd(d,compute_uv=False) if n else np.array([])
        rank=int(np.linalg.matrix_rank(d)) if n else 0;relative_rank=int((sv>sv[0]*1e-12).sum()) if n else 0
        condition=float(np.linalg.cond(d)) if n else np.inf
        groups=np.bincount(y.astype(int),minlength=2);status='READY_NO_CLINICAL_EFFECT_COMPUTED'
        if n<50 or groups.min()<10:status='INSUFFICIENT_MATCHED_CASE_OR_RESPONSE_GROUP'
        elif not np.isfinite(d).all() or (sd<=1e-12).any() or rank!=19:status='NONFINITE_CONSTANT_OR_NONFULL_DESIGN'
        elif relative_rank!=rank:status='NUMERICALLY_AMBIGUOUS_RANK'
        elif n-rank<10 or condition>=1e8:status='LOW_DF_OR_ILL_CONDITIONED'
        elif np.max(np.sum(d*np.linalg.pinv(d).T,axis=1))>=.999999:status='EXTREME_LEVERAGE'
        rec={k:getattr(row,k) for k in ['record','context','program','author_subtype_A','author_subtype_B','family_A','family_B']}
        rec.update(eligible_cases=n,response0=int(groups[0]),response1=int(groups[1]),design_rank=rank,design_columns=19,residual_df=n-rank,condition=condition,status=status)
        if status=='READY_NO_CLINICAL_EFFECT_COMPUTED':
            ma=masked.get((row.context,row.program,'typeA'));mb=masked.get((row.context,row.program,'typeB'));ok=ma is not None and mb is not None
            path=designs/(row.record+'.npz')
            np.savez_compressed(path,D=d,Y=y,raw_A=x,raw_B=v,C_raw=cov,raw_means=means,raw_sd=sd,cases=a.index.to_numpy(str),unitA=a.unit_index.to_numpy(int),unitB=b.unit_index.to_numpy(int),
                                masked_A=ma.reindex(a.unit_index).to_numpy() if ok else np.full(n,np.nan),masked_B=mb.reindex(b.unit_index).to_numpy() if ok else np.full(n,np.nan))
            hashes.append(dict(record=row.record,sha256=sha(path)));rec['shared_gene_mask_available']=ok
        records.append(rec)
    t=pd.DataFrame(records);assert len(t)==572 and t.record.tolist()==oldscope.record.tolist()
    t.to_csv(out/'private_all572_scope.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(hashes).to_csv(out/'private_ready_design_SHA256.tsv',sep='\t',index=False)
    terminal=dict(status='COMPLETE_ALL572_MATCHED_CASE_FULL_CONTROL_DESIGNS_NO_CLINICAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),planned_models=len(t),status_counts=t.status.value_counts().to_dict(),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
