"""Freeze every compatible two-end programme/subtype combination before effects."""
import hashlib,json,os,itertools
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary','response']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'NSCLC_context_native_programme_projection_20261016_v1'
    assert json.loads((P/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL56_SCOPE_')
    C=R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1'
    assert json.loads((R/'NSCLC_frozen_signature_clinical_results_20261016_v1/private_independent_audit.json').read_text())['status'].startswith('PASS_ALL_PLANNED_CLINICAL_SCOPE_')
    O=R/'NSCLC_two_end_coordination_preflight_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    scores_path=P/'private_native_author_subtype_programme_scores.tsv.gz';scores=pd.read_csv(scores_path,sep='\t')
    params_path=R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    params=pd.read_csv(params_path,sep='\t')
    scope=pd.read_csv(P/'private_all56_projection_scope.tsv',sep='\t')
    active=scope[scope.status.str.startswith('PARTIAL_FROZEN')]
    parents=active.groupby(['context','program']).filter(lambda x:len(x)==2)
    assert parents.context.nunique()==4 and len(parents)==16
    c=pd.read_csv(C/'private_all243_clinical_eligibility.tsv',sep='\t').set_index('sampleID',drop=False)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        input_sha256=dict(scores=sha(scores_path),parameters=sha(params_path),case_eligibility=sha(C/'private_all243_clinical_eligibility.tsv')),
        scope='All4 compatible parent contexts x2 MCPs x all native author subtype cartesian products, no filter by previous clinical/LR significance.',
        question='Do fixed programmes in distinct native cell lineages co-vary within the same cases beyond pathology group, treatment, abundance and quality? Coordination is not directed communication.',
        units='One author case per pair, both native subtype cell counts>=20, previous complete clinical/recorded-PD1 cohort; retain all absent/insufficient pairs.',
        nuisance=COMMON+['A_'+x for x in UNIT]+['B_'+x for x in UNIT],
        effect='Rank both raw fixed-weight signatures; residualize ranks against all17 standardized raw nuisance variables plus intercept; report residual correlation. Pathology response retained as a control.',
        test='Studentized HC3 t for ranked-B regressed on ranked-A and full nuisance. 20000 null HC2 Rademacher case-block wild draws; not exact treatment randomization.',
        bootstrap='2000 case resamples with X/Y ranks recomputed in each draw, all nuisance kept; require full nuisance rank and1900 finite rho draws for CI.',
        gates='>=50 paired cases, >=10 each pathological response group, full nuisance rank18 and full rank19, residualdf>=10, condition<1e8, leverage<.999999.',
        family='Every planned pair is in one BH family; unestimable p1 in denominator but scientific result NA.',
        sensitivity='No-pathology-control partial rank, fully controlled linear signatures and old-residual rank, plus shared-gene-removed scores where both retain>=80% original absolute weights and>=95% remaining weights measured; point descriptions only.',
        masking='For shared-gene sensitivity remove all intersecting frozen A/B genes from both. Never reweight or substitute genes; if mask changes programme too much mark unavailable.',
        interpretation='New exploratory source/author-subtype coordination, not original fine-identity acceptance, independent validation of post-selected clinical results, clinical prediction or causality.',
        finite_stop='One full family and required audit; no subtype/covariate/score/rank changes after effects.',seed=2026101605,wild_draws=20000,bootstrap_draws=2000,min_valid_bootstraps=1900)
    (O/'plan_before_any_two_end_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    N=R/'NSCLC_full_native_gene_units_20261016_v1';genes=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene)
    M=np.load(N/'mean_logCP10K.npy',mmap_mode='r')
    masked={};mask_scope=[]
    for (context,program),slots in parents.groupby(['context','program'],sort=True):
        pp={typ:params[params.context.eq(context)&params.type.eq(typ)&params.program.eq(program)] for typ in ['typeA','typeB']}
        shared=set(pp['typeA'].gene)&set(pp['typeB'].gene)
        for typ,p in pp.items():
            available=genes.get_indexer(p.gene)>=0;keep=~p.gene.isin(shared).to_numpy();w=p.weight.to_numpy()
            remaining=float(np.abs(w[keep]).sum()/np.abs(w).sum());native=float(np.abs(w[keep&available]).sum()/np.abs(w[keep]).sum()) if keep.any() else 0.
            ok=remaining>=.8 and native>=.95 and int((keep&available).sum())>=20
            f=scores[scores.context.eq(context)&scores.type.eq(typ)&scores.program.eq(program)]
            if ok:
                chosen=p.loc[keep&available];columns=genes.get_indexer(chosen.gene);ix=f.unit_index.to_numpy(int)
                z=(np.asarray(M[np.ix_(ix,columns)])-chosen.gene_mean.to_numpy())/chosen.gene_sd.to_numpy()
                masked[(int(context),program,typ)]=pd.Series(z@chosen.weight.to_numpy(),index=f.unit_index)
            mask_scope.append(dict(context=context,program=program,type=typ,shared_frozen_genes=len(shared),original_genes=len(p),nonshared_original_L1_fraction=remaining,
                nonshared_native_weight_fraction=native,nonshared_score_available=ok))
    pd.DataFrame(mask_scope).to_csv(O/'private_all16_parent_shared_gene_masks.tsv',sep='\t',index=False)
    design=O/'private_pair_designs';design.mkdir(mode=0o700);records=[]
    for (context,program),slots in parents.groupby(['context','program'],sort=True):
        a=scores[scores.context.eq(context)&scores.program.eq(program)&scores.type.eq('typeA')]
        b=scores[scores.context.eq(context)&scores.program.eq(program)&scores.type.eq('typeB')]
        sa=sorted(a.sub_cell_type.unique());sb=sorted(b.sub_cell_type.unique())
        for nameA,nameB in itertools.product(sa,sb):
            fa=a[a.sub_cell_type.eq(nameA)].set_index('sampleID');fb=b[b.sub_cell_type.eq(nameB)].set_index('sampleID')
            common=sorted(set(fa.index)&set(fb.index));aa=fa.loc[common];bb=fb.loc[common];cc=c.loc[common]
            use=aa.n_cells.ge(20).to_numpy()&bb.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool)
            aa=aa.loc[use];bb=bb.loc[use];cc=cc.loc[use];n=len(aa)
            assert aa.index.tolist()==bb.index.tolist()==cc.index.tolist() and aa.index.is_unique
            assert not set(aa.unit_index)&set(bb.unit_index)
            cov=np.column_stack([cc[COMMON].to_numpy(float),aa[UNIT].to_numpy(float),bb[UNIT].to_numpy(float)])
            sd=cov.std(0) if n else np.zeros(17)
            Z=np.column_stack([np.ones(n),(cov-cov.mean(0))/np.where(sd>1e-12,sd,1)]) if n else np.empty((0,18))
            x=aa.raw_weighted_score.to_numpy(float);y=bb.raw_weighted_score.to_numpy(float)
            from scipy.stats import rankdata
            xr=rankdata(x);xr=(xr-xr.mean())/xr.std() if n and xr.std()>0 else np.zeros(n)
            D=np.column_stack([Z,xr]);rank=int(np.linalg.matrix_rank(D)) if n else 0;rankZ=int(np.linalg.matrix_rank(Z)) if n else 0
            condition=float(np.linalg.cond(D)) if n else np.inf;groups=np.bincount(cc.response.to_numpy(int),minlength=2)
            status='READY_FOR_FROZEN_TWO_END_TEST'
            if n<50 or groups.min()<10:status='INSUFFICIENT_PAIRED_CASE_OR_RESPONSE_COVERAGE'
            elif (sd<=1e-12).any() or rankZ!=18 or rank!=19:status='CONSTANT_OR_COLLINEAR_FULL_NUISANCE_OR_SIGNATURE'
            elif n-rank<10 or condition>=1e8:status='LOW_DF_OR_ILL_CONDITIONED'
            elif np.ptp(y)==0:status='NO_TARGET_SIGNATURE_VARIATION'
            elif np.max(np.sum(D*np.linalg.pinv(D).T,axis=1))>=.999999:status='EXTREME_LEVERAGE'
            record='PAIR'+str(len(records)+1).zfill(4)
            row=dict(record=record,context=context,program=program,author_subtype_A=nameA,author_subtype_B=nameB,
                family_A=slots[slots.type.eq('typeA')].source_family.iloc[0],family_B=slots[slots.type.eq('typeB')].source_family.iloc[0],
                shared_cases_any=len(common),eligible_paired_cases=n,response0=int(groups[0]),response1=int(groups[1]),nuisance_rank=rankZ,design_rank=rank,
                design_columns=19,residual_df=n-rank,condition=condition,status=status)
            if status=='READY_FOR_FROZEN_TWO_END_TEST':
                am=masked.get((int(context),program,'typeA'));bm=masked.get((int(context),program,'typeB'))
                masked_ok=am is not None and bm is not None
                np.savez_compressed(design/(record+'.npz'),X=x,Y=y,Z=Z,covariates=cov,oldA=aa.score.to_numpy(float),oldB=bb.score.to_numpy(float),
                    maskedA=am.reindex(aa.unit_index).to_numpy() if masked_ok else np.full(n,np.nan),maskedB=bm.reindex(bb.unit_index).to_numpy() if masked_ok else np.full(n,np.nan),
                    cases=aa.index.to_numpy(str),unitA=aa.unit_index.to_numpy(int),unitB=bb.unit_index.to_numpy(int),response=cc.response.to_numpy(int))
                row['shared_gene_mask_available']=masked_ok
            records.append(row)
    t=pd.DataFrame(records);assert len(t)==572 and not t.duplicated(['context','program','author_subtype_A','author_subtype_B']).any()
    t.to_csv(O/'private_all572_planned_pair_designs.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL572_TWO_END_AUTHOR_SUBTYPE_DESIGNS_FROZEN_NO_ASSOCIATION_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        planned_family=len(t),status_counts=t.status.value_counts().to_dict(),scope_sha256=sha(O/'private_all572_planned_pair_designs.tsv'),
        mask_sha256=sha(O/'private_all16_parent_shared_gene_masks.tsv'),plan_sha256=sha(O/'plan_before_any_two_end_effect.json'),new_effect_tests=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
