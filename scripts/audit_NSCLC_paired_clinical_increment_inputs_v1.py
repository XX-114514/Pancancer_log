"""Independent merge-based all-scope input audit, without clinical effects."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.linalg import qr

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);o=r/'NSCLC_paired_programme_clinical_increment_preflight_20261016_v1'
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_NSCLC_paired_clinical_increment_v1/launch_receipt.json').read_text())
    deadline=time.monotonic()+3600
    while not (o/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();st=p.read_text().rsplit(')',1)[1].split()
        assert st[19]==str(receipt['proc_start_ticks']) and st[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((o/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(o/n)==h
    plan=json.loads((o/'plan_before_any_paired_clinical_increment_effect.json').read_text());assert plan['family']==572
    p=r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv';assert sha(p)==plan['original_scope_sha256']
    parent=pd.read_csv(p,sep='\t');scope=pd.read_csv(o/'private_all572_scope.tsv',sep='\t');assert len(scope)==len(parent)==572
    keys=['record','context','program','author_subtype_A','author_subtype_B','family_A','family_B']
    assert scope[keys].equals(parent[keys]) and not scope.duplicated(keys[1:5]).any()
    sp=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz'
    cp=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv'
    pp=r/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    for k,p in [('scores',sp),('case_eligibility',cp),('parameters',pp)]:assert sha(p)==plan['input_sha256'][k]
    scores=pd.read_csv(sp,sep='\t');c=pd.read_csv(cp,sep='\t');assert c.sampleID.is_unique
    rawclinical=pd.read_csv(r/'NSCLC_all_case_clinical_scope_20261016_v1/private_all243_author_case_fields.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('sampleID')
    expected=rawclinical.pathological_response.map({'pCR':1.,'MPR':1.,'non-MPR':0.}).reindex(c.sampleID)
    assert np.allclose(expected,c.response,equal_nan=True)
    params=pd.read_csv(pp,sep='\t');nr=r/'NSCLC_full_native_gene_units_20261016_v1'
    genes=pd.Index(pd.read_csv(nr/'private_native_gene_axis.tsv',sep='\t').gene);m=np.load(nr/'mean_logCP10K.npy',mmap_mode='r');masked={}
    for context,program in scope[['context','program']].drop_duplicates().itertuples(index=False):
        aa=params[params.context.eq(context)&params.program.eq(program)&params.type.eq('typeA')];bb=params[params.context.eq(context)&params.program.eq(program)&params.type.eq('typeB')]
        shared=set(aa.gene).intersection(bb.gene)
        for typ,p in [('typeA',aa),('typeB',bb)]:
            selected=~p.gene.isin(shared);present=genes.get_indexer(p.gene)>=0;w=np.abs(p.weight.to_numpy())
            ok=selected.any() and w[selected].sum()/w.sum()>=.8 and w[selected&present].sum()/w[selected].sum()>=.95 and (selected&present).sum()>=20
            if ok:
                p=p.loc[selected&present];f=scores[scores.context.eq(context)&scores.program.eq(program)&scores.type.eq(typ)]
                values=((np.asarray(m[np.ix_(f.unit_index.to_numpy(int),genes.get_indexer(p.gene))])-p.gene_mean.to_numpy())/p.gene_sd.to_numpy())@p.weight.to_numpy()
                masked[(context,program,typ)]=pd.Series(values,index=f.unit_index)
    hashes=pd.read_csv(o/'private_ready_design_SHA256.tsv',sep='\t').set_index('record').sha256;checked=0;maximum=0.
    for row in scope.itertuples():
        f=scores[scores.context.eq(row.context)&scores.program.eq(row.program)]
        a=f[f.type.eq('typeA')&f.sub_cell_type.eq(row.author_subtype_A)].drop(columns=['context','program','type','sub_cell_type'])
        b=f[f.type.eq('typeB')&f.sub_cell_type.eq(row.author_subtype_B)].drop(columns=['context','program','type','sub_cell_type'])
        merged=a.merge(b,on='sampleID',suffixes=('_A','_B'),validate='one_to_one').merge(c,on='sampleID',validate='one_to_one')
        t=merged[merged.n_cells_A.ge(20)&merged.n_cells_B.ge(20)&merged.clinical_eligible].sort_values('sampleID')
        n=len(t);assert n==row.eligible_cases and t.sampleID.is_unique
        y=t.response.to_numpy(float);assert int((y==0).sum())==row.response0 and int((y==1).sum())==row.response1
        cov=t[COMMON+[k+'_A' for k in UNIT]+[k+'_B' for k in UNIT]].to_numpy(float)
        raw=np.column_stack([cov,t.raw_weighted_score_A.to_numpy(),t.raw_weighted_score_B.to_numpy()])
        sd=raw.std(0) if n else np.zeros(18);mean=raw.mean(0) if n else np.zeros(18)
        d=np.column_stack([np.ones(n),(raw-mean)/np.where(sd>1e-12,sd,1)])
        rank=int(np.linalg.matrix_rank(d)) if n else 0;condition=float(np.linalg.cond(d)) if n else np.inf
        groups=np.bincount(y.astype(int),minlength=2);status='READY_NO_CLINICAL_EFFECT_COMPUTED'
        if n<50 or groups.min()<10:status='INSUFFICIENT_MATCHED_CASE_OR_RESPONSE_GROUP'
        elif not np.isfinite(d).all() or (sd<=1e-12).any() or rank!=19:status='NONFINITE_CONSTANT_OR_NONFULL_DESIGN'
        elif np.sum(np.linalg.svd(d,compute_uv=False)>np.linalg.svd(d,compute_uv=False)[0]*1e-12)!=rank:status='NUMERICALLY_AMBIGUOUS_RANK'
        elif n-rank<10 or condition>=1e8:status='LOW_DF_OR_ILL_CONDITIONED'
        elif np.max(np.sum(d*np.linalg.pinv(d).T,axis=1))>=.999999:status='EXTREME_LEVERAGE'
        assert row.status==status and row.design_rank==rank
        if status!='READY_NO_CLINICAL_EFFECT_COMPUTED':assert row.record not in hashes;continue
        path=o/'private_designs'/(row.record+'.npz');assert sha(path)==hashes.loc[row.record];z=np.load(path,allow_pickle=False)
        assert np.array_equal(z['cases'],t.sampleID.to_numpy(str)) and np.array_equal(z['Y'],y)
        assert np.allclose(z['D'],d,atol=1e-12) and np.allclose(z['C_raw'],cov,atol=1e-12)
        assert np.allclose(z['raw_sd'],sd) and np.allclose(z['raw_means'],mean)
        # QR directly verifies complete 19-column geometry without fitting the clinical outcome.
        q,rr,pivot=qr(d,mode='economic',pivoting=True);assert np.min(np.abs(np.diag(rr)))>1e-10
        assert np.allclose(np.sum(q*q,axis=1),np.sum(d*np.linalg.pinv(d).T,axis=1),atol=1e-9)
        ma=masked.get((row.context,row.program,'typeA'));mb=masked.get((row.context,row.program,'typeB'))
        for key,series,unit in [('masked_A',ma,t.unit_index_A),('masked_B',mb,t.unit_index_B)]:
            v=series.reindex(unit).to_numpy() if ma is not None and mb is not None else np.full(n,np.nan)
            assert np.allclose(z[key],v,atol=1e-9,equal_nan=True)
            if np.isfinite(v).any():maximum=max(maximum,float(np.nanmax(np.abs(z[key]-v))))
        checked+=1
    assert checked==len(hashes)
    proof=dict(status='PASS_ALL572_MATCHED_CASE_OUTCOME_FULL_CONTROLS_ORIGINAL_PROGRAMMES_AND_NONSHARED_NATIVE_READOUTS_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
               audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(o/'private_terminal.json'),ready_designs_checked=checked,maximum_nonshared_score_difference=maximum,effects_computed=False)
    path=o/'private_independent_input_audit.json';assert not path.exists();path.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
