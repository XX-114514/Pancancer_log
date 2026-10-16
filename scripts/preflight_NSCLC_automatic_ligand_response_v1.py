"""Freeze all original pair orientations x prior slots before any LR effects.

No filtering by clinical/coordination significance and no effect calculation.
Inputs remain private. Shared sign slots map to one unique ligand hypothesis.
"""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_prior_joint_receptor_RNA_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_prior_joint_receptor_RNA_v1/launch_receipt.json').read_text());deadline=time.monotonic()+18000
    while not (I/'private_independent_audit.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(15)
    assert json.loads((I/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL_UNIT_GENE_MARGINALS_')
    O=R/'NSCLC_automatic_ligand_response_preflight_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    P=R/'NSCLC_all_frozen_receiver_prior_ranking_20261016_v1';prior_path=P/'private_all960_prospective_prior_slots.tsv'
    scopes=R/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv'
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        input_sha256={'prior_slots':sha(prior_path),'original_pair_scope':sha(scopes),'joint_RNA_audit':sha(I/'private_independent_audit.json')},
        scope='Every original572 author-subtype pair, BOTH receiver orientations, original40 signed prior slots per receiver programme. No filter by any prior clinical or coordination effect/P/q.',
        no_manual_LR_selection=True,no_lower_rank_replacement=True,
        deduplication='Same pair/orientation/ligand appears in both prior signs only once as an inferential hypothesis. Original logical slots and sign provenance retained, empty candidate slots preserved as no hypothesis.',
        cases='Existing223 complete recorded-PD1 case eligibility, matched sender/receiver author subtype each>=20 cells. Cases not all formally reviewed independent patients.',
        native_ligand='Single literal native ligand gene only; original mean of per-cell log1pCP10K; >=10 matched cases have any sender ligand RNA. True measured zeros retained, missing genes not imputed.',
        receptors='All original candidate literal reference/marginal-qualified interfaces with same-cell pooled1% RNA. Availability covariate is maximum native per-case same-cell positive fraction across these fixed receptors, never effect-selected.',
        target='Raw original frozen receiver gene score after removing ALL original sender programme genes plus current ligand and ALL qualified receptor components. Keep frozen weights/scales, no relearning. Require>=80% original L1 remains, >=95% remaining weight native, >=20 native genes.',
        nuisance=COMMON+['sender_'+s for s in UNIT]+['receiver_'+s for s in UNIT]+['receiver_own_ligand','max_qualified_receptor_joint_fraction'],
        main_geometry='Ranks of ligand/receiver scores; all18 raw standardized nuisance variables+intercept, full20-column design. No pathology outcome conditioning. >=50 cases, df>=10, condition<1e8, leverage<.999999; no automatic covariate removal.',
        sensitivity='Preserve sender raw frozen programme as an optional additional nuisance input for separately predeclared descriptive sensitivity; not a second selected main model.',
        effect_gate='NO correlations, coefficients, P/BH or randomization here. Independent full input/design audit plus a separate precision/complete-family inference contract required before effects.',
        boundary='Exploratory ligand-to-receiver-programme candidates. A shared ligand-target prior cannot establish a particular receptor, signaling direction, protein complex, clinical prediction or causal mechanism.')
    (O/'plan_before_any_ligand_response_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    parents=pd.read_csv(scopes,sep='\t');assert len(parents)==572
    prior=pd.read_csv(prior_path,sep='\t',keep_default_na=False)
    interfaces=pd.read_csv(I/'private_all_frozen_interfaces_joint_RNA_qualification.tsv',sep='\t')
    interfaces=interfaces[interfaces.native_joint1pct_qualified]
    joint=pd.read_csv(I/'private_all_native_unit_receptor_joint_counts.tsv.gz',sep='\t')
    joint_lookup={r:f.set_index('unit_index').joint_positive_cells for r,f in joint.groupby('receptor')}
    N=R/'NSCLC_full_native_gene_units_20261016_v1';native=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene)
    M=np.load(N/'mean_logCP10K.npy',mmap_mode='r')
    params=pd.read_csv(R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv',sep='\t')
    scores=pd.read_csv(R/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',sep='\t')
    c=pd.read_csv(R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',sep='\t').set_index('sampleID')
    source={(int(ctx),typ,program):f for (ctx,typ,program),f in scores.groupby(['context','type','program'])}
    params_lookup={(int(ctx),typ,program):f for (ctx,typ,program),f in params.groupby(['context','type','program'])}
    mask_cache={};logical=[];rows=[];design=O/'private_frozen_designs';design.mkdir(mode=0o700)
    for parent in parents.itertuples():
        for receive,send in [('typeA','typeB'),('typeB','typeA')]:
            signature=f'{parent.context}_{receive}_{parent.program}'
            ps=prior[prior.signature.eq(signature)].sort_values(['programme_sign','prospective_position']);assert len(ps)==40
            rf=source[(parent.context,receive,parent.program)];sf=source[(parent.context,send,parent.program)]
            rn=parent.author_subtype_A if receive=='typeA' else parent.author_subtype_B
            sn=parent.author_subtype_B if send=='typeB' else parent.author_subtype_A
            rf=rf[rf.sub_cell_type.eq(rn)].set_index('sampleID');sf=sf[sf.sub_cell_type.eq(sn)].set_index('sampleID')
            cases=sorted(set(rf.index)&set(sf.index));rr=rf.loc[cases];ss=sf.loc[cases];cc=c.loc[cases]
            use=rr.n_cells.ge(20).to_numpy()&ss.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool)
            rr=rr.loc[use];ss=ss.loc[use];cc=cc.loc[use];n=len(rr)
            mapping={}
            for slot in ps.itertuples():
                ligand=str(slot.ligand);hypothesis=''
                if ligand:
                    if ligand not in mapping:
                        hypothesis='LRD'+str(len(rows)+1).zfill(5);mapping[ligand]=hypothesis
                        row=dict(record=hypothesis,parent_pair=parent.record,context=parent.context,program=parent.program,receiver_type=receive,receiver_signature=signature,
                            sender_author_subtype=sn,receiver_author_subtype=rn,ligand=ligand,prior_signs='|'.join(sorted(set(ps.loc[ps.ligand.eq(ligand),'programme_sign']))),
                            matched_native_cases=n,status='READY_FOR_INDEPENDENT_DESIGN_AUDIT',target_retained_original_L1=np.nan,target_native_remaining_L1=np.nan,qualified_receptors='',design_rank=0,design_columns=20)
                        edges=interfaces[interfaces.signature.eq(signature)&interfaces.ligand.eq(ligand)]
                        receptors=sorted(set(edges.receptor));row['qualified_receptors']='|'.join(receptors)
                        x=None;y=None;Z=None
                        if ligand not in native:row['status']='LITERAL_LIGAND_NATIVE_GENE_NOT_MEASURED'
                        elif not receptors:row['status']='NO_FROZEN_REFERENCE_AND_JOINT_RNA_QUALIFIED_INTERFACE'
                        else:
                            cachekey=(parent.context,parent.program,receive,ligand)
                            if cachekey not in mask_cache:
                                pp=params_lookup[(parent.context,receive,parent.program)];sp=params_lookup[(parent.context,send,parent.program)]
                                drop=set(sp.gene)|{ligand}|{g for r in receptors for g in r.split('_')}
                                keep=~pp.gene.isin(drop);measured=native.get_indexer(pp.gene)>=0;w=np.abs(pp.weight.to_numpy())
                                retained=float(w[keep].sum()/w.sum());native_frac=float(w[keep&measured].sum()/w[keep].sum()) if keep.any() else 0.
                                valid=retained>=.8 and native_frac>=.95 and int((keep&measured).sum())>=20
                                pp=pp.loc[keep&measured];mask_cache[cachekey]=(valid,retained,native_frac,pp)
                            valid,retained,measured,pp=mask_cache[cachekey];row.update(target_retained_original_L1=retained,target_native_remaining_L1=measured)
                            if not valid:row['status']='RECEIVER_MASK_CHANGES_PROGRAMME_OR_NATIVE_COVERAGE_TOO_MUCH'
                            elif n<50:row['status']='INSUFFICIENT_MATCHED_CASES'
                            else:
                                ui=rr.unit_index.to_numpy(int);si=ss.unit_index.to_numpy(int);column=native.get_loc(ligand)
                                x=np.asarray(M[si,column]);own=np.asarray(M[ui,column])
                                y=np.asarray(M[np.ix_(ui,native.get_indexer(pp.gene))])@(pp.weight/pp.gene_sd).to_numpy()-float((pp.weight*pp.gene_mean/pp.gene_sd).sum())
                                rec=np.max(np.column_stack([joint_lookup[r].reindex(ui).to_numpy(float)/rr.n_cells.to_numpy(float) for r in receptors]),axis=1)
                                cov=np.column_stack([cc[COMMON].to_numpy(float),ss[UNIT].to_numpy(float),rr[UNIT].to_numpy(float),own,rec])
                                assert cov.shape==(n,18) and np.isfinite(cov).all()
                                sd=cov.std(0);Z=np.column_stack([np.ones(n),(cov-cov.mean(0))/np.where(sd>1e-12,sd,1)])
                                xx=rankdata(x);xx=(xx-xx.mean())/xx.std() if xx.std()>0 else np.zeros(n);D=np.column_stack([Z,xx]);rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D))
                                row.update(design_rank=rank,condition=condition,residual_df=n-rank,sender_RNA_positive_cases=int(np.sum(x>0)))
                                if np.sum(x>0)<10:row['status']='INSUFFICIENT_CASES_WITH_NATIVE_SENDER_LIGAND_RNA'
                                elif (sd<=1e-12).any() or np.linalg.matrix_rank(Z)!=19 or rank!=20:row['status']='CONSTANT_OR_COLLINEAR_REQUIRED_DESIGN'
                                elif n-rank<10 or condition>=1e8:row['status']='LOW_DF_OR_ILL_CONDITIONED'
                                elif np.ptp(y)<=1e-12:row['status']='NO_MASKED_RECEIVER_SIGNATURE_VARIATION'
                                elif np.max(np.sum(D*np.linalg.pinv(D).T,axis=1))>=.999999:row['status']='EXTREME_LEVERAGE'
                                if row['status']=='READY_FOR_INDEPENDENT_DESIGN_AUDIT':
                                    np.savez_compressed(design/(hypothesis+'.npz'),X=x,Y=y,Z=Z,covariates=cov,cases=rr.index.to_numpy(str),sender_units=si,receiver_units=ui,
                                        sender_programme=ss.raw_weighted_score.to_numpy(float),receiver_full_programme=rr.raw_weighted_score.to_numpy(float),response=cc.response.to_numpy(int))
                        rows.append(row)
                    else:hypothesis=mapping[ligand]
                logical.append(dict(parent_pair=parent.record,receiver_type=receive,receiver_signature=signature,programme_sign=slot.programme_sign,prior_position=int(slot.prospective_position),ligand=ligand,hypothesis=hypothesis))
        if len(rows)%500<80:(O/'private_progress.json').write_text(json.dumps(dict(parent_pair=parent.record,unique_hypotheses=len(rows),logical_slots=len(logical))))
    logical=pd.DataFrame(logical);table=pd.DataFrame(rows);assert len(logical)==572*2*40 and table.record.is_unique
    assert set(logical.loc[logical.hypothesis.ne(''),'hypothesis'])==set(table.record)
    logical.to_csv(O/'private_all45760_logical_slots.tsv.gz',sep='\t',index=False,compression='gzip')
    table.to_csv(O/'private_all_unique_ligand_response_designs.tsv',sep='\t',index=False,na_rep='NA')
    result=dict(status='ALL_ORIGINAL_PAIR_ORIENTATIONS_AUTOMATIC_PRIOR_INPUTS_FROZEN_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        logical_slots=len(logical),unique_hypotheses=len(table),status_counts=table.status.value_counts().to_dict(),
        logical_scope_sha256=sha(O/'private_all45760_logical_slots.tsv.gz'),unique_scope_sha256=sha(O/'private_all_unique_ligand_response_designs.tsv'),
        new_effects_or_P_or_BH=False,independent_design_audit_pending=True,precision_contract_before_effects_required=True)
    (O/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
