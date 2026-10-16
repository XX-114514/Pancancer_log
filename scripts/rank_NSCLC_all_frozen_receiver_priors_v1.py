"""Rank fixed priors for all measured receiver programmes, without LR effects.

All expression/weight/ranking tables remain in the private server run root.
Custom AP is not official NicheNet activity, a P value or a causal LR test.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,pyreadr
from sklearn.metrics import average_precision_score,roc_auc_score

PRIOR_SHA={
    'lr_network_human_21122021.rds':'47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094',
    'ligand_target_matrix_nsga2r_final.rds':'699fce17ff65e2511277359696306bb130aa833bd7f30dd9290ad9efd9dc9c5d',
}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def rds(p):
    x=pyreadr.read_r(str(p));assert len(x)==1
    return next(iter(x.values()))

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    N=R/'NSCLC_full_native_gene_units_20261016_v1';P=R/'NSCLC_context_native_programme_projection_20261016_v1'
    assert json.loads((P/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL56_SCOPE_')
    priors=R/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for f,h in PRIOR_SHA.items():assert sha(priors/f)==h
    O=R/'NSCLC_all_frozen_receiver_prior_ranking_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),prior_sha256=PRIOR_SHA,
        scope='All24 measurable programme slots, all1226 prior ligands per signed gene set. No clinical, coordination or LR q-based programme filter.',
        expression_background='All243 native cases and all original author subtypes in each compatible lineage; no clinical outcome read. Native cells pooled only for background/receptor detection, never independent inference units.',
        family_definitions={'Macrophage':'Mphi author prefix','CD8':'CD8T_ prefix','Treg':'CD4T_Treg_ prefix'},
        background='>=10% pooled native cells detect gene, intersect original per-slot discovery gene universe, unique native symbols and fixed prior target genes.',
        foreground='Original frozen nonzero-weight genes split by original sign. At least10 expressed prior genes, no DEG or outcome-based gene reselection.',
        ranking='sklearn discrete average precision minus foreground prevalence, rawAP and AUROC descriptive. Not official NicheNet activity or P/FDR.',
        candidates='Before sender-response effects, freeze top20 positive AP-excess ligands per sign among ligands having at least one literal flat-prior receptor detected in>=10% pooled receiver cells; lexical ties. Retain20 slots even if unavailable. No rank21 replacement after qualification.',
        reference='All literal shortlist ligand edges in existing4624 curated table; all receptor components native and each detected>=10% pooled lineage cells reported. Marginal component detection does NOT establish same-cell complex.',
        followup='This is a prior/input inventory, not a new LR inferential family. Native per-case coverage, joint receptor detection, controls and full effect family must be frozen before sender-response effects.',
        limitations='All source cases previously explored; new ranking is exploratory, ligand-target prior unsigned and not receptor-specific. No directed, clinical-predictive or mechanistic validation.')
    (O/'plan_before_any_ranking_or_sender_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    units=pd.read_csv(N/'private_all_native_units_and_quality.tsv',sep='\t')
    assert units.unit_index.tolist()==list(range(len(units)))
    genes=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene);assert genes.is_unique
    hits=np.load(N/'positive_cell_counts.npy',mmap_mode='r');assert hits.shape==(len(units),len(genes))
    pools={};poolrows=[];totalrows={}
    prefixes={'Macrophage':'M\u03c6_','CD8':'CD8T_','Treg':'CD4T_Treg_'}
    for family,prefix in prefixes.items():
        selected=units[units.sub_cell_type.str.startswith(prefix)];ix=selected.unit_index.to_numpy(int)
        total=int(selected.n_cells.sum());assert total>0
        count=np.zeros(len(genes),dtype=np.int64)
        for start in range(0,len(genes),512):
            end=min(start+512,len(genes));count[start:end]=np.asarray(hits[np.ix_(ix,np.arange(start,end))]).sum(axis=0,dtype=np.int64)
        assert (count<=total).all();pools[family]=count*10>=total;totalrows[family]=count
        poolrows.extend(dict(source_family=family,gene=g,positive_cells=int(n),native_cells=total,expressed10pct=bool(ok)) for g,n,ok in zip(genes,count,pools[family]))
    pd.DataFrame(poolrows).to_csv(O/'private_all_lineage_native_detection.tsv.gz',sep='\t',index=False,compression='gzip')
    targets=rds(priors/'ligand_target_matrix_nsga2r_final.rds');lr=rds(priors/'lr_network_human_21122021.rds')
    targets.index=targets.index.astype(str);targets.columns=targets.columns.astype(str)
    assert targets.index.is_unique and targets.columns.is_unique and len(targets.columns)==1226
    assert np.isfinite(targets.to_numpy()).all() and (targets.to_numpy()>=0).all()
    scope=pd.read_csv(P/'private_all56_projection_scope.tsv',sep='\t');scope=scope[scope.status.str.startswith('PARTIAL_FROZEN')].copy();assert len(scope)==24
    params_path=R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    params=pd.read_csv(params_path,sep='\t');ranking=[];shortlist=[];backgrounds=[];coverage=[]
    receptor_by_family={}
    for family,detected in pools.items():
        expressed=set(genes[detected]);receptor_by_family[family]={}
        for ligand in targets.columns:
            receptors=sorted(set(lr.loc[lr['from'].astype(str).eq(ligand),'to'].astype(str)))
            receptor_by_family[family][ligand]=[r for r in receptors if r in expressed]
    for r in scope.sort_values(['context','type','program']).itertuples():
        signature=f'{r.context}_{r.type}_{r.program}'
        discovery=set((V/'state_after_DIALOGUE_inputs_20261016_v3'/f'context_{r.context:04d}'/(r.type+'_discovery_genes.tsv')).read_text().splitlines())
        expressed=set(genes[pools[r.source_family]]);bg=sorted(discovery&expressed&set(targets.index))
        p=params[params.context.eq(r.context)&params.type.eq(r.type)&params.program.eq(r.program)]
        assert len(p)>0 and p.gene.is_unique;weights=p.set_index('gene').weight
        backgrounds.extend(dict(signature=signature,source_family=r.source_family,gene=g,programme_sign='positive' if weights.get(g,0)>0 else 'negative' if weights.get(g,0)<0 else 'background') for g in bg)
        B=targets.loc[bg].to_numpy(float)
        for sign in ['positive','negative']:
            fg=set(p.loc[p.weight.gt(0) if sign=='positive' else p.weight.lt(0),'gene'])
            y=np.array([g in fg for g in bg],int);n=int(y.sum());ok=n>=10 and n<len(bg)
            coverage.append(dict(signature=signature,source_family=r.source_family,programme_sign=sign,original_signed_genes=len(fg),expressed_prior_genes=n,background_genes=len(bg),ranking_available=ok))
            rows=[]
            for j,ligand in enumerate(targets.columns):
                compatible=receptor_by_family[r.source_family][ligand]
                ap=float(average_precision_score(y,B[:,j])) if ok else np.nan
                rows.append(dict(signature=signature,context=r.context,type=r.type,program=r.program,source_family=r.source_family,programme_sign=sign,ligand=ligand,
                    foreground_genes=n,background_genes=len(bg),average_precision=ap,AP_excess_over_prevalence=ap-n/len(bg) if ok else np.nan,
                    AUROC=float(roc_auc_score(y,B[:,j])) if ok else np.nan,receiver_flat_receptors='|'.join(compatible),receiver_compatible=bool(compatible),ranking_available=ok))
            table=pd.DataFrame(rows);eligible=table[table.receiver_compatible & table.AP_excess_over_prevalence.gt(0)].sort_values(['AP_excess_over_prevalence','ligand'],ascending=[False,True]).head(20)
            selected=set(eligible.ligand);table['prospective_top20']=table.ligand.isin(selected);ranking.extend(table.to_dict('records'))
            chosen=eligible.to_dict('records')
            for position in range(20):
                slot=dict(signature=signature,context=r.context,type=r.type,program=r.program,source_family=r.source_family,programme_sign=sign,prospective_position=position+1)
                slot.update(chosen[position] if position<len(chosen) else dict(ligand='',status='NO_PRIOR_CANDIDATE_UNDER_FROZEN_RULE'))
                if position<len(chosen):slot['status']='FROZEN_PRIOR_SLOT_REQUIRES_NATIVE_INTERFACE_QUALIFICATION'
                shortlist.append(slot)
        (O/'private_progress.json').write_text(json.dumps(dict(last_signature=signature,finished_signatures=len(coverage)//2,total_signatures=24))+'\n')
    pd.DataFrame(backgrounds).to_csv(O/'private_all24_native_prior_backgrounds.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(coverage).to_csv(O/'private_all48_signed_programme_coverage.tsv',sep='\t',index=False)
    ranks=pd.DataFrame(ranking);slots=pd.DataFrame(shortlist);assert len(ranks)==24*2*1226 and len(slots)==960
    ranks.to_csv(O/'private_all_ligand_rankings.tsv.gz',sep='\t',index=False,na_rep='NA',compression='gzip')
    slots.to_csv(O/'private_all960_prospective_prior_slots.tsv',sep='\t',index=False,na_rep='NA')
    refpath=R/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv';ref=pd.read_csv(refpath,sep='\t',keep_default_na=False);assert len(ref)==4624
    accepted={'simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'};edges=[]
    selected=slots[slots.ligand.ne('')][['signature','source_family','ligand']].drop_duplicates()
    for row in selected.itertuples():
        family=row.source_family;expressed=set(genes[pools[family]])
        for edge in ref[ref.ligand.eq(row.ligand)].itertuples():
            parts=str(edge.receptor).split('_');measured=all(g in genes for g in parts);marginal=measured and all(g in expressed for g in parts)
            edges.append(dict(signature=row.signature,source_family=family,ligand=row.ligand,receptor=edge.receptor,reference_status=edge.reference_status,
                receptor_components_measured=measured,all_component_pooled_marginals10pct=marginal,reference_and_marginal_compatible=edge.reference_status in accepted and marginal,
                same_cell_RNA_verified=False,complex_assembly_verified=False))
    pd.DataFrame(edges).to_csv(O/'private_all_shortlist_reference_interfaces.tsv',sep='\t',index=False)
    files=['private_all_lineage_native_detection.tsv.gz','private_all24_native_prior_backgrounds.tsv.gz','private_all48_signed_programme_coverage.tsv','private_all_ligand_rankings.tsv.gz','private_all960_prospective_prior_slots.tsv','private_all_shortlist_reference_interfaces.tsv']
    report=dict(status='ALL24_RECEIVER_PRIOR_RANKINGS_AND960_SLOTS_FROZEN_NO_SENDER_RESPONSE_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(__file__),parameters_sha256=sha(params_path),curated_reference_sha256=sha(refpath),ranking_records=len(ranks),fixed_slots=len(slots),
        slots_with_candidate=int(slots.ligand.ne('').sum()),unique_ligands=int(slots.loc[slots.ligand.ne(''),'ligand'].nunique()),reference_interfaces=len(edges),
        files={f:sha(O/f) for f in files},new_P_or_BH=False,directed_or_clinical_interaction_validated=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
