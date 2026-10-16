"""Sender-agnostic prior ranking for the fixed C1QC programme; no effect tests.

Custom AP ranks the frozen programme gene set against native expressed
background using the pinned NicheNet prior. It is not the official NicheNet
activity implementation and does not identify causal or receptor-specific axes.
"""
import hashlib,json,os
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import numpy as np,pandas as pd,pyreadr
from sklearn.metrics import average_precision_score,roc_auc_score
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs

PRIOR_SHA={
    'lr_network_human_21122021.rds':'47c971d2fbba4ecd0ba7485d1846a74054432a0d1a97ea3e6a79ae27d0da8094',
    'ligand_target_matrix_nsga2r_final.rds':'699fce17ff65e2511277359696306bb130aa833bd7f30dd9290ad9efd9dc9c5d',
}


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8388608),b''):h.update(block)
    return h.hexdigest()


def read_rds(path):
    objects=pyreadr.read_r(str(path));assert len(objects)==1
    matrix=next(iter(objects.values()));assert isinstance(matrix,pd.DataFrame)
    return matrix


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    priors=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for name,expected in PRIOR_SHA.items():assert sha(priors/name)==expected
    out=root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),prior_sha256=PRIOR_SHA,
        receiver='All original baseline Macro_C1QC multigene candidate patients with the frozen partial CRC macrophage programme.',
        receiver_selected_after_clinical_results=True,not_independent_clinical_validation=True,
        no_new_P_BH_or_expression_download=True,source='2026-10-14 immune atlas v1; original P3 v20261010 identities/9CM unchanged.',
        background='All native measured genes in >=10% of pooled eligible baseline C1QC cells, intersect unique physical prior-target columns.',
        gene_set='All nonzero original programme weights, positive/negative signs separately. No DEG, clinical-expression or effect-based reselection.',
        gene_aliases='Unique stored stable-ID aliases only; prior-target aliases sharing a physical column prefer original programme symbol, then canonical symbol, then lexical tie.',
        metric='Custom sklearn average precision minus positive-background prevalence; report rawAP, AUROC, Pearson and full ranking. Not official NicheNet activity.',
        receptor_eligibility='At least one literal flat NicheNet LR edge to a source-measured receptor detected in >=10% receiver cells. Flat prior does not prove assembled complex or direct binding.',
        prior_ranking='All prior ligands retained. Primary prospective follow-up top20 per gene-set sign among receptor-eligible ligands with positive AP excess; lexical ligand ties.',
        minimum_gene_set=10,top_ligands_per_sign=20,
        sender_scope='All observed baseline major lineages except unresolved/T_NK boundary. Remove C1QC receiver cells from Myeloid sender; >=20 cells in >=20 receiver-matched patients.',
        future_effects='Freeze candidate ligand x every eligible sender group before effects; measured zeros retained, unavailable measurement missing. Main required controls/identifiability predeclared separately before any correlations.',
        limitations='Unsigned target prior does not predict activation sign, distinguish alternative receptors, validate clinical value, or prove a sender exists.')
    (out/'plan_before_native_readout_or_ranking.json').write_text(json.dumps(plan,indent=2)+'\n')
    g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
    receiver_table=pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv',sep='\t')
    receiver=receiver_table[receiver_table.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id').sort_index()
    assert receiver.index.is_unique
    pre=cells.timepoint.eq('pre')&cells.patient_id.isin(receiver.index)&cells.library_id.isin(eligible.library_id)
    receiver_mask=pre&cells.major.eq('Myeloid')&cells.subtype.eq('Macro_C1QC')&cells.evidence.eq('multigene_subtype_candidate')
    selected=cells[receiver_mask];count=selected.groupby('patient_id').size()
    assert set(count.index)==set(receiver.index) and count.ge(20).all()
    assert np.array_equal(count.loc[receiver.index],receiver.native_subtype_cells)
    positions=np.flatnonzero(receiver_mask.to_numpy());matrix=ad.io.sparse_dataset(g['layers/counts'])
    hits=np.zeros(len(var),np.int64)
    for start in range(0,len(positions),2000):
        block=matrix[source_rows[positions[start:start+2000]],:].tocsr();block.eliminate_zeros()
        assert np.isfinite(block.data).all() and (block.data>0).all()
        hits+=np.asarray((block>0).sum(axis=0)).ravel()
    expressed=measured&(hits*10>=len(positions))
    lookup=defaultdict(set)
    for column,(name,row) in enumerate(var.iterrows()):
        if not measured[column] or str(name).startswith('AMBIGUOUS::'):continue
        aliases={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':aliases.update(str(row.get('source_gene_aliases','')).split('|'))
        for alias in aliases:
            if alias not in ['','nan','None']:lookup[alias].add(column)
    weight_path=root/'HNSCC_frozen_MCP_gene_axis_20261014_v2/Macnode18_private_frozen_gene_columns.tsv'
    weights=pd.read_csv(weight_path,sep='\t')
    assert weights.gene.is_unique
    native_programme={str(r.gene):int(r.native_column) for r in weights.itertuples() if r.native_column>=0}
    for gene,column in native_programme.items():assert lookup[gene]=={column}
    lr=read_rds(priors/'lr_network_human_21122021.rds')
    targets=read_rds(priors/'ligand_target_matrix_nsga2r_final.rds')
    targets.index=targets.index.astype(str);targets.columns=targets.columns.astype(str)
    assert targets.index.is_unique and targets.columns.is_unique
    assert np.isfinite(targets.to_numpy()).all() and (targets.to_numpy()>=0).all()
    choices=defaultdict(list)
    for gene in targets.index:
        candidates=lookup.get(gene,set())
        if len(candidates)==1:
            column=next(iter(candidates))
            if expressed[column]:choices[column].append(gene)
    rows=[]
    for column,aliases in choices.items():
        canonical=str(var.iloc[column].get('gene_symbol',''))
        gene=sorted(aliases,key=lambda symbol:(symbol not in native_programme,symbol!=canonical,symbol))[0]
        matching=weights[weights.native_column.eq(column)]
        assert len(matching)<=1
        sign='positive' if len(matching) and matching.weight.iloc[0]>0 else 'negative' if len(matching) and matching.weight.iloc[0]<0 else 'background'
        rows.append(dict(gene=gene,native_column=column,receiver_detected_cells=int(hits[column]),programme_sign=sign))
    background=pd.DataFrame(rows).sort_values('gene').reset_index(drop=True)
    background.to_csv(out/'private_native_prior_background_and_programme.tsv',sep='\t',index=False)
    B=targets.loc[background.gene].to_numpy(float)
    coverage=[];rankings=[];receptor_rows=[]
    eligible_ligands={}
    for ligand in targets.columns:
        receptors=sorted(set(lr.loc[lr['from'].astype(str).eq(ligand),'to'].astype(str)))
        compatible=[]
        for receptor in receptors:
            candidates=lookup.get(receptor,set())
            if len(candidates)==1:
                column=next(iter(candidates))
                if expressed[column]:compatible.append(receptor)
        eligible_ligands[ligand]=compatible
        receptor_rows.append(dict(ligand=ligand,flat_prior_receptors='|'.join(receptors),receiver_expressed_receptors='|'.join(compatible),
            receiver_compatible=bool(compatible),assembled_complex_verified=False))
    for sign in ['positive','negative']:
        gene_set=background.programme_sign.eq(sign).to_numpy(int);n=int(gene_set.sum())
        total=int((weights.weight.gt(0) if sign=='positive' else weights.weight.lt(0)).sum())
        coverage.append(dict(programme_sign=sign,original_programme_genes=total,expressed_prior_background_genes=n,background_genes=len(background),
            gene_fraction_in_prior_expressed_background=n/total if total else None,ranking_available=n>=10 and n<len(background)))
        if n<10 or n==len(background):continue
        for j,ligand in enumerate(targets.columns):
            vector=B[:,j];ap=float(average_precision_score(gene_set,vector));auc=float(roc_auc_score(gene_set,vector))
            pearson=float(np.corrcoef(gene_set,vector)[0,1]) if np.std(vector)>0 else np.nan
            rankings.append(dict(programme_sign=sign,ligand=ligand,programme_genes_in_background=n,background_genes=len(background),
                average_precision=ap,AP_excess_over_prevalence=ap-n/len(background),AUROC=auc,Pearson=pearson,
                receiver_prior_compatible=bool(eligible_ligands[ligand]),receiver_expressed_receptors='|'.join(eligible_ligands[ligand]),
                full_complex_binding_verified=False))
    ranking=pd.DataFrame(rankings)
    if len(ranking):
        ranking=ranking.sort_values(['programme_sign','AP_excess_over_prevalence','ligand'],ascending=[True,False,True])
        ranking['all_prior_rank']=ranking.groupby('programme_sign').cumcount()+1
        ranking['prospective_top20']=False
        for sign,group in ranking.groupby('programme_sign'):
            valid=group[group.receiver_prior_compatible&group.AP_excess_over_prevalence.gt(0)]
            ranking.loc[valid.head(20).index,'prospective_top20']=True
    ranking.to_csv(out/'private_all_ligand_prior_rankings.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(receptor_rows).to_csv(out/'private_all_ligand_receptor_compatibility.tsv',sep='\t',index=False)
    sender_counts=[];sender_summary=[]
    donor_cells=cells[pre];sender_cells=donor_cells[~donor_cells.major.isin(['Unresolved','T_NK_boundary'])&~receiver_mask.loc[donor_cells.index]]
    for major,group in sender_cells.groupby('major',observed=True):
        counts=group.groupby('patient_id').size().reindex(receiver.index,fill_value=0)
        name=str(major)+'_excluding_C1QC' if str(major)=='Myeloid' else str(major)
        for patient,value in counts.items():sender_counts.append(dict(sender_group=name,original_major=str(major),patient_id=patient,n_cells=int(value)))
        sender_summary.append(dict(sender_group=name,patients_with_any_cells=int(counts.gt(0).sum()),patients_with_atleast20_cells=int(counts.ge(20).sum()),
            prospective_sender_eligible=bool(counts.ge(20).sum()>=20),receiver_cells_removed=major=='Myeloid'))
    pd.DataFrame(sender_counts).to_csv(out/'private_all_sender_patient_capture.tsv',sep='\t',index=False)
    pd.DataFrame(sender_summary).to_csv(out/'private_sender_capture_summary.tsv',sep='\t',index=False)
    report=dict(status='HNSCC_FIXED_RECEIVER_ALL_PRIOR_LIGANDS_RANKED_NO_ASSOCIATION_TESTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),native_receiver_cells=len(positions),receiver_patients=len(receiver),
        prior_target_ligands=len(targets.columns),receiver_prior_background_genes=len(background),programme_coverage=coverage,
        prospective_top20_slots=int(ranking.prospective_top20.sum()) if len(ranking) else 0,sender_capture=sender_summary,
        prior_ranking_sha256=sha(out/'private_all_ligand_prior_rankings.tsv'),programme_weights_sha256=sha(weight_path),
        metric_is_custom_AP_not_official_NicheNet_activity=True,not_receptor_specific_or_causal=True,
        clinical_tests_added=False,sender_receiver_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    if len(ranking):print(ranking[ranking.prospective_top20].head(8)[['programme_sign','ligand','average_precision','AP_excess_over_prevalence','receiver_expressed_receptors']].to_json(orient='records'))


if __name__=='__main__':main()
