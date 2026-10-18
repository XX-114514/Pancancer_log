"""All-prior native RNA qualification and fixed automatic, exploratory ranking."""
import gzip
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import PRIORS, read_rds


def pieces(name):
    # A complex is measured through every constituent, not a fictitious
    # feature bearing the full concatenated complex name.
    return set(str(name).split('_'))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,64*1024**2)
    acceptance=root/'immune_RNA_MCP_original_node_decomposition_audit_20261018_v1/private_terminal.json'
    assert sha(acceptance)=='f109cecf33be6ad8a937bce8c50d567aa1426f1a1e1ea05a0c20335ce4e123d7'
    prior=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for name,digest in PRIORS.items():assert sha(prior/name)==digest
    refpath=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(refpath)=='df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    out=root/'immune_RNA_MCP_all_ligand_qualification_20261018_v2'
    assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        accepted_composition_followup_sha256=sha(acceptance),prior_sha256=PRIORS,reference_sha256=sha(refpath),
        scope='All1226 prior ligands x2 programmes x3 receiver views x2 signed sets=14712 ranking positions. Full old frozen19018 selected cells, native ligand/receptor RNA, no clinical or new effect family.',
        native='Exact source-measured unique symbol/stable-ID aliases. All complex components. Native fullRNA denominator checked; measured zero retained, missing gene unknown. Fit-excluded34 markers may be assayed but are not added to programme training.',
        qualification='At least one literal accepted reviewed reference edge; all ligand/receptor components measurable in all3 sources. Discovery pooled selected-cell receptor and at least one other sender view detection>=1pct,10pct sensitivity.',
        background='Original5000 per-view discovery candidates intersect unsigned prior targets. Remove ligand and ALL flat/reference receptor components and physical aliases before signed gene sets and continuous vectors.',
        ranking='All signed nonzero genes, minimum10 signed genes and200 background. Custom AP minus prevalence; AUROC and coefficient Pearson/Spearman sensitivity. Not official NicheNet activity, no P/FDR or inferred inhibition.',
        shortlist='Top10 AP-excess positive, RNA/reference-eligible ligands per programme/view/sign, lexical ties. All accepted receptors/senders retained; no rank11 replacement.',
        limits='Candidate generation after aggregate source results. No sender->receiver causality, receptor assembly, clinical or independent-cohort validation. All learned candidates and weights remain server-private.')
    (out/'plan_before_native_assay_or_automatic_ranking.json').write_text(json.dumps(plan,indent=2)+'\n')
    target=read_rds(prior/'ligand_target_matrix_nsga2r_final.rds')
    target.index=target.index.astype(str);target.columns=target.columns.astype(str)
    assert target.index.is_unique and target.columns.is_unique and len(target.columns)==1226
    assert np.isfinite(target.to_numpy()).all() and (target.to_numpy()>=0).all()
    flat=read_rds(prior/'lr_network_human_21122021.rds');flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    ref=pd.read_csv(refpath,sep='\t',keep_default_na=False);assert len(ref)==4624
    ligands=sorted(target.columns.tolist());exclusions={};requested=set()
    for ligand in ligands:
        receptors=set(flat.loc[flat['from'].eq(ligand),'to'])|set(ref.loc[ref.ligand.eq(ligand),'receptor'])
        symbols=pieces(ligand)
        for receptor in receptors:symbols.update(pieces(receptor))
        exclusions[ligand]=symbols;requested.update(symbols)
    genes=sorted(requested);ng=len(genes);gi={g:i for i,g in enumerate(genes)}
    src=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1'
    assert sha(src/'private_terminal.json')=='2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656'
    md=pd.read_csv(src/'private_selected_native_cell_depth.tsv',sep='\t')
    assert len(md)==19018 and md.global_row.is_unique and md.global_row.is_monotonic_increasing
    units=md[['fit_sample','dataset','partition','view']].drop_duplicates().sort_values(['dataset','fit_sample','view']).reset_index(drop=True)
    assert len(units)==198
    key=pd.MultiIndex.from_frame(units[['fit_sample','view']]);u=key.get_indexer(pd.MultiIndex.from_frame(md[['fit_sample','view']]))
    assert (u>=0).all()
    sources=sorted(md.dataset.unique());views=sorted(md.view.unique())
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    order=list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)
    symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    columns={};maps={}
    for source in sources:
        names=defaultdict(set)
        for col in np.flatnonzero(measured[:,order.index(source)]):
            if canonical[col].startswith('AMBIGUOUS::') or identity[col]=='ambiguous_source_feature':continue
            labels={symbol[col],canonical[col]}
            if identity[col]=='stable_Ensembl_ID':labels.update(aliases[col].split('|'))
            for name in labels:
                if name in requested:names[name].add(int(col))
        maps[source]={name:next(iter(pos)) for name,pos in names.items() if len(pos)==1}
        columns[source]=np.array([maps[source].get(g,-1) for g in genes],int)
    gene_columns=pd.DataFrame(dict(gene=genes,**columns));gene_columns.to_csv(out/'private_all_prior_native_gene_columns.tsv',sep='\t',index=False)
    hits=np.zeros((198,ng),np.int64);sums=np.zeros((198,ng));cells=np.bincount(u,minlength=198)
    counts=z['layers/counts'];nc,nf=map(int,counts.attrs['shape']);rows=md.global_row.to_numpy(int);read=0
    for start in range(0,nc,8192):
        stop=min(nc,start+8192);a,b=np.searchsorted(rows,[start,stop])
        if b>a:
            ptr=np.asarray(counts['indptr'][start:stop+1],np.int64)
            x=csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,nf))[rows[a:b]-start].copy()
            x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert np.isfinite(x.data).all() and (x.data>0).all() and (x.data==np.floor(x.data)).all()
            total=np.asarray(x.sum(axis=1)).ravel();assert np.array_equal(total,md.fullRNA.iloc[a:b].to_numpy())
            assert np.array_equal(np.diff(x.indptr),md.full_detected_genes.iloc[a:b].to_numpy())
            for source in md.dataset.iloc[a:b].unique():
                take=np.flatnonzero(md.dataset.iloc[a:b].to_numpy()==source);col=columns[source];valid=np.flatnonzero(col>=0)
                v=x[take][:,col[valid]].astype(float);unique,inverse=np.unique(u[a:b][take],return_inverse=True)
                agg=csr_matrix((np.ones(len(take)),(inverse,np.arange(len(take)))),shape=(len(unique),len(take)))
                positive=(v>0).astype(np.int64);hits[np.ix_(unique,valid)]+=np.asarray((agg@positive).toarray(),np.int64)
                v.data=np.log1p(v.data*np.repeat(10000/total[take],np.diff(v.indptr)))
                sums[np.ix_(unique,valid)]+=np.asarray((agg@v).toarray())
            read+=b-a
        if start//8192%32==0 or stop==nc:
            (out/'private_native_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),original_rows_scanned=int(stop),selected_cells_read=int(read)))+'\n')
    assert read==19018 and np.all(hits<=cells[:,None]) and cells.min()>=20
    means=sums/cells[:,None]
    for i,source in enumerate(units.dataset):means[i,columns[source]<0]=np.nan
    np.savez_compressed(out/'private_all198_native_ligand_receptor_assays.npz',genes=np.array(genes),positive_cells=hits,mean_logCP10K=means,cells=cells)
    units.assign(selected_cells=cells).to_csv(out/'private_all198_original_case_view_units.tsv',sep='\t',index=False)
    pooled={};totals={}
    for view in views:
        take=units.view.eq(view)&units.partition.eq('reviewed_case_discovery')
        pooled[view]=hits[take].sum(axis=0);totals[view]=int(cells[take].sum())
    def available(parts):
        return all(g in gi and all(columns[s][gi[g]]>=0 for s in sources) for g in parts)
    def detect(parts,view,denominator):
        return available(parts) and all(int(pooled[view][gi[g]])*denominator>=totals[view] for g in parts)
    accepted={'simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'}
    edges=[]
    for view in views:
        for r in ref[ref.ligand.isin(ligands)].itertuples(index=False):
            rp=pieces(r.receptor);lp=pieces(r.ligand)
            eligible= r.reference_status in accepted and detect(rp,view,100)
            senders=[s for s in views if s!=view and detect(lp,s,100)] if eligible else []
            strict=[s for s in senders if detect(rp,view,10) and detect(lp,s,10)]
            edges.append(dict(receiver_view=view,ligand=r.ligand,receptor=r.receptor,reference_status=r.reference_status,
                receptor_components_measured=available(rp),ligand_components_measured=available(lp),
                receiver_RNA1pct=detect(rp,view,100),receiver_RNA10pct=detect(rp,view,10),
                eligible_senders='|'.join(senders),eligible_senders10pct='|'.join(strict),candidate_qualified=bool(senders),
                assembled_complex_or_signalling_verified=False))
    edge=pd.DataFrame(edges);edge.to_csv(out/'private_all_reference_RNA_qualifications.tsv',sep='\t',index=False)
    weights=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    native_map=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    ranks=[];backgrounds=[]
    for view in views:
        candidate=Path(str(src/'private_inputs'/(view+'_discovery'))+'_genes.tsv').read_text().splitlines();assert len(candidate)==5000
        bg=[g for g in candidate if g in target.index];prior_values=target.loc[bg,ligands].to_numpy(float)
        physical=native_map.loc[bg,sources].to_numpy(int)
        for program in ['MCP1','MCP2']:
            w=weights[weights.view.eq(view)&weights.program.eq(program)].set_index('gene')
            coefficient=w.weight.reindex(bg,fill_value=0).to_numpy(float)
            backgrounds.append(dict(view=view,program=program,candidate_genes=5000,prior_background_genes=len(bg),
                original_positive_genes=int((coefficient>0).sum()),original_negative_genes=int((coefficient<0).sum())))
            for j,ligand in enumerate(ligands):
                omit=np.zeros(len(bg),bool)
                for si,source in enumerate(sources):
                    physical_exclusions={maps[source][p] for p in exclusions[ligand] if p in maps[source]}
                    omit|=np.isin(physical[:,si],list(physical_exclusions))
                keep=~omit;values=prior_values[keep,j];coef=coefficient[keep]
                native=edge[edge.receiver_view.eq(view)&edge.ligand.eq(ligand)&edge.candidate_qualified]
                for direction,sign in [('positive',1),('negative',-1)]:
                    y=(sign*coef>0).astype(int);n=int(y.sum());status='AVAILABLE'
                    if len(y)<200 or n<10 or n==len(y) or np.ptp(values)==0:status='INSUFFICIENT_EFFECTIVE_PRIOR_OR_SIGNED_SET'
                    ap=float(average_precision_score(y,values)) if status=='AVAILABLE' else np.nan
                    pearson=float(np.corrcoef(sign*coef,values)[0,1]) if status=='AVAILABLE' and np.ptp(coef)>0 else np.nan
                    spearman=float(spearmanr(sign*coef,values).statistic) if status=='AVAILABLE' and np.ptp(coef)>0 else np.nan
                    ranks.append(dict(view=view,program=program,direction=direction,ligand=ligand,
                        background_genes=len(y),signed_genes=n,excluded_background_genes=int(omit.sum()),status=status,
                        custom_AP=ap,prevalence=n/len(y) if len(y) else np.nan,AP_excess=ap-n/len(y) if len(y) else np.nan,
                        AUROC=float(roc_auc_score(y,values)) if status=='AVAILABLE' else np.nan,
                        coefficient_Pearson=pearson,coefficient_Spearman=spearman,
                        qualified_receptor_edges=len(native),qualified_senders='|'.join(sorted(set(s for val in native.eligible_senders for s in val.split('|')))),
                        native_reference_qualified=bool(len(native)),not_official_NicheNet_activity=True))
    ranking=pd.DataFrame(ranks);assert len(ranking)==14712
    ranking=ranking.sort_values(['view','program','direction','AP_excess','ligand'],ascending=[True,True,True,False,True],kind='stable')
    ranking['all_prior_rank']=ranking.groupby(['view','program','direction']).cumcount()+1
    selected=[]
    for key,group in ranking.groupby(['view','program','direction'],sort=True):
        valid=group[group.status.eq('AVAILABLE')&group.native_reference_qualified&group.AP_excess.gt(0)].head(10).copy()
        valid['fixed_qualified_rank']=np.arange(1,len(valid)+1);selected.append(valid)
    shortlist=pd.concat(selected,ignore_index=True);assert len(shortlist)<=120
    ranking.to_csv(out/'private_all14712_automatic_prior_rankings.tsv',sep='\t',index=False)
    shortlist.to_csv(out/'private_frozen_top10_per_programme_view_direction.tsv',sep='\t',index=False)
    pd.DataFrame(backgrounds).to_csv(out/'private_all6_receiver_prior_backgrounds.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    assert sum(p['bytes'] for p in files.values())<64*1024**2
    (out/'private_terminal.json').write_text(json.dumps(dict(
        status='COMPLETE_ALL_PRIOR_NATIVE_RNA_AND_AUTOMATIC_SHORTLIST_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        native_cells=19018,case_view_units=198,prior_ligands=1226,ranking_positions=14712,
        frozen_shortlist_positions=len(shortlist),source_effects_or_clinical_models_fitted=False,
        no_identity_or_original_model_change=True),indent=2)+'\n')
    print('ALL_PRIOR_NATIVE_RNA_AND_AUTOMATIC_CANDIDATE_POSITIONS_COMPLETE_REVIEW_PENDING')


if __name__=='__main__':main()
