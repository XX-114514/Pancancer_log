"""Independent raw native assay replay and complete automatic ranking audit."""
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata.io import read_elem, sparse_dataset

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import PRIORS, read_rds


def close(x,y):
    assert (pd.isna(x) and pd.isna(y)) or abs(x-y)<1e-9


def precision_metrics(labels,values):
    """Explicit tied-score AP and Mann-Whitney AUC; no sklearn producer code."""
    n=int(labels.sum());order=np.argsort(-values,kind='stable')
    score=values[order];positive=np.cumsum(labels[order])
    ends=np.r_[np.flatnonzero(score[1:]!=score[:-1]),len(score)-1]
    at=positive[ends];new=np.diff(np.r_[0,at])
    ap=float(np.sum(new*at/(ends+1))/n)
    ranks=pd.Series(values).rank().to_numpy()
    auc=float((ranks[labels.astype(bool)].sum()-n*(n+1)/2)/(n*(len(labels)-n)))
    return ap,auc


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    source=root/'immune_RNA_MCP_all_ligand_qualification_20261018_v2'
    terminal=json.loads((source/'private_terminal.json').read_text())
    assert terminal['source_sha256']=='153213c6c03ef0bbf5f7e579f52615f079c9874c275ac2e53b65baa6100370d5'
    assert terminal['ranking_positions']==14712
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    out=root/'immune_RNA_MCP_all_ligand_qualification_audit_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_native_and_ranking_audit.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        producer_terminal_sha256=sha(source/'private_terminal.json'),
        scope='All19018 cells/all198 case-view native ligand/receptor measurements using alternate unitwise CSR access; independently reconstructed physical gene maps, all prior/reference qualifications, all14712 metrics via explicit tied-score AP/AUC and entire top10 family. No candidate/weight output.'),indent=2)+'\n')
    assays=np.load(source/'private_all198_native_ligand_receptor_assays.npz',allow_pickle=False)
    genes=assays['genes'].tolist();index={g:i for i,g in enumerate(genes)};requested=set(genes)
    units=pd.read_csv(source/'private_all198_original_case_view_units.tsv',sep='\t')
    stored=pd.read_csv(source/'private_all_prior_native_gene_columns.tsv',sep='\t').set_index('gene')
    md=pd.read_csv(root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_selected_native_cell_depth.tsv',sep='\t')
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r')
    order=list(np.asarray(read_elem(z['uns/source_order'])).astype(str));measured=np.asarray(read_elem(z['varm/measured_in_source']),bool)
    symbol=np.asarray(read_elem(z['var/gene_symbol'])).astype(str);canonical=np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    identity=np.asarray(read_elem(z['var/feature_identity_status'])).astype(str);aliases=np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    sources=sorted(units.dataset.unique());maps={};columns={}
    for dataset in sources:
        possible=defaultdict(set)
        for j in np.flatnonzero(measured[:,order.index(dataset)]):
            if identity[j]=='ambiguous_source_feature' or canonical[j].startswith('AMBIGUOUS::'):continue
            labels={symbol[j],canonical[j]}
            if identity[j]=='stable_Ensembl_ID':labels|=set(aliases[j].split('|'))
            for g in labels&requested:possible[g].add(int(j))
        maps[dataset]={g:next(iter(cols)) for g,cols in possible.items() if len(cols)==1}
        columns[dataset]=np.array([maps[dataset].get(g,-1) for g in genes],int)
        assert np.array_equal(stored.loc[genes,dataset],columns[dataset])
    counts=sparse_dataset(z['layers/counts']);hits=np.zeros_like(assays['positive_cells']);means=np.full_like(assays['mean_logCP10K'],np.nan)
    read=0;maximum_native_error=0.
    for row,r in enumerate(units.itertuples(index=False)):
        cells=md[md.fit_sample.eq(r.fit_sample)&md.view.eq(r.view)].sort_values('global_row')
        assert len(cells)==r.selected_cells==assays['cells'][row]
        x=counts[cells.global_row.to_numpy(int),:].tocsr();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
        assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(),cells.fullRNA.to_numpy())
        assert np.array_equal(np.diff(x.indptr),cells.full_detected_genes.to_numpy())
        valid=columns[r.dataset]>=0;v=x[:,columns[r.dataset][valid]].astype(float)
        hits[row,valid]=np.asarray((v>0).sum(axis=0)).ravel()
        # Independent row-by-row normalization, not producer repeated sparse factors.
        for k,total in enumerate(cells.fullRNA.to_numpy(float)):
            a,b=v.indptr[k:k+2];v.data[a:b]=np.log1p(v.data[a:b]*10000/total)
        means[row,valid]=np.asarray(v.sum(axis=0)).ravel()/len(cells)
        assert np.array_equal(hits[row],assays['positive_cells'][row])
        error=float(np.nanmax(np.abs(means[row]-assays['mean_logCP10K'][row])))
        maximum_native_error=max(maximum_native_error,error);assert error<1e-9
        assert np.array_equal(np.isnan(means[row]),np.isnan(assays['mean_logCP10K'][row]))
        read+=len(cells)
    assert read==19018
    views=sorted(units.view.unique());pools={};totals={}
    for view in views:
        take=units.view.eq(view)&units.partition.eq('reviewed_case_discovery')
        pools[view]=hits[take].sum(axis=0);totals[view]=int(assays['cells'][take].sum())
    def parts(g):return str(g).split('_')
    def measurable(names):return all(g in index and all(columns[s][index[g]]>=0 for s in sources) for g in names)
    def detectable(names,view,threshold):
        return measurable(names) and all(pools[view][index[g]]*threshold>=totals[view] for g in names)
    edge=pd.read_csv(source/'private_all_reference_RNA_qualifications.tsv',sep='\t',keep_default_na=False)
    accepted={'simplex_reference_supported','all_constituents_supported_complex_assembly_unverified'}
    for r in edge.itertuples(index=False):
        rp,lp=parts(r.receptor),parts(r.ligand)
        send=[v for v in views if v!=r.receiver_view and detectable(lp,v,100)] if r.reference_status in accepted and detectable(rp,r.receiver_view,100) else []
        strict=[v for v in send if detectable(rp,r.receiver_view,10) and detectable(lp,v,10)]
        assert r.eligible_senders=='|'.join(send) and r.eligible_senders10pct=='|'.join(strict)
        assert r.candidate_qualified==bool(send) and r.receptor_components_measured==measurable(rp) and r.ligand_components_measured==measurable(lp)
    prior=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    for name,digest in PRIORS.items():assert sha(prior/name)==digest
    target=read_rds(prior/'ligand_target_matrix_nsga2r_final.rds');target.index=target.index.astype(str);target.columns=target.columns.astype(str)
    flat=read_rds(prior/'lr_network_human_21122021.rds');flat['from']=flat['from'].astype(str);flat['to']=flat['to'].astype(str)
    ref=pd.read_csv(root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv',sep='\t',keep_default_na=False)
    expected_ref=set(map(tuple,ref[ref.ligand.isin(target.columns)][['ligand','receptor','reference_status']].to_numpy()))
    for view in views:
        assert set(map(tuple,edge[edge.receiver_view.eq(view)][['ligand','receptor','reference_status']].to_numpy()))==expected_ref
    ranking=pd.read_csv(source/'private_all14712_automatic_prior_rankings.tsv',sep='\t',keep_default_na=False)
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1'
    assert sha(projection/'private_terminal.json')=='06c5637d9acf562a9266436da9f814af403fceed330a6bfd951ff45dc1d920a3'
    for name,digest in json.loads((projection/'private_terminal.json').read_text())['files'].items():assert sha(projection/name)==digest
    weights=pd.read_csv(projection/'private_R_projection/private_frozen_gene_weights.tsv',sep='\t')
    weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    native=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    checked=0
    for view in views:
        names=(root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'/(view+'_discovery_genes.tsv')).read_text().splitlines()
        bg=[g for g in names if g in target.index]
        native_columns=native.loc[bg,sources].to_numpy(int)
        for program in ['MCP1','MCP2']:
            w=weights[weights.view.eq(view)&weights.program.eq(program)].set_index('gene').weight.reindex(bg,fill_value=0)
            group=ranking[ranking.view.eq(view)&ranking.program.eq(program)].set_index(['ligand','direction'])
            for ligand in sorted(target.columns):
                strings={ligand}|set(flat.loc[flat['from'].eq(ligand),'to'])|set(ref.loc[ref.ligand.eq(ligand),'receptor'])
                exclude=set(g for text in strings for g in parts(text))
                physical={s:{maps[s][g] for g in exclude if g in maps[s]} for s in sources}
                keep=np.ones(len(bg),bool)
                for si,s in enumerate(sources):keep&=~np.isin(native_columns[:,si],list(physical[s]))
                vals=target.loc[np.array(bg)[keep],ligand].to_numpy(float);coef=w.to_numpy()[keep]
                e=edge[edge.receiver_view.eq(view)&edge.ligand.eq(ligand)&edge.candidate_qualified]
                for direction,sign in [('positive',1),('negative',-1)]:
                    r=group.loc[(ligand,direction)];assert isinstance(r,pd.Series)
                    y=(sign*coef>0).astype(int);n=int(y.sum());available=len(y)>=200 and n>=10 and n<len(y) and np.ptp(vals)>0
                    assert r.status==('AVAILABLE' if available else 'INSUFFICIENT_EFFECTIVE_PRIOR_OR_SIGNED_SET')
                    assert r.background_genes==len(y) and r.signed_genes==n and r.excluded_background_genes==int((~keep).sum())
                    assert r.native_reference_qualified==bool(len(e)) and r.qualified_receptor_edges==len(e)
                    if available:
                        ap,auc=precision_metrics(y,vals);close(float(r.custom_AP),ap);close(float(r.AUROC),auc);close(float(r.AP_excess),ap-n/len(y))
                        a=sign*coef-np.mean(sign*coef);b=vals-vals.mean();close(float(r.coefficient_Pearson),float(a@b/np.sqrt((a@a)*(b@b))))
                        a=pd.Series(sign*coef).rank().to_numpy();b=pd.Series(vals).rank().to_numpy();a-=a.mean();b-=b.mean()
                        close(float(r.coefficient_Spearman),float(a@b/np.sqrt((a@a)*(b@b))))
                    checked+=1
    assert checked==14712
    # Convert only numerical missing fields, keeping all candidate identifiers private.
    for name in ['AP_excess','custom_AP']:ranking[name]=pd.to_numeric(ranking[name],errors='coerce')
    expected=[]
    for _,g in ranking.groupby(['view','program','direction'],sort=True):
        g=g.sort_values(['AP_excess','ligand'],ascending=[False,True],kind='stable')
        assert g.all_prior_rank.tolist()==list(range(1,1227))
        g=g[g.status.eq('AVAILABLE')&g.native_reference_qualified&g.AP_excess.gt(0)].head(10)
        expected.extend(map(tuple,g[['view','program','direction','ligand']].to_numpy()))
    selected=pd.read_csv(source/'private_frozen_top10_per_programme_view_direction.tsv',sep='\t')
    assert expected==list(map(tuple,selected[['view','program','direction','ligand']].to_numpy()))
    (out/'private_terminal.json').write_text(json.dumps(dict(
        status='ACCEPTED_ALL_NATIVE_ASSAYS_REFERENCE_QUALIFICATIONS_AND14712_AUTOMATIC_RANKINGS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        producer_terminal_sha256=sha(source/'private_terminal.json'),native_cells_replayed=read,
        native_case_view_units=198,prior_ranking_positions=checked,reference_qualifications=len(edge),
        frozen_shortlist_positions=len(selected),maximum_native_mean_error=maximum_native_error,
        clinical_or_directional_effect_model=False),indent=2)+'\n')
    print('ALL_NATIVE_ASSAYS_AND_AUTOMATIC_RANKINGS_INDEPENDENTLY_ACCEPTED')


if __name__=='__main__':main()
