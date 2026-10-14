"""All 50 Hallmarks within all macrophage subtypes, avoiding gene overlap.
Descriptive only. Frozen programme genes removed from signature AND rank axis.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import h5py,numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem
from scipy.stats import rankdata


def rho(x,y,cov):
    X=np.column_stack((np.ones(len(x)),np.apply_along_axis(rankdata,0,cov))) if cov.shape[1] else np.ones((len(x),1))
    if len(x)-X.shape[1]<10 or np.linalg.matrix_rank(X)!=X.shape[1]:return np.nan
    x=rankdata(x);y=rankdata(y)
    ex=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ey=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
    return float(np.corrcoef(ex,ey)[0,1]) if min(ex.std(),ey.std())>1e-12 else np.nan


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    inputs=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1'
    refs=root/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3'
    out=root/'HNSCC_full50_orthogonal_function_diagnostics_20261014_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      scope='all50 Hallmarks * all3 observed Macro_ candidates *2 variants=300 descriptive records',
      min_ref_measured_fraction=.8,min_ref_measured_genes=10,min_residual_df=10,
      scoring='UCell CPU average-tie scalar formula, maxrank1500; native counts; same previously frozen subtype cell samples',
      overlap='remove all measurable frozen T and Mac programme genes from every Hallmark and rank background',
      variants=['unadjusted','native_RNA_detected_cellcount_rank_residual'],leave_one_patient=True,
      no_new_clinical_or_function_P_or_BH=True,not_state_annotation_or_mechanistic_validation=True)
    (out/'plan_before_function_scores.json').write_text(json.dumps(plan,indent=2)+'\n')
    contract=json.loads((atlas/'config/source_contract_v2_private.json').read_text())
    rec=next(r for r in contract['records'] if r['dataset']=='GSE286827')
    g=h5py.File(rec['path'],'r') if str(rec['path']).endswith('.h5ad') else zarr.open_group(rec['path'],mode='r')
    native=read_elem(g['obs']);var=read_elem(g[rec['var_key']]);source_rows=np.load(rec['row_positions_ref'])
    names=pd.Index('GSE286827::'+native.index[source_rows].astype(str))
    selected=pd.read_csv(inputs/'private_selected_cells.tsv',sep='\t');assert selected.cell.is_unique
    local=names.get_indexer(selected.cell);assert (local>=0).all()
    matrix=ad.io.sparse_dataset(g[rec['counts_key']]);raw=matrix[source_rows[local],:].tocsr();raw.eliminate_zeros()
    assert np.isfinite(raw.data).all() and (raw.data>=0).all()
    order=list(read_elem(g['uns'])['source_order'])
    measured=np.asarray(g['varm/measured_in_source'])[:,order.index(rec['new_parent_dataset'])].astype(bool)
    frozen=set()
    for endpoint in ['Tnode9','Macnode18']:
        mapping=pd.read_csv(root/'HNSCC_frozen_MCP_gene_axis_20261014_v2'/(endpoint+'_private_frozen_gene_columns.tsv'),sep='\t')
        frozen.update(mapping.loc[mapping.native_column.ge(0),'native_column'].astype(int))
    keep=measured&~var.index.astype(str).str.startswith('AMBIGUOUS::')
    keep[list(frozen)]=False;background=np.flatnonzero(keep);column_rank={int(c):i for i,c in enumerate(background)}
    aliases={}
    for c,(name,row) in enumerate(var.iterrows()):
        if not measured[c] or str(name).startswith('AMBIGUOUS::'):continue
        known={str(name),str(row.get('gene_symbol',''))}
        if row.get('feature_identity_status')=='stable_Ensembl_ID':known.update(str(row.get('source_gene_aliases','')).split('|'))
        for name in known:
            if name not in ['','nan','None']:aliases.setdefault(name,set()).add(c)
    provenance=json.loads((refs/'reference_provenance.json').read_text());assert len(provenance['records'])==50
    coverage=[];signatures={}
    for r in provenance['records']:
        path=refs/(r['name']+'.json');assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
        genes=json.loads(path.read_text())[r['name']]['geneSymbols'];assert len(set(genes))==len(genes)
        ix=[];excluded=0;denom=0
        for gene in genes:
            candidates=aliases.get(gene,set())
            if len(candidates)==1 and next(iter(candidates)) in frozen:excluded+=1;continue
            denom+=1
            if len(candidates)==1 and next(iter(candidates)) in column_rank:ix.append(column_rank[next(iter(candidates))])
        ix=np.array(sorted(set(ix)),int);fraction=len(ix)/denom if denom else 0
        eligible=fraction>=.8 and len(ix)>=10
        coverage.append(dict(reference=r['name'],original_genes=len(genes),frozen_programme_genes_removed=excluded,
          retained_reference_genes=denom,measured_nonoverlapping_genes=len(ix),measured_fraction=fraction,score_measurable=eligible))
        if eligible:signatures[r['name']]=ix
    X=raw[:,background].tocsr();cell_scores=np.empty((len(selected),len(signatures)),float)
    signature_names=list(signatures)
    for i in range(X.shape[0]):
        lo,hi=X.indptr[i:i+2];indices=X.indices[lo:hi];values=X.data[lo:hi]
        ranks=np.full(X.shape[1],1500,dtype=float)
        ranks[indices]=np.minimum(rankdata(-values,method='average'),1500)
        for j,name in enumerate(signature_names):
            ix=signatures[name];n=len(ix);U=float(ranks[ix].sum()-n*(n+1)/2)
            cell_scores[i,j]=1-U/(n*1500-n*(n+1)/2)
    assert np.isfinite(cell_scores).all() and (cell_scores>=-1e-9).all() and (cell_scores<=1+1e-9).all()
    # Full sparse and dense independent order comparisons for first/last cells.
    for i in sorted(set([0,len(selected)-1])):
        values=X[i].toarray().ravel();nonzero=np.flatnonzero(values)
        for j,name in enumerate(signature_names):
            total=0.;ix=signatures[name];n=len(ix)
            for column in ix:
                v=values[column]
                rank=(np.count_nonzero(values[nonzero]>v)+(np.count_nonzero(values[nonzero]==v)+1)/2) if v>0 else 1500
                total+=min(rank,1500)
            direct=1-(total-n*(n+1)/2)/(n*1500-n*(n+1)/2)
            assert abs(direct-cell_scores[i,j])<1e-10
    # Match each selected cell to original library via canonical native rows.
    mapping=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    library_field='original_library_id'
    assert library_field in native
    selected['library_id']=native.iloc[source_rows[local]][library_field].astype(str).to_numpy()
    assert selected.library_id.isin(mapping.library_id).all()
    values=pd.DataFrame(cell_scores,columns=signature_names)
    values['library_id']=selected.library_id.to_numpy();values['subtype_candidate']=selected.subtype_candidate.to_numpy()
    means=values.groupby(['library_id','subtype_candidate']).mean()
    programme=pd.read_csv(inputs/'private_patient_subtype_scores.tsv',sep='\t').set_index(['library_id','subtype_candidate'])
    assert set(means.index)==set(programme.index)
    records=[]
    for subtype in sorted(programme.index.get_level_values('subtype_candidate').unique()):
        md=programme.xs(subtype,level='subtype_candidate').sort_index();functions=means.xs(subtype,level='subtype_candidate').loc[md.index]
        cov=np.column_stack((md.mean_logRNA,md.mean_log_genes,np.log1p(md.native_subtype_cells)))
        x=md.score.to_numpy()
        for coverage_row in coverage:
            name=coverage_row['reference']
            for variant,cv in [('unadjusted',np.empty((len(md),0))),('native_RNA_detected_cellcount_rank_residual',cov)]:
                row=dict(subtype_candidate=subtype,reference=name,variant=variant,n_patients=len(md),rho=np.nan,
                  LOO_valid=0,LOO_min=np.nan,LOO_max=np.nan,status='reference_measurement_insufficient')
                if name in signatures:
                    y=functions[name].to_numpy();r=rho(x,y,cv);loo=[]
                    for i in range(len(x)):
                        mask=np.arange(len(x))!=i;v=rho(x[mask],y[mask],cv[mask])
                        if np.isfinite(v):loo.append(v)
                    row.update(rho=r,status='descriptive_estimated' if np.isfinite(r) else 'rank_or_variation_insufficient',
                      LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan)
                records.append(row)
    result=pd.DataFrame(records);assert len(result)==300
    result.to_csv(out/'all300_private_descriptive_functions.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(coverage).to_csv(out/'all50_reference_coverage.tsv',sep='\t',index=False)
    means.to_csv(out/'private_patient_subtype_function_means.tsv',sep='\t')
    report=dict(status='ALL50_NONOVERLAPPING_FUNCTION_DIAGNOSTICS_COMPLETE_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),
      references=50,measurable_references=len(signatures),descriptive_rows=300,original_labels_changed=False,
      trained_genes_removed_from_rank_axis=True,scalar_UCell_parity=True)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
    # Pre-specified functional categories from prior complete CRC diagnostics,
    # followed by top absolute descriptive correlations, never a P-value filter.
    focus=result[result.subtype_candidate.eq('Macro_C1QC')&result.variant.eq('native_RNA_detected_cellcount_rank_residual')]
    summaries=focus[focus.reference.str.contains('INTERFERON|TNFA|IL6|TGF_BETA|PROTEIN_SECRETION',regex=True)]
    print(summaries.to_json(orient='records'))
    print(focus.assign(absolute_rho=focus.rho.abs()).sort_values('absolute_rho',ascending=False).head(8).drop(columns='absolute_rho').to_json(orient='records'))


if __name__=='__main__':main()
