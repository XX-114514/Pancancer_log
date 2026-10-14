"""Complete orthogonal functional interpretation of both frozen endpoints.
All current transfer sources and HNSCC included, no new P/BH or refitting.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import anndata as ad
import h5py,numpy as np,pandas as pd,zarr
from anndata._io.specs import read_elem
from scipy.stats import rankdata


def seed(key):return 20261014+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)


def partial(x,y,cov):
    N=np.column_stack((np.ones(len(x)),*[rankdata(cov[:,i]) for i in range(cov.shape[1])]))
    if len(x)-N.shape[1]<10 or np.linalg.matrix_rank(N)<N.shape[1]:return np.nan
    ex=rankdata(x);ey=rankdata(y)
    ex-=N@np.linalg.lstsq(N,ex,rcond=None)[0];ey-=N@np.linalg.lstsq(N,ey,rcond=None)[0]
    return float(np.corrcoef(ex,ey)[0,1]) if min(ex.std(),ey.std())>1e-12 else np.nan


def score_functions(raw,aliases,excluded,reference):
    background=np.array([i for i in range(raw.shape[1]) if i not in excluded],int)
    ranks_by_column={int(c):i for i,c in enumerate(background)}
    coverage=[];sets={}
    for name,genes in reference.items():
        columns=[];denom=0;removed=0
        for gene in genes:
            ix=aliases.get(gene,set())
            if len(ix)==1 and next(iter(ix)) in excluded:removed+=1;continue
            denom+=1
            if len(ix)==1 and next(iter(ix)) in ranks_by_column:columns.append(ranks_by_column[next(iter(ix))])
        columns=np.array(sorted(set(columns)),int);fraction=len(columns)/denom if denom else 0
        eligible=fraction>=.8 and len(columns)>=10
        coverage.append(dict(reference=name,original_genes=len(genes),model_overlap_removed=removed,
          measured_remaining_genes=len(columns),measurement_fraction=fraction,eligible=eligible))
        if eligible:sets[name]=columns
    X=raw[:,background].tocsr();X.eliminate_zeros();names=list(sets);scores=np.empty((X.shape[0],len(names)))
    for i in range(X.shape[0]):
        lo,hi=X.indptr[i:i+2];rank=np.full(X.shape[1],1500,dtype=float)
        rank[X.indices[lo:hi]]=np.minimum(rankdata(-X.data[lo:hi],method='average'),1500)
        for j,name in enumerate(names):
            ix=sets[name];n=len(ix);scores[i,j]=1-(rank[ix].sum()-n*(n+1)/2)/(n*1500-n*(n+1)/2)
    assert np.isfinite(scores).all() and (scores>=-1e-10).all() and (scores<=1+1e-10).all()
    for i in sorted(set([0,len(scores)-1])) if len(scores) else []:
        x=X[i].toarray().ravel();nonzero=x[x>0]
        for j,name in enumerate(names):
            n=len(sets[name]);total=0.
            for column in sets[name]:
                v=x[column]
                r=(np.count_nonzero(nonzero>v)+(np.count_nonzero(nonzero==v)+1)/2) if v>0 else 1500
                total+=min(r,1500)
            assert abs(scores[i,j]-(1-(total-n*(n+1)/2)/(n*1500-n*(n+1)/2)))<1e-10
    return scores,names,coverage


def hnscc_inputs(root,project):
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    rec=next(r for r in json.loads((atlas/'config/source_contract_v2_private.json').read_text())['records'] if r['dataset']=='GSE286827')
    g=zarr.open_group(rec['path'],mode='r');obs=read_elem(g['obs']);var=read_elem(g[rec['var_key']])
    rows=np.load(rec['row_positions_ref']);names=pd.Index('GSE286827::'+obs.index[rows].astype(str))
    eligible=pd.read_csv(root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_eligible_libraries_clinical.tsv',sep='\t')
    chunks=[]
    for chunk in pd.read_csv(atlas/'tables/all_core_annotation_private.tsv.gz',sep='\t',index_col=0,dtype=str,keep_default_na=False,chunksize=200000):
        take=chunk[chunk.dataset.eq('GSE286827')]
        if len(take):chunks.append(take[['library_id','patient_id','timepoint','major','subcluster']])
    cells=pd.concat(chunks).loc[names]
    provisional=pd.read_csv(atlas/'tables/subcluster_catalog_private.tsv',sep='\t')
    refined=pd.read_csv(atlas/'tables/subcluster_catalog_refined_private.tsv',sep='\t')
    bridge=provisional[['subcluster','node_id']].merge(refined[['node_id','subtype_candidate','subtype_evidence_status']],on='node_id',validate='one_to_one').set_index('subcluster')
    cells['subtype']=cells.subcluster.map(bridge.subtype_candidate);cells['evidence']=cells.subcluster.map(bridge.subtype_evidence_status)
    source_order=list(read_elem(g['uns'])['source_order']);measured=np.asarray(g['varm/measured_in_source'])[:,source_order.index(rec['new_parent_dataset'])].astype(bool)
    return g,rows,cells,var,eligible,measured


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    crc=root/'CRC_fine_pair_DIALOGUE_20261014_v1';prior=crc/'private_existing_source_transfer_observed_part_v1'
    out=root/'frozen_MCP_two_end_function_review_20261014_v1';assert not out.exists();out.mkdir()
    sources=['GSE148071','GSE161529','GSE182109','GSE183904','GSE286827']
    endpoints=['Tnode9','Macnode18'];variants=['unadjusted','two_native_RNA_quality_means','RNA_quality_and_cell_count']
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      sources=sources,endpoints=endpoints,program='MCP1 only supported fixed gene programme; MCP2 remains unsupported',
      function_scope='all50 *5 sources *2 endpoints *3 variants=1500 descriptive records',
      signature_and_rank_axis_exclusion='all measured genes in both original frozen endpoint programmes',
      min_ref_measured_fraction=.8,min_measured_genes=10,maxrank=1500,ties='average',min_residual_df=10,
      cell_draw='exact same deterministic draws as previous frozen source or HNSCC projection, never outcome-based',
      HNSCC_scope='new-atlas multigene CD8 and Macro_ candidates; no accepted old172 mapping',
      other_source_scope='original fixed CD8_T_cell_NKG7_9 and Macrophage_IFI6_18',
      no_new_P_BH_or_clinical_tests=True,no_weights_or_signature_refit=True,all_unestimable_records_retained=True)
    (out/'plan_before_results.json').write_text(json.dumps(plan,indent=2)+'\n')
    refroot=root/'CRC_fine_pair_DIALOGUE_code_20261014_v1/hallmark50_public_reference_v3'
    provenance=json.loads((refroot/'reference_provenance.json').read_text());assert len(provenance['records'])==50
    reference={}
    for r in provenance['records']:
        p=refroot/(r['name']+'.json');assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
        reference[r['name']]=json.loads(p.read_text())[r['name']]['geneSymbols']
    weights={e:pd.read_csv(crc/'private_transfer_fixed_weights_v1'/(e+'_MCP1_genes.tsv'),sep='\t') for e in endpoints}
    original_scores=pd.read_csv(prior/'private_patient_scores.tsv',sep='\t')
    hnscc_scores=pd.read_csv(root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_library_frozen_scores.tsv',sep='\t')
    records=[];cover=[];private_means=[];scope=[]
    for source in sources:
        if source=='GSE286827':
            g,source_rows,cells,var,eligible,measured=hnscc_inputs(root,project)
            matrix=ad.io.sparse_dataset(g['layers/counts'])
            lookup={}
            for i,(name,row) in enumerate(var.iterrows()):
                if not measured[i] or str(name).startswith('AMBIGUOUS::'):continue
                known={str(name),str(row.get('gene_symbol',''))}
                if row.get('feature_identity_status')=='stable_Ensembl_ID':known.update(str(row.get('source_gene_aliases','')).split('|'))
                for name in known:
                    if name not in ['','nan','None']:lookup.setdefault(name,set()).add(i)
            excluded=set(np.flatnonzero(~measured|var.index.astype(str).str.startswith('AMBIGUOUS::')))
        else:
            a=ad.read_h5ad(root/'state_singlecell_inputs_20261014_v1/private_native_sources'/(source+'.h5ad'))
            cells=a.obs;matrix=a.X;source_rows=np.arange(a.n_obs)
            eligible=pd.read_csv(prior/(source+'_private_selected_libraries.tsv'),sep='\t')
            lookup={str(name):{i} for i,name in enumerate(a.var_names)};excluded=set()
        # Remove both trained signatures, including those not present in the
        # measured-score subset, wherever their native physical columns exist.
        for w in weights.values():
            for gene in w.gene:
                if len(lookup.get(gene,set()))==1:excluded.update(lookup[gene])
        for endpoint,label in zip(endpoints,['CD8_T_cell_NKG7_9','Macrophage_IFI6_18']):
            if source=='GSE286827':
                major,start=('T','CD8') if endpoint=='Tnode9' else ('Myeloid','Macro_')
                chosen=cells[cells.library_id.isin(eligible.library_id)&cells.major.eq(major)&cells.subtype.str.startswith(start)&cells.evidence.eq('multigene_subtype_candidate')]
                rng=np.random.default_rng(20261014+(0 if endpoint=='Tnode9' else 1));positions=[]
                # Consume the original draw over all pre/post eligible libraries
                # before retaining baseline, so selected baseline cells match.
                for library,group in chosen.groupby('library_id',sort=True):
                    ix=cells.index.get_indexer(group.index);draw=rng.choice(ix,min(50,len(ix)),replace=False)
                    if group.timepoint.eq('pre').all():positions.extend(draw)
                positions=np.sort(np.asarray(positions,int));md=cells.iloc[positions]
                block=md.patient_id.astype(str).to_numpy()
            else:
                chosen=cells[cells.sample_uid.isin(eligible.sample_uid)&cells.subCluster.eq(label)]
                rng=np.random.default_rng(seed(source+endpoint+'cells'));positions=[]
                for _,group in chosen.groupby('analysis_block',observed=True,sort=True):
                    ix=cells.index.get_indexer(group.index);positions.extend(rng.choice(ix,min(50,len(ix)),replace=False))
                positions=np.sort(np.asarray(positions,int));md=cells.iloc[positions];block=md.analysis_block.astype(str).to_numpy()
            if not len(positions):
                for name in reference:
                    for variant in variants:records.append(dict(source=source,endpoint=endpoint,reference=name,variant=variant,n_patients=0,rho=np.nan,LOO_valid=0,LOO_min=np.nan,LOO_max=np.nan,status='no_reviewed_scored_patient_pair'))
                continue
            raw=matrix[source_rows[positions],:].tocsr();total=np.asarray(raw.sum(1)).ravel();detected=np.asarray((raw>0).sum(1)).ravel()
            assert (total>0).all() and np.isfinite(raw.data).all() and (raw.data>=0).all()
            w=weights[endpoint];indices=np.array([next(iter(lookup[gene])) if len(lookup.get(gene,set()))==1 else -1 for gene in w.gene]);valid=indices>=0
            assert np.abs(w.weight[valid]).sum()/np.abs(w.weight).sum()>=.95
            q=np.load(crc/'private_inputs'/(endpoint+'_frozen_projection.npz'))
            quality=(np.column_stack((np.log1p(total),np.log1p(detected)))-q['quality_mean'])/q['quality_sd']
            qc=pd.read_csv(crc/'private_transfer_fixed_weights_v1'/(endpoint+'_MCP1_quality.tsv'),sep='\t').set_index('term').value
            value=np.log1p(raw[:,indices[valid]].toarray()*10000/total[:,None])
            programme=((value-w.loc[valid,'mean'].to_numpy())/w.loc[valid,'sd'].to_numpy())@w.weight.to_numpy()[valid]-qc['(Intercept)']-quality[:,0]*qc.cellQ-quality[:,1]*qc.geneQ
            frame=pd.DataFrame(dict(block=block,programme=programme,logRNA=np.log1p(total),loggenes=np.log1p(detected)))
            means=frame.groupby('block').mean();means['cell_count']=frame.groupby('block').size()
            if source=='GSE286827':expected=hnscc_scores[hnscc_scores.endpoint.eq(endpoint)&hnscc_scores.timepoint.eq('pre')].set_index('patient_id').score
            else:expected=original_scores[original_scores.source.eq(source)&original_scores.endpoint.eq(endpoint)&original_scores.program.eq('MCP1')].set_index('analysis_block').score
            assert set(expected.index)==set(means.index) and np.allclose(expected.loc[means.index],means.programme,atol=1e-9,rtol=1e-9)
            # Native group cell counts, not the capped draw size.
            if source=='GSE286827':counts=chosen[chosen.timepoint.eq('pre')].groupby('patient_id').size()
            else:counts=chosen.groupby('analysis_block',observed=True).size()
            means['native_cell_count']=counts.loc[means.index].to_numpy()
            functions,names,cov=score_functions(raw,lookup,excluded,reference)
            for row in cov:cover.append(dict(source=source,endpoint=endpoint,**row))
            f=pd.DataFrame(functions,columns=names);f['block']=block;fm=f.groupby('block').mean().loc[means.index]
            for patient,row in means.iterrows():private_means.append(dict(source=source,endpoint=endpoint,block=patient,**row.to_dict(),**fm.loc[patient].to_dict()))
            x=means.programme.to_numpy();covs=[np.empty((len(means),0)),means[['logRNA','loggenes']].to_numpy(),np.column_stack((means[['logRNA','loggenes']].to_numpy(),np.log1p(means.native_cell_count)))]
            for name in reference:
                for variant,cv in zip(variants,covs):
                    r=np.nan;loo=[];status='reference_measurement_insufficient'
                    if name in fm:
                        y=fm[name].to_numpy();r=partial(x,y,cv);status='descriptive_estimated' if np.isfinite(r) else 'rank_df_or_variation_insufficient'
                        for i in range(len(x)):
                            mask=np.arange(len(x))!=i;v=partial(x[mask],y[mask],cv[mask])
                            if np.isfinite(v):loo.append(v)
                    records.append(dict(source=source,endpoint=endpoint,reference=name,variant=variant,n_patients=len(means),rho=r,LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan,status=status))
            scope.append(dict(source=source,endpoint=endpoint,cells=len(md),patient_blocks=len(means),fixed_original_score_replay=True,scalar_UCell_parity=True))
        if source!='GSE286827':del a
    result=pd.DataFrame(records);assert len(result)==1500 and not result.duplicated(['source','endpoint','reference','variant']).any()
    result.to_csv(out/'all1500_private_function_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(cover).to_csv(out/'private_reference_measurement_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(private_means).to_csv(out/'private_patient_programme_function_means.tsv',sep='\t',index=False)
    report=dict(status='COMPLETE_TWO_END_FULL50_SOURCE_FUNCTION_DIAGNOSTICS_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),scope=scope,records=1500,original_labels_changed=False,new_expression_download=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    focus=result[result.variant.eq('two_native_RNA_quality_means')&result.reference.str.contains('INTERFERON|ALLOGRAFT|IL6|TNFA',regex=True)]
    print(focus.to_json(orient='records'))


if __name__=='__main__':main()
