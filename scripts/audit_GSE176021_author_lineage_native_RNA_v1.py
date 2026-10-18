"""Full independent native replay: scalar fixed models and every gene/cell summary."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def close(x,y,tol=1e-9):assert np.allclose(x,y,rtol=0,atol=tol,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);src=root/'GSE176021_author_lineage_RNA_scope_20261019_v2';producer=root/'GSE176021_author_lineage_native_RNA_20261019_v1';out=root/'GSE176021_author_lineage_native_RNA_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    t=wait(root,out,'measure_GSE176021_author_lineage_native_RNA_v1',producer/'private_terminal.json',{'COMPLETE_ALL_AUTHOR_SCOPE_NATIVE_RNA_AND824_ANCHORS_INDEPENDENT_REPLAY_PENDING'})
    for name,v in t['files'].items():assert sha(producer/name)==v['sha256']
    audit=root/'GSE176021_author_lineage_RNA_scope_audit_20261019_v1/private_terminal.json';assert sha(audit)=='fe5d399430fb1ced7637bb17200ef2327c794ae808d4653a77cbd1238504b982'
    pools=pd.read_csv(src/'private_all_eligible_pools.tsv',sep='\t');membership=np.load(src/'private_all_scope_cell_memberships.npz');mr=membership['global_row'].astype(np.int64);pi=membership['pool_index'];selected=membership['selected'];rows=np.unique(mr);n=len(rows);nu=len(pools);gm=pd.read_csv(src/'private_all904_native_gene_roles_and_mapping.tsv',sep='\t');genes=gm.gene.tolist();physical=np.unique(gm.loc[gm.native_column.ge(0),'native_column'].to_numpy(int));columns=np.full(904,-1,int);valid=np.flatnonzero(gm.native_column.to_numpy(int)>=0);columns[valid]=np.searchsorted(physical,gm.native_column.to_numpy(int)[valid]);stats=np.load(producer/'private_all_gene_pool_sufficient_statistics.npz');assert stats['genes'].tolist()==genes
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';params=pd.read_csv(native/'private_all12_frozen_models.tsv',sep='\t').set_index(['view','program','variant']);weights=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv',sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')];stress={'FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'};models={};inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    for vi,view in enumerate(['R_FOXP3_capture','E_effector_capture']):
        frozen=np.load(inputs/(view+'_frozen_projection.npz'));pieces={}
        for pr,program in enumerate(['MCP1','MCP2']):
            w=weights[weights.view.eq(view)&weights.program.eq(program)]
            for st,variant in enumerate(['full','fixed_stress_complement']):
                use=w if st==0 else w[~w.gene.str.upper().isin(stress)];p=params.loc[(view,program,variant)];ix=pd.Index(genes).get_indexer(use.gene);assert (ix>=0).all() and (columns[ix]>=0).all();pieces[pr*2+st]={'columns':columns[ix],'b':use.weight.to_numpy()/use.discovery_sd.to_numpy(),'center':float(p.center),'coef':p[['intercept','cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)}
        models[vi]={'mu':frozen['quality_mean'],'sd':frozen['quality_sd'],'pieces':pieces,'training':pd.read_csv(inputs/(view+'_discovery_metadata.tsv'),sep='\t')}
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');raw=z['layers/counts'];nf=int(raw.attrs['shape'][1]);allraw=np.zeros((nu,904),np.int64);alldet=allraw.copy();alllog=np.zeros((nu,904));selraw=allraw.copy();seldet=allraw.copy();sellog=alllog.copy();RNA=np.empty(n,np.int64);det=np.empty(n,np.int64);qv=np.empty((len(mr),4));sv=np.empty((len(mr),4));uv=np.empty_like(sv);nnz=0;hashes={k:hashlib.sha256() for k in ['native_data_int64','native_indices_int64','native_row_endpoints_int64','global_row_int64','RNA_int64','detected_int64']}
    for lo in range(0,n,512):
        hi=min(n,lo+512);rr=rows[lo:hi];first=int(rr[0]);last=int(rr[-1])+1;ptr=np.asarray(raw['indptr'][first:last+1],np.int64);x=csr_matrix((np.asarray(raw['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(raw['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(last-first,nf))[rr-first].copy();x.sum_duplicates();x.eliminate_zeros();x.sort_indices();assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all();x=x.astype(np.int64);RNA[lo:hi]=np.asarray(x.sum(axis=1)).ravel();det[lo:hi]=np.diff(x.indptr);assert (RNA[lo:hi]>0).all()
        for k,v in [('native_data_int64',x.data),('native_indices_int64',x.indices),('native_row_endpoints_int64',x.indptr[1:].astype(np.int64)+nnz),('global_row_int64',rr),('RNA_int64',RNA[lo:hi]),('detected_int64',det[lo:hi])]:hashes[k].update(np.asarray(v,dtype='<i8').tobytes())
        nnz+=x.nnz;dense=x[:,physical].toarray();logged=np.log1p(dense.astype(float)/RNA[lo:hi,None]*1e4);a,b=np.searchsorted(mr,[first,last]);positions=np.searchsorted(rr,mr[a:b])
        for pool in set(pi[a:b].tolist()):
            local=np.flatnonzero(pi[a:b]==pool);mids=a+local;xi=positions[local];p=pools.iloc[pool];md=models[int(p.view_index)];d=dense[xi][:,columns[valid]];l=logged[xi][:,columns[valid]];chosen=selected[mids]
            allraw[pool,valid]+=np.add.reduce(d,axis=0);alldet[pool,valid]+=np.count_nonzero(d,axis=0);alllog[pool,valid]+=np.add.reduce(l,axis=0);selraw[pool,valid]+=np.add.reduce(d[chosen],axis=0);seldet[pool,valid]+=np.count_nonzero(d[chosen],axis=0);sellog[pool,valid]+=np.add.reduce(l[chosen],axis=0)
            q=np.column_stack([np.log1p(RNA[lo+xi]),np.log1p(det[lo+xi]),np.repeat(p.original_capture_cells/p.annotated_scope_cells,len(xi)),np.repeat(np.log1p(p.original_capture_cells),len(xi))]);q=(q-md['mu'])/md['sd'];qv[mids]=q
            for k,m in md['pieces'].items():
                signal=np.sum(logged[xi][:,m['columns']]*m['b'],axis=1)-m['center'];c=m['coef'];nuisance=c[0]+q[:,0]*c[1]+q[:,1]*c[2]+q[:,2]*c[4]+q[:,3]*c[5];uv[mids,k]=signal;sv[mids,k]=signal-nuisance
        if lo//512%32==0 or hi==n:(out/'private_progress.json').write_text(json.dumps({'phase':'INDEPENDENT_ALL_NATIVE_AUTHOR_SCOPE_REPLAY','unique_cells_replayed':hi,'unique_cells':n,'observed_UTC':datetime.now(timezone.utc).isoformat()})+'\n')
    dig=json.loads((producer/'private_all_unique_native_canonical_digests.json').read_text());assert dig['canonical_nonzero_values']==nnz and dig['digests']=={k:v.hexdigest() for k,v in hashes.items()}
    for arr in [allraw,alldet,selraw,seldet]:arr[:,gm.native_column.to_numpy(int)<0]=-1
    for arr in [alllog,sellog]:arr[:,gm.native_column.to_numpy(int)<0]=np.nan
    for key,arr in [('all_raw_sum',allraw),('all_detected_cells',alldet),('selected_raw_sum',selraw),('selected_detected_cells',seldet)]:assert np.array_equal(stats[key],arr)
    for key,arr,N in [('all_logCP10K_sum',alllog,pools.cells.to_numpy()),('selected_logCP10K_sum',sellog,pools.selected_cells.to_numpy())]:close(stats[key]/N[:,None],arr/N[:,None],1e-10)
    saved=np.load(producer/'private_selected_cell_quality_scores.npz');uix=np.searchsorted(rows,mr[selected]);assert np.array_equal(saved['global_row'],mr[selected]) and np.array_equal(saved['pool_index'],pi[selected]) and np.array_equal(saved['fullRNA'],RNA[uix]) and np.array_equal(saved['detected'],det[uix]);close(saved['quality'],qv[selected]);close(saved['scores'],sv[selected]);maxerr=float(np.max(np.abs(saved['scores']-sv[selected])))
    table=pd.read_csv(producer/'private_all_author_scope_programme_measurements.tsv',sep='\t');assert len(table)==len(pools)*4;idx=table.set_index(['pool_index','program','variant']);domain=pd.read_csv(producer/'private_all_training_domain_positions.tsv',sep='\t').set_index(['pool_index','covariate']);assert len(domain)==len(pools)*4
    for pool,p in pools.iterrows():
        ids=np.flatnonzero(pi==pool);ss=ids[selected[ids]];md=models[int(p.view_index)];own=(np.array([p.cells/p.annotated_scope_cells,np.log1p(p.cells)])-md['mu'][2:])/md['sd'][2:]
        for j,name in enumerate(['cellQ','geneQ','fraction','logN']):
            rr=domain.loc[(pool,name)];tr=md['training'][name];v=qv[ids,j];close([rr['min'],rr['max'],rr.training_min,rr.training_max],[v.min(),v.max(),tr.min(),tr.max()]);assert rr.outside_training_cells==int(((v<tr.min())|(v>tr.max())).sum()) and rr.cells==len(ids)
        for k,m in md['pieces'].items():
            pr='MCP'+str(k//2+1);variant='full' if k%2==0 else 'fixed_stress_complement';r=idx.loc[(pool,pr,variant)];delta=-np.sum((own-qv[ids[0],2:])*m['coef'][4:]);close([r.score,r.all_captured_mean_sensitivity,r.own_population_background_score,r.all_own_population_background_score,r.unresidualized_weighted_RNA_score,r.RNA_quality_only_background_score],[sv[ss,k].mean(),sv[ids,k].mean(),sv[ss,k].mean()+delta,sv[ids,k].mean()+delta,uv[ss,k].mean(),sv[ss,k].mean()+np.sum(qv[ids[0],2:]*m['coef'][4:])]);assert r['scope']==p.scope and r.sampleID==p.sampleID and r['view']==p['view'] and r.full_view_cells==len(ids) and r.selected_cells==len(ss)
            close(r[['cellQ','geneQ','fraction','logN']].to_numpy(float),qv[ss].mean(axis=0));close(r[['all_cellQ','all_geneQ','all_fraction','all_logN']].to_numpy(float),qv[ids].mean(axis=0));close([r.population_fraction,r.population_logN],[p.cells/p.annotated_scope_cells,np.log1p(p.cells)])
    old=pd.read_csv(native/'private_all3488_programme_measurements.tsv',sep='\t').set_index(['sampleID','view','program','variant']);anchors=table[table.scope.eq('legacy_capture_anchor')];assert len(anchors)==824
    for r in anchors.itertuples(index=False):
        x=old.loc[(r.sampleID,r.view,r.program,r.variant)];close([r.score,r.all_captured_mean_sensitivity],[x.score,x.all_captured_mean_sensitivity])
    receipt={'status':'ACCEPTED_COMPLETE_ALL_NATIVE_AUTHOR_SCOPES_SCORES_GENE_STATS_DOMAINS_AND824_ANCHORS','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(producer/'private_terminal.json'),'scope_acceptance_sha256':sha(audit),'unique_native_cells':n,'all_memberships':len(mr),'selected_memberships':int(selected.sum()),'canonical_native_nonzero_values':int(nnz),'gene_pool_positions':len(pools)*904,'programme_measurements':len(table),'legacy_anchors':824,'maximum_selected_score_error':maxerr,'new_association_or_clinical_effect':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json'),'maximum_error':maxerr}))


if __name__=='__main__':main()
