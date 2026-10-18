"""Full native RNA for all frozen author-lineage scopes, exact old anchors."""
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

STRESS={'FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'}
VIEWS=['R_FOXP3_capture','E_effector_capture']
QFIELDS=['cellQ','geneQ','fraction','logN']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,64*1024**2);src=root/'GSE176021_author_lineage_RNA_scope_20261019_v2';audit=root/'GSE176021_author_lineage_RNA_scope_audit_20261019_v1/private_terminal.json';assert sha(audit)=='fe5d399430fb1ced7637bb17200ef2327c794ae808d4653a77cbd1238504b982'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,v in terminal['files'].items():assert sha(src/name)==v['sha256']
    assert terminal['conservative_native_uncompressed_bound_bytes']<64*1024**2
    out=root/'GSE176021_author_lineage_native_RNA_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'scope_acceptance_sha256':sha(audit),'resource_guard':guard,'full_scope':'All549 frozen pools,all401966 unique original cells/all741561 memberships,all904 predefined genes. Original scopes/full-stress/both programmes and original frozen transformations. No new weights or expression download.','background':'Primary original capture fraction/logN retained. Own-population composition adjustment and RNA-only/unresidualized weighted signal stored as separately labelled sensitivities. No silent change in learned quality meaning.','paired_protocol_addendum':'All6 main questions already frozen on10 primary shared pairs. Author-confirmed capture has9 complete pairs; full10-case sensitivity marked unavailable, separately report matched9-case C/B point-only comparison, no new P. Coactivity point/LOO only for all scope/tissue/program contexts; no new significance screen.','storage':'Full per-gene sufficient sums/detections for all and hash100 members; selected individual qualities/scores plus full native canonical digests. Raw-count matrices never saved. Independent full-native replay required before effect review.'}
    (out/'plan_before_all_new_author_RNA.json').write_text(json.dumps(plan,indent=2)+'\n')
    pools=pd.read_csv(src/'private_all_eligible_pools.tsv',sep='\t');lib=pd.read_csv(src/'private_all_source_materials.tsv',sep='\t',keep_default_na=False);members=np.load(src/'private_all_scope_cell_memberships.npz');mr=members['global_row'].astype(np.int64);pi=members['pool_index'];selected=members['selected'];rows=np.unique(mr);npools=len(pools);n=len(rows);ng=904
    gm=pd.read_csv(src/'private_all904_native_gene_roles_and_mapping.tsv',sep='\t');genes=gm.gene.tolist();col=gm.native_column.to_numpy(int);valid=np.flatnonzero(col>=0);physical=col[valid];assert len(genes)==ng and len(set(physical))==len(physical)
    native=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';params=pd.read_csv(native/'private_all12_frozen_models.tsv',sep='\t').set_index(['view','program','variant']);wpath=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wpath)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7';weights=pd.read_csv(wpath,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')];models={};records=[]
    train=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'
    for vi,view in enumerate(VIEWS):
        frozen=np.load(train/(view+'_frozen_projection.npz'));B=np.zeros((ng,4));centers=[];coefs=[]
        for pr,program in enumerate(['MCP1','MCP2']):
            w=weights[weights.view.eq(view)&weights.program.eq(program)]
            for st,variant in enumerate(['full','fixed_stress_complement']):
                v=w if st==0 else w[~w.gene.str.upper().isin(STRESS)];b=v.weight.to_numpy()/v.discovery_sd.to_numpy();ii=pd.Index(genes).get_indexer(v.gene);assert (ii>=0).all() and (col[ii]>=0).all();B[ii,2*pr+st]=b;record=params.loc[(view,program,variant)];centers.append(float(record.center));coefs.append(record[['intercept','cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float));assert abs(float(b@v.discovery_mean.to_numpy())-record.center)<1e-10
                records.append({'view':view,'program':program,'variant':variant,'active_genes':len(v),'center':record.center})
        models[vi]={'B':B,'center':np.array(centers),'coef':np.array(coefs),'mu':frozen['quality_mean'],'sd':frozen['quality_sd'],'training':pd.read_csv(train/(view+'_discovery_metadata.tsv'),sep='\t')}
    z=zarr.open_group(str(root.parents[1]/'runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'),mode='r');nativecounts=z['layers/counts'];nf=int(nativecounts.attrs['shape'][1]);allraw=np.zeros((npools,ng),np.int64);alldet=np.zeros_like(allraw);alllog=np.zeros((npools,ng));selraw=allraw.copy();seldet=alldet.copy();sellog=alllog.copy();RNA=np.zeros(n,np.int64);det=np.zeros(n,np.int32);quality=np.empty((len(mr),4));scores=np.empty((len(mr),4));unres=np.empty_like(scores)
    hashes={k:hashlib.sha256() for k in ['native_data_int64','native_indices_int64','native_row_endpoints_int64','global_row_int64','RNA_int64','detected_int64']};nnz=0
    for lo in range(0,n,2048):
        hi=min(n,lo+2048);rr=rows[lo:hi];start=int(rr[0]);stop=int(rr[-1])+1;ptr=np.asarray(nativecounts['indptr'][start:stop+1],np.int64);x=csr_matrix((np.asarray(nativecounts['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(nativecounts['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,nf))[rr-start].copy();x.sum_duplicates();x.eliminate_zeros();x.sort_indices();assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.rint(x.data)).all();x=x.astype(np.int64);RNA[lo:hi]=np.asarray(x.sum(axis=1)).ravel();det[lo:hi]=np.diff(x.indptr);assert (RNA[lo:hi]>0).all()
        for k,v in [('native_data_int64',x.data),('native_indices_int64',x.indices),('native_row_endpoints_int64',x.indptr[1:]+nnz),('global_row_int64',rr),('RNA_int64',RNA[lo:hi]),('detected_int64',det[lo:hi])]:hashes[k].update(np.asarray(v,dtype='<i8').tobytes())
        nnz+=x.nnz;dense=x[:,physical].toarray();log=np.log1p(dense.astype(float)*10000/RNA[lo:hi,None]);ml,mh=np.searchsorted(mr,[int(rr[0]),int(rr[-1])+1]);poolids=pi[ml:mh];where=np.searchsorted(rr,mr[ml:mh])
        for pool in np.unique(poolids):
            local=np.flatnonzero(poolids==pool);mid=ml+local;xi=where[local];take=selected[mid];p=pools.iloc[int(pool)];model=models[int(p.view_index)];d=dense[xi];l=log[xi]
            allraw[pool,valid]+=d.sum(axis=0);alldet[pool,valid]+=(d>0).sum(axis=0);alllog[pool,valid]+=l.sum(axis=0);selraw[pool,valid]+=d[take].sum(axis=0);seldet[pool,valid]+=(d[take]>0).sum(axis=0);sellog[pool,valid]+=l[take].sum(axis=0)
            q=np.c_[np.log1p(RNA[lo+xi]),np.log1p(det[lo+xi]),np.full(len(xi),p.original_capture_cells/p.annotated_scope_cells),np.full(len(xi),np.log1p(p.original_capture_cells))];q=(q-model['mu'])/model['sd'];quality[mid]=q;D=np.c_[np.ones(len(xi)),q[:,:2],np.zeros(len(xi)),q[:,2:]];unres[mid]=l@model['B'][valid]-model['center'];scores[mid]=unres[mid]-D@model['coef'].T
        if lo//2048%16==0 or hi==n:(out/'private_progress.json').write_text(json.dumps({'phase':'ALL_UNIQUE_NATIVE_AUTHOR_AND_ANCHOR_CELLS','unique_cells_read':hi,'unique_cells':n,'observed_UTC':datetime.now(timezone.utc).isoformat()})+'\n')
    assert np.isfinite(quality).all() and np.isfinite(scores).all()
    records=[];domain=[]
    for pool,p in pools.iterrows():
        ids=np.flatnonzero(pi==pool);chosen=ids[selected[ids]];md=models[int(p.view_index)];pop=np.array([p.cells/p.annotated_scope_cells,np.log1p(p.cells)]);pop=(pop-md['mu'][2:])/md['sd'][2:];oldq=quality[ids[0],2:];offset=-(pop-oldq)@md['coef'][:,4:].T
        for j,f in enumerate(QFIELDS):
            tr=md['training'][f];q=quality[ids,j];domain.append({'pool_index':pool,'covariate':f,'min':float(q.min()),'max':float(q.max()),'training_min':float(tr.min()),'training_max':float(tr.max()),'outside_training_cells':int(((q<tr.min())|(q>tr.max())).sum()),'cells':len(ids)})
        for pr,program in enumerate(['MCP1','MCP2']):
            for st,variant in enumerate(['full','fixed_stress_complement']):
                k=2*pr+st;records.append({'pool_index':pool,'scope':p.scope,'sampleID':p.sampleID,'view':p['view'],'program':program,'variant':variant,'score':float(scores[chosen,k].mean()),'all_captured_mean_sensitivity':float(scores[ids,k].mean()),'own_population_background_score':float(scores[chosen,k].mean()+offset[k]),'all_own_population_background_score':float(scores[ids,k].mean()+offset[k]),'unresidualized_weighted_RNA_score':float(unres[chosen,k].mean()),'RNA_quality_only_background_score':float(scores[chosen,k].mean()+oldq@md['coef'][k,4:]),'full_view_cells':len(ids),'selected_cells':len(chosen),'original_capture_cells':int(p.original_capture_cells),'all_material_cells':int(p.annotated_scope_cells),'population_fraction':float(p.cells/p.annotated_scope_cells),'population_logN':float(np.log1p(p.cells)),**{f:float(quality[chosen,j].mean()) for j,f in enumerate(QFIELDS)},**{'all_'+f:float(quality[ids,j].mean()) for j,f in enumerate(QFIELDS)}})
    measures=pd.DataFrame(records);old=pd.read_csv(native/'private_all3488_programme_measurements.tsv',sep='\t');anchors=measures[measures.scope.eq('legacy_capture_anchor')].merge(old,on=['sampleID','view','program','variant'],suffixes=('_new','_old'),validate='one_to_one');assert len(anchors)==824
    maxerr=0.
    for field in ['score','all_captured_mean_sensitivity',*QFIELDS,*['all_'+f for f in QFIELDS]]:
        err=float(np.max(np.abs(anchors[field+'_new']-anchors[field+'_old'])));assert err<1e-8;maxerr=max(maxerr,err)
    for x in [allraw,alldet,selraw,seldet]:x[:,col<0]=-1
    for x in [alllog,sellog]:x[:,col<0]=np.nan
    np.savez_compressed(out/'private_all_gene_pool_sufficient_statistics.npz',genes=np.array(genes),all_raw_sum=allraw,all_logCP10K_sum=alllog,all_detected_cells=alldet,selected_raw_sum=selraw,selected_logCP10K_sum=sellog,selected_detected_cells=seldet)
    uindex=np.searchsorted(rows,mr[selected]);np.savez_compressed(out/'private_selected_cell_quality_scores.npz',global_row=mr[selected].astype(np.int32),pool_index=pi[selected],fullRNA=RNA[uindex],detected=det[uindex],quality=quality[selected],scores=scores[selected]);measures.to_csv(out/'private_all_author_scope_programme_measurements.tsv',sep='\t',index=False);pd.DataFrame(domain).to_csv(out/'private_all_training_domain_positions.tsv',sep='\t',index=False);anchors.to_csv(out/'private_all824_legacy_capture_anchor_checks.tsv',sep='\t',index=False)
    (out/'private_all_unique_native_canonical_digests.json').write_text(json.dumps({'unique_cells':n,'memberships':len(mr),'canonical_nonzero_values':int(nnz),'digests':{k:v.hexdigest() for k,v in hashes.items()}},indent=2)+'\n')
    files={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()};size=sum(v['bytes'] for v in files.values());assert size<64*1024**2
    t={'status':'COMPLETE_ALL_AUTHOR_SCOPE_NATIVE_RNA_AND824_ANCHORS_INDEPENDENT_REPLAY_PENDING','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'scope_acceptance_sha256':sha(audit),'files':files,'unique_native_cells':n,'membership_rows':len(mr),'selected_memberships':int(selected.sum()),'pools':npools,'measurements':len(measures),'maximum_legacy_anchor_error':maxerr,'actual_bytes':size,'new_association_or_clinical_effect':False}
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps({'status':t['status'],'terminal_sha256':sha(out/'private_terminal.json'),'bytes':size}))


if __name__=='__main__':main()
