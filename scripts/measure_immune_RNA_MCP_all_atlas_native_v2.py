"""Complete prospective atlas RNA measurement, compact sufficient statistics."""
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.io import mmread
from scipy.sparse import csr_matrix
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']
STRESS={'FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'}
QFIELDS=['cellQ','geneQ','fraction','logN']
KEYS=['sampleID','view','program','variant']


def bound(path, expected):
    assert sha(path/'private_terminal.json')==expected
    t=json.loads((path/'private_terminal.json').read_text())
    for n,v in t['files'].items():assert sha(path/n)==(v['sha256'] if isinstance(v,dict) else v)
    return t


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);cap=128*1024**2;guard=resource(root,cap)
    opp=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';bound(opp,'2a37f7af2f73ac5ceddcb065626df411ca6f016271fce7f0830a1bd658f83cd1')
    contract=root/'immune_RNA_MCP_all_atlas_projection_contract_20261019_v1';bound(contract,'663c56f3aa8f01049b3d082b62ded4d4f51922cd5f714d56f5387a61ee950aca')
    acceptance=root/'immune_RNA_MCP_atlas_review_and_contract_audit_20261019_v1/private_terminal.json';assert sha(acceptance)=='9bc95d9880e360ed102c2200739366c50ef373d7ad4e48ee216d3ac1e79ad8aa'
    lib=pd.read_csv(opp/'private_all400_material_RNA_capture_scope.tsv',sep='\t',keep_default_na=False)
    axes=pd.read_csv(contract/'private_all2400_frozen_projection_axis_scope.tsv',sep='\t');assert axes.selected_for_projection.sum()==1744
    selected_axes=set(axes.loc[axes.selected_for_projection,['sampleID','view','program']].itertuples(index=False,name=None))
    roles=pd.read_csv(contract/'private_frozen_all_programme_LR_and_diagnostic_gene_roles.tsv',sep='\t');genes=roles.gene.tolist();assert len(genes)==904
    lib['material_index']=np.arange(400);lookup=pd.Index(lib.sampleID)
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');raw=cache['raw_gene_counts'];markers=cache['genes'].tolist()
    T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t');tc=lookup.get_indexer(T.dataset.astype(str)+'::'+T.library_id.astype(str))[cache['original_T_unit_index']]
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    gate=lambda names:np.count_nonzero(raw[:,[markers.index(g) for g in names]]>0,axis=1)
    core=gate(panels['pan_T'])>=2;reg=(gate(panels['regulatory_associated'])>=2)&(raw[:,markers.index('FOXP3')]>0);eff=gate(['GZMB','PRF1','CTSW','NKG7'])>=2
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1';units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t');mapping=np.load(oldscope/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    um=lookup.get_indexer(units.dataset.astype(str)+'::'+units.library_id.astype(str));um[~units.major.eq('Myeloid')]=-1
    rr=np.flatnonzero(mapping>=0);mc=um[mapping[rr]]
    candidates=[(cache['global_original_T_row_index'][core&reg&~eff],tc[core&reg&~eff]),(cache['global_original_T_row_index'][core&eff&~reg],tc[core&eff&~reg]),(rr[mc>=0],mc[mc>=0])]
    cellparts=[];poolrows=[]
    for m,material in lib.iterrows():
        for vi,view in enumerate(VIEWS):
            flags=[(material.sampleID,view,p) in selected_axes for p in ['MCP1','MCP2']]
            if not any(flags):continue
            rows=candidates[vi][0][candidates[vi][1]==m];assert len(rows)==material[view]>=20
            pool=len(poolrows);chosen=set(sorted(rows,key=lambda r:hashlib.sha256(('2026101802|'+view+'|'+str(int(r))).encode()).hexdigest())[:100])
            cellparts.append(pd.DataFrame(dict(global_row=rows,pool_index=pool,selected=[int(r) in chosen for r in rows])))
            poolrows.append(dict(pool_index=pool,material_index=m,sampleID=material.sampleID,dataset=material.dataset,view=view,view_index=vi,full_view_cells=len(rows),selected_cells=len(chosen),all_material_cells=material.annotated_scope_cells,MCP1=flags[0],MCP2=flags[1]))
    pools=pd.DataFrame(poolrows);cells=pd.concat(cellparts,ignore_index=True).sort_values('global_row').reset_index(drop=True);assert cells.global_row.is_unique and pools.sampleID.nunique()==355
    n=len(cells);nu=len(pools);ng=len(genes);ns=int(cells.selected.sum())
    # Worst-case dense arrays plus conservative text/ZIP overhead. No selected
    # raw-count matrix is written; both implementations replay every native row.
    dense_bound=n*11+ns*112+nu*ng*48+nu*16384+904*4096+8*1024**2
    assert dense_bound<cap
    out=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,parent256MiB_contract_sha256=sha(contract/'private_terminal.json'),acceptance_sha256=sha(acceptance),output_cap_bytes=cap,conservative_uncompressed_bound_bytes=int(dense_bound),
        resource_revision='Prospective128MiB compact output replaces unexecuted256MiB reservation only. All scientific axes/genes/materials/selection/models unchanged; no deletion/reserve reduction. Same full raw replay required.',
        scope='All1744 complete selected axes,355materials; every captured cell, original hash100 primary, allcapture means and own qualities, all fixed394LR components/43diagnostics/34capture in same native read. No effects/outcomes.',
        storage='Full captured row/pool/selection membership; fixed selected cell qualities/scores; complete primary/allcapture per-pool per-gene raw sums/logCP10K sums/detection counts with source missingness; canonical full-row native counts/depth digests. All individual full-cell scores replayed from original counts and frozen model via independent audit.',
        no_new_expression_download_or_gene_fit_or_clinical_effect=True)
    (out/'plan_before_any_new_native_expression_compact128.json').write_text(json.dumps(plan,indent=2)+'\n')
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');read=lambda k:np.asarray(read_elem(z[k])).astype(str)
    source=read('obs/dataset');library=read('obs/library_id');allmaterial=lookup.get_indexer(np.char.add(np.char.add(source,'::'),library));assert (allmaterial>=0).all()
    assert np.array_equal(np.bincount(allmaterial,minlength=400),lib.annotated_scope_cells.to_numpy(int))
    order=read('uns/source_order').tolist();measured=np.asarray(read_elem(z['varm/measured_in_source']),bool);symbol=read('var/gene_symbol');canonical=read('var/canonical_gene_identity');identity=read('var/feature_identity_status');alias=read('var/source_gene_aliases')
    datasets=sorted(lib.dataset.unique());genemap=np.full((len(datasets),ng),-1,int);maprows=[]
    for si,dataset in enumerate(datasets):
        names=defaultdict(set)
        for col in np.flatnonzero(measured[:,order.index(dataset)]):
            if identity[col]=='ambiguous_source_feature' or canonical[col].startswith('AMBIGUOUS::'):continue
            labels={symbol[col],canonical[col]}
            if identity[col]=='stable_Ensembl_ID':labels.update(alias[col].split('|'))
            for gene in labels&set(genes):names[gene].add(int(col))
        for gi,gene in enumerate(genes):
            matches=names[gene];col=next(iter(matches)) if len(matches)==1 else -1;genemap[si,gi]=col
            maprows.append(dict(dataset=dataset,gene=gene,native_column=col,status='UNIQUE_SOURCE_MEASURED' if col>=0 else 'MULTIPLE_NATIVE_MATCHES' if matches else 'UNMEASURABLE',roles=roles.iloc[gi].roles))
    originalmaps=pd.read_csv(opp/'private_all_active_gene_source_native_mapping.tsv.gz',sep='\t')
    for r in originalmaps.itertuples(index=False):assert genemap[datasets.index(r.dataset),genes.index(r.gene)]==r.native_column
    marker_contract={r['dataset']:r for r in json.loads((old/'private_source_gene_measurement_contract.json').read_text())}
    for si,dataset in enumerate(datasets):
        c=marker_contract[dataset];assert c['requested_genes']==markers
        assert np.array_equal(genemap[si,[genes.index(g) for g in markers]],c['physical_columns'])
    physical=np.unique(genemap[genemap>=0]);upos=np.full_like(genemap,-1);upos[genemap>=0]=np.searchsorted(physical,genemap[genemap>=0])
    wpath=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wpath)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7'
    weights=pd.read_csv(wpath,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    traininputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs';models={};records=[]
    for vi,view in enumerate(VIEWS):
        md=pd.read_csv(traininputs/(view+'_discovery_metadata.tsv'),sep='\t');D=np.c_[np.ones(len(md)),md[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]
        names=Path(str(traininputs/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        with gzip.open(str(traininputs/(view+'_discovery'))+'_tpm.mtx.gz','rb') as h:tr=mmread(h).tocsr()
        frozen=np.load(traininputs/(view+'_frozen_projection.npz'));B=np.zeros((ng,4));centers=[];coefs=[]
        for pi,program in enumerate(['MCP1','MCP2']):
            w=weights[weights.view.eq(view)&weights.program.eq(program)]
            for st,variant in enumerate(['full','fixed_stress_complement']):
                wv=w if st==0 else w[~w.gene.str.upper().isin(STRESS)];beta=wv.weight.to_numpy()/wv.discovery_sd.to_numpy();mu=float(beta@wv.discovery_mean.to_numpy());k=2*pi+st
                B[[genes.index(g) for g in wv.gene],k]=beta;centers.append(mu);tidx=pd.Index(names).get_indexer(wv.gene);assert (tidx>=0).all()
                coeff=np.linalg.lstsq(D,np.asarray(tr[tidx].T.dot(beta)).ravel()-mu,rcond=None)[0];coefs.append(coeff)
                records.append(dict(view=view,program=program,variant=variant,active_genes=len(wv),center=mu,**{key:float(value) for key,value in zip(['intercept','cellQ','geneQ','sourceCC','fraction','logN'],coeff)}))
        models[vi]=dict(B=B,center=np.array(centers),coef=np.array(coefs),mu=frozen['quality_mean'],sd=frozen['quality_sd'],training=md)
    counts=z['layers/counts'];nc,nf=map(int,counts.attrs['shape']);assert len(mapping)==nc
    rows=cells.global_row.to_numpy(int);pi=cells.pool_index.to_numpy(int);selected=cells.selected.to_numpy(bool)
    allraw=np.zeros((nu,ng),np.int64);alllog=np.zeros((nu,ng));alldet=np.zeros((nu,ng),np.int64);selraw=allraw.copy();sellog=alllog.copy();seldet=alldet.copy()
    RNA=np.empty(n,np.int64);det=np.empty(n,np.int32);qs=np.empty((n,4));values=np.full((n,4),np.nan)
    hashes={k:hashlib.sha256() for k in ['native_data_int64','native_indices_int64','native_row_endpoints_int64','row_ids_int64','fullRNA_int64','detected_int64']};nnz=0
    for start in range(0,nc,8192):
        stop=min(nc,start+8192);lo,hi=np.searchsorted(rows,[start,stop])
        if hi>lo:
            ptr=np.asarray(counts['indptr'][start:stop+1],np.int64);x=csr_matrix((np.asarray(counts['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(counts['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,nf))[rows[lo:hi]-start].copy();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.rint(x.data)).all() and x.data.max()<2**53
            x=x.astype(np.int64);RNA[lo:hi]=np.asarray(x.sum(axis=1)).ravel();det[lo:hi]=np.diff(x.indptr);assert (RNA[lo:hi]>0).all()
            for key,a in [('native_data_int64',x.data),('native_indices_int64',x.indices),('native_row_endpoints_int64',x.indptr[1:].astype(np.int64)+nnz),('row_ids_int64',rows[lo:hi]),('fullRNA_int64',RNA[lo:hi]),('detected_int64',det[lo:hi])]:hashes[key].update(np.asarray(a,dtype='<i8').tobytes())
            nnz+=x.nnz;rawunion=x[:,physical].toarray();logged=np.log1p(rawunion.astype(float)*10000/RNA[lo:hi,None]);pids=pi[lo:hi]
            for pool in np.unique(pids):
                local=np.flatnonzero(pids==pool);ids=lo+local;p=pools.iloc[pool];vi=int(p.view_index);si=datasets.index(p.dataset);gi=np.flatnonzero(upos[si]>=0);columns=upos[si,gi];d=rawunion[local][:,columns];l=logged[local][:,columns];take=selected[ids]
                allraw[pool,gi]+=d.sum(axis=0);alllog[pool,gi]+=l.sum(axis=0);alldet[pool,gi]+=np.count_nonzero(d,axis=0)
                selraw[pool,gi]+=d[take].sum(axis=0);sellog[pool,gi]+=l[take].sum(axis=0);seldet[pool,gi]+=np.count_nonzero(d[take],axis=0)
                model=models[vi];q=np.c_[np.log1p(RNA[ids]),np.log1p(det[ids]),np.repeat(p.full_view_cells/p.all_material_cells,len(ids)),np.repeat(np.log1p(p.full_view_cells),len(ids))];q=(q-model['mu'])/model['sd'];qs[ids]=q
                sourceCC=float(p.dataset=='GSE294109');design=np.c_[np.ones(len(ids)),q[:,:2],np.repeat(sourceCC,len(ids)),q[:,2:]]
                B=model['B'][gi];score=l@B-model['center']-design@model['coef'].T
                for pr in range(2):
                    if getattr(p,'MCP'+str(pr+1)):values[ids,2*pr:2*pr+2]=score[:,2*pr:2*pr+2]
                if vi<2:
                    tix=np.searchsorted(cache['global_original_T_row_index'],rows[ids]);assert np.array_equal(cache['global_original_T_row_index'][tix],rows[ids]);assert np.array_equal(rawunion[local][:,upos[si,[genes.index(g) for g in markers]]],raw[tix]);assert np.array_equal(RNA[ids],cache['full_native_total_counts'][tix]) and np.array_equal(det[ids],cache['full_native_detected_genes'][tix])
        if start//8192%16==0 or stop==nc:(out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),phase='FULL_NATIVE_ROWS',native_rows_scanned=stop,native_rows=nc,captured_cells_read=int(hi),captured_cells=n))+'\n')
    assert np.isfinite(qs).all()
    for pool,p in pools.iterrows():
        missing=genemap[datasets.index(p.dataset)]<0
        for a in [allraw,alldet,selraw,seldet]:a[pool,missing]=-1
        for a in [alllog,sellog]:a[pool,missing]=np.nan
    case=[];domain=[]
    for pool,p in pools.iterrows():
        idx=np.flatnonzero(pi==pool);sub=idx[selected[idx]];assert len(idx)==p.full_view_cells and len(sub)==p.selected_cells
        for j,f in enumerate(QFIELDS):
            lower=models[int(p.view_index)]['training'][f].min();upper=models[int(p.view_index)]['training'][f].max();v=qs[idx,j]
            domain.append(dict(pool_index=pool,view=p['view'],covariate=f,training_min=lower,training_max=upper,new_min=float(v.min()),new_max=float(v.max()),outside_training_cells=int(((v<lower)|(v>upper)).sum()),captured_cells=len(v)))
        for pr,program in enumerate(['MCP1','MCP2']):
            if not getattr(p,program):continue
            for st,variant in enumerate(['full','fixed_stress_complement']):
                k=2*pr+st;assert np.isfinite(values[idx,k]).all()
                case.append(dict(sampleID=p.sampleID,dataset=p.dataset,pool_index=pool,view=p['view'],program=program,variant=variant,score=float(values[sub,k].mean()),all_captured_mean_sensitivity=float(values[idx,k].mean()),full_view_cells=len(idx),selected_cells=len(sub),all_material_cells=p.all_material_cells,**{f:float(qs[sub,j].mean()) for j,f in enumerate(QFIELDS)},**{'all_'+f:float(qs[idx,j].mean()) for j,f in enumerate(QFIELDS)}))
    case=pd.DataFrame(case);assert len(case)==3488
    official=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t');official=official[official.fit_id.eq('K2_seed20261018_source_adjusted')].copy();assert len(official)==396
    partition=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_metadata_only_case_source_partition.tsv',sep='\t').set_index('fit_sample').sampleID;official['sampleID']=official['sample'].map(partition);assert official.sampleID.notna().all()
    anchors=case[case.variant.eq('full')].merge(official[['sampleID','view','program','score']],on=['sampleID','view','program'],suffixes=('_new','_accepted'),validate='one_to_one');assert len(anchors)==396
    maximum=float(np.max(np.abs(anchors.score_new-anchors.score_accepted)));assert maximum<1e-8
    paired=pd.read_csv(root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1/private_all624_material_programme_measurements.tsv',sep='\t')
    pairedanchors=case.merge(paired[KEYS+['score','all_captured_mean_sensitivity']],on=KEYS,suffixes=('_new','_accepted'),validate='one_to_one');assert len(pairedanchors)==624
    for f in ['score','all_captured_mean_sensitivity']:assert np.max(np.abs(pairedanchors[f+'_new']-pairedanchors[f+'_accepted']))<1e-8
    np.savez_compressed(out/'private_all_captured_membership.npz',global_row=rows,pool_index=pi.astype(np.int16),selected=selected)
    np.savez_compressed(out/'private_selected_cells_quality_scores.npz',global_row=rows[selected],pool_index=pi[selected].astype(np.int16),fullRNA=RNA[selected],detected=det[selected],quality=qs[selected],scores=values[selected])
    np.savez_compressed(out/'private_all_gene_measurement_sufficient_statistics.npz',genes=np.asarray(genes),all_raw_sum=allraw,all_logCP10K_sum=alllog,all_detected_cells=alldet,selected_raw_sum=selraw,selected_logCP10K_sum=sellog,selected_detected_cells=seldet)
    pools.to_csv(out/'private_all_capture_pools.tsv',sep='\t',index=False);lib.to_csv(out/'private_all400_materials.tsv',sep='\t',index=False);axes.to_csv(out/'private_all2400_axis_scope.tsv',sep='\t',index=False);pd.DataFrame(maprows).to_csv(out/'private_all_source_auxiliary_native_gene_mapping.tsv.gz',sep='\t',index=False)
    pd.DataFrame(records).to_csv(out/'private_all12_frozen_models.tsv',sep='\t',index=False);case.to_csv(out/'private_all3488_programme_measurements.tsv',sep='\t',index=False);pd.DataFrame(domain).to_csv(out/'private_all_capture_training_domain_descriptions.tsv',sep='\t',index=False);anchors.to_csv(out/'private_all396_official_anchor_checks.tsv',sep='\t',index=False);pairedanchors.to_csv(out/'private_all624_paired_anchor_checks.tsv',sep='\t',index=False)
    (out/'private_canonical_full_native_digests.json').write_text(json.dumps(dict(captured_cells=n,canonical_nonzero_values=int(nnz),digests={k:v.hexdigest() for k,v in hashes.items()}),indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<cap
    terminal=dict(status='COMPLETE_ALL_SELECTED_ATLAS_NATIVE_RNA_FULL_PROGRAMMES_AND_AUXILIARY_NO_EFFECTS_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,selected_complete_axes=1744,materials=355,captured_cells=n,primary_selected_cells=ns,capture_pools=nu,auxiliary_genes=ng,canonical_native_nonzero_values=int(nnz),maximum_original_anchor_error=maximum,no_new_association_or_clinical_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),captured_cells=n,selected_cells=ns,capture_pools=nu,output_bytes=sum(v['bytes'] for v in files.values()))))


if __name__=='__main__':main()
