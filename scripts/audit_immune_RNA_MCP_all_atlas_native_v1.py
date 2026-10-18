"""Entire atlas assay replay: independent membership, batches, QR and sums."""
import gzip
import hashlib
import json
import os
import time
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


def close(a,b,atol=1e-8):
    assert np.allclose(a,b,rtol=0,atol=atol,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    producer=root/'immune_RNA_MCP_all_atlas_native_measurement_20261019_v2'
    launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2/measure_immune_RNA_MCP_all_atlas_native_v2/launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['source_sha256']=='8cf59df7d20485ed960029282f7d1b66038afa45caaf198ece1ffe6ff41e37e2'
    out=root/'immune_RNA_MCP_all_atlas_native_measurement_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_native_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,launch_sha256=sha(launchpath),scope='Global-row assignment rather than per-material candidate grouping; global hash sort; required-name-only native map; captured1024-row batches; all canonical native/depth digests, exact raw/detection sufficient statistics and mean-log quantities, selected individual qualities/scores, all material full/stress/own-quality means, domains and every396official/624paired anchor. QR rather than least-squares coefficients.',no_new_association_or_clinical_effects=True),indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not (producer/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'PRODUCER_MISSING_WITHOUT_TERMINAL';st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z';assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(phase='WAITING_SAME_LIVE_NATIVE_PRODUCER',pid=launch['pid'],startticks=st[19],observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n');time.sleep(15)
    term=json.loads((producer/'private_terminal.json').read_text());assert term['status']=='COMPLETE_ALL_SELECTED_ATLAS_NATIVE_RNA_FULL_PROGRAMMES_AND_AUXILIARY_NO_EFFECTS_AUDIT_PENDING'
    for n,v in term['files'].items():assert sha(producer/n)==v['sha256']
    opp=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';lib=pd.read_csv(opp/'private_all400_material_RNA_capture_scope.tsv',sep='\t',keep_default_na=False);lib['material_index']=np.arange(400)
    pd.testing.assert_frame_equal(pd.read_csv(producer/'private_all400_materials.tsv',sep='\t',keep_default_na=False),lib)
    contract=root/'immune_RNA_MCP_all_atlas_projection_contract_20261019_v1';axes=pd.read_csv(contract/'private_all2400_frozen_projection_axis_scope.tsv',sep='\t');pd.testing.assert_frame_equal(pd.read_csv(producer/'private_all2400_axis_scope.tsv',sep='\t'),axes)
    accept=root/'immune_RNA_MCP_atlas_review_and_contract_audit_20261019_v1/private_terminal.json';assert sha(accept)=='9bc95d9880e360ed102c2200739366c50ef373d7ad4e48ee216d3ac1e79ad8aa'
    pools=pd.read_csv(producer/'private_all_capture_pools.tsv',sep='\t');assert np.array_equal(pools.pool_index,np.arange(len(pools)))
    keys={(r.sampleID,r.view,r.program) for r in axes[axes.selected_for_projection].itertuples(index=False)}
    wanted=np.full((400,3),-1,int);records=[]
    for m,material in lib.iterrows():
        for vi,view in enumerate(VIEWS):
            pr=[(material.sampleID,view,p) in keys for p in ['MCP1','MCP2']]
            if not any(pr):continue
            pool=len(records);wanted[m,vi]=pool;records.append((pool,m,material.sampleID,material.dataset,view,vi,int(material[view]),min(100,int(material[view])),int(material.annotated_scope_cells),*pr))
    pd.testing.assert_frame_equal(pools,pd.DataFrame(records,columns=pools.columns),check_dtype=False)
    lookup={s:i for i,s in enumerate(lib.sampleID)}
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');raw=cache['raw_gene_counts'];markers=cache['genes'].tolist();T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t')
    tcase=np.array([lookup[str(r.dataset)+'::'+str(r.library_id)] for r in T.itertuples(index=False)])[cache['original_T_unit_index']]
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    positive=lambda names:np.count_nonzero(raw[:,[markers.index(g) for g in names]]>=1,axis=1)
    core=positive(panels['pan_T'])>=2;reg=(positive(panels['regulatory_associated'])>=2)&(positive(['FOXP3'])==1);eff=positive(['GZMB','PRF1','CTSW','NKG7'])>=2
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1';mapping=np.load(oldscope/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r');units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t')
    assignment=np.full(len(mapping),-1,np.int32)
    for vi,take in enumerate([core&reg&~eff,core&eff&~reg]):
        rr=cache['global_original_T_row_index'][take];assignment[rr]=wanted[tcase[take],vi]
    unitpool=np.array([wanted[lookup[str(r.dataset)+'::'+str(r.library_id)],2] if r.major=='Myeloid' else -1 for r in units.itertuples(index=False)])
    rr=np.flatnonzero(mapping>=0);pp=unitpool[mapping[rr]];take=pp>=0;assert (assignment[rr[take]]<0).all();assignment[rr[take]]=pp[take]
    rows=np.flatnonzero(assignment>=0);pi=assignment[rows];assert len(rows)==885670
    pool_views=pools['view'].to_numpy(str)
    digest=np.array([hashlib.sha256(('2026101802|'+pool_views[p]+'|'+str(r)).encode()).hexdigest() for r,p in zip(rows,pi)],dtype='S64')
    order=np.lexsort((digest,pi));ordered_p=pi[order];first=np.r_[0,np.flatnonzero(ordered_p[1:]!=ordered_p[:-1])+1];offset=np.arange(len(rows))-np.repeat(first,np.diff(np.r_[first,len(rows)]));chosen=np.zeros(len(rows),bool);chosen[order[offset<100]]=True
    saved=np.load(producer/'private_all_captured_membership.npz');assert np.array_equal(saved['global_row'],rows) and np.array_equal(saved['pool_index'],pi) and np.array_equal(saved['selected'],chosen)
    assert np.array_equal(np.bincount(pi),pools.full_view_cells) and np.array_equal(np.bincount(pi[chosen]),pools.selected_cells)
    roles=pd.read_csv(contract/'private_frozen_all_programme_LR_and_diagnostic_gene_roles.tsv',sep='\t');genes=roles.gene.tolist();ng=len(genes);nu=len(pools)
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');read=lambda key:np.asarray(read_elem(z[key])).astype(str)
    source=read('obs/dataset');libraries=read('obs/library_id');materialkeys=np.char.add(np.char.add(source,'::'),libraries);obs_material=pd.Index(lib.sampleID).get_indexer(materialkeys);assert (obs_material>=0).all() and np.array_equal(np.bincount(obs_material,minlength=400),lib.annotated_scope_cells)
    datasets=sorted(lib.dataset.unique());sorder=read('uns/source_order').tolist();measured=np.asarray(read_elem(z['varm/measured_in_source']),bool);symbols=read('var/gene_symbol');canonical=read('var/canonical_gene_identity');identity=read('var/feature_identity_status');aliases=read('var/source_gene_aliases')
    mp=pd.read_csv(producer/'private_all_source_auxiliary_native_gene_mapping.tsv.gz',sep='\t').set_index(['dataset','gene']);gm=np.full((len(datasets),ng),-1,int);requested=set(genes)
    for si,dataset in enumerate(datasets):
        matches={g:set() for g in genes}
        for col in np.flatnonzero(measured[:,sorder.index(dataset)]):
            if canonical[col].startswith('AMBIGUOUS::') or identity[col]=='ambiguous_source_feature':continue
            names={canonical[col],symbols[col]}
            if identity[col]=='stable_Ensembl_ID':names|=set(aliases[col].split('|'))
            for g in names&requested:matches[g].add(int(col))
        for gi,g in enumerate(genes):
            found=matches[g];col=next(iter(found)) if len(found)==1 else -1;gm[si,gi]=col;r=mp.loc[(dataset,g)]
            assert r.native_column==col and r.status==('UNIQUE_SOURCE_MEASURED' if col>=0 else 'MULTIPLE_NATIVE_MATCHES' if found else 'UNMEASURABLE') and r.roles==roles.iloc[gi].roles
    assert len(mp)==len(datasets)*ng
    wpath=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wpath)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7';weights=pd.read_csv(wpath,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs';models={};params=pd.read_csv(producer/'private_all12_frozen_models.tsv',sep='\t').set_index(['view','program','variant'])
    for vi,view in enumerate(VIEWS):
        md=pd.read_csv(inputs/(view+'_discovery_metadata.tsv'),sep='\t');D=np.column_stack([np.ones(len(md)),md[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]);Q,R=np.linalg.qr(D,mode='reduced')
        names=Path(str(inputs/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        with gzip.open(str(inputs/(view+'_discovery'))+'_tpm.mtx.gz','rb') as h:tr=mmread(h).tocsr()
        frozen=np.load(inputs/(view+'_frozen_projection.npz'));models[vi]=dict(mu=frozen['quality_mean'],sd=frozen['quality_sd'],training=md,pieces={})
        for pr,program in enumerate(['MCP1','MCP2']):
            w=weights[weights.view.eq(view)&weights.program.eq(program)];keep=~w.gene.str.upper().isin(STRESS).to_numpy()
            for st,variant in enumerate(['full','fixed_stress_complement']):
                use=np.ones(len(w),bool) if st==0 else keep;wv=w[use];b=wv.weight.to_numpy()/wv.discovery_sd.to_numpy();center=float(np.sum(b*wv.discovery_mean.to_numpy()));tid=pd.Index(names).get_indexer(wv.gene);assert (tid>=0).all();y=np.asarray(tr[tid].T.multiply(b).sum(axis=1)).ravel()-center;coef=np.linalg.solve(R,Q.T@y)
                p=params.loc[(view,program,variant)];close(p[['intercept','cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float),coef);close(p.center,center);assert p.active_genes==len(wv)
                models[vi]['pieces'][2*pr+st]=dict(gene_index=np.array([genes.index(g) for g in wv.gene]),b=b,center=center,coef=coef)
    native=z['layers/counts'];nc,nf=map(int,native.attrs['shape']);assert len(mapping)==nc
    allraw=np.zeros((nu,ng),np.int64);alllog=np.zeros((nu,ng));alldet=np.zeros((nu,ng),np.int64);selraw=allraw.copy();sellog=alllog.copy();seldet=alldet.copy();n=len(rows)
    RNA=np.empty(n,np.int64);det=np.empty(n,np.int64);quality=np.empty((n,4));scores=np.full((n,4),np.nan)
    hashnames=['native_data_int64','native_indices_int64','native_row_endpoints_int64','row_ids_int64','fullRNA_int64','detected_int64'];hashes={key:hashlib.sha256() for key in hashnames};nnz=0
    physical=np.unique(gm[gm>=0]);columns=np.full_like(gm,-1);columns[gm>=0]=np.searchsorted(physical,gm[gm>=0])
    for lo in range(0,n,1024):
        hi=min(n,lo+1024);batch=rows[lo:hi];start=int(batch[0]);stop=int(batch[-1])+1;ptr=np.asarray(native['indptr'][start:stop+1],np.int64)
        x=csr_matrix((np.asarray(native['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(native['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(stop-start,nf))[batch-start].copy();x.sum_duplicates();x.eliminate_zeros();x.sort_indices();assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.floor(x.data)).all() and x.data.max()<2**53
        x=x.astype(np.int64);RNA[lo:hi]=np.asarray(x.sum(axis=1)).ravel();det[lo:hi]=np.diff(x.indptr)
        for key,a in [('native_data_int64',x.data),('native_indices_int64',x.indices),('native_row_endpoints_int64',x.indptr[1:].astype(np.int64)+nnz),('row_ids_int64',batch),('fullRNA_int64',RNA[lo:hi]),('detected_int64',det[lo:hi])]:hashes[key].update(np.asarray(a,dtype='<i8').tobytes())
        nnz+=x.nnz;dense=x[:,physical].toarray();log=np.log1p(dense.astype(float)*10000/RNA[lo:hi,None])
        for pool in set(pi[lo:hi].tolist()):
            local=np.flatnonzero(pi[lo:hi]==pool);ids=lo+local;p=pools.iloc[pool];si=datasets.index(p.dataset);vi=int(p.view_index);valid=np.flatnonzero(gm[si]>=0);u=columns[si,valid];d=dense[local][:,u];l=log[local][:,u];take=chosen[ids]
            allraw[pool,valid]+=np.add.reduce(d,axis=0);alllog[pool,valid]+=np.add.reduce(l,axis=0);alldet[pool,valid]+=(d>0).sum(axis=0)
            selraw[pool,valid]+=np.add.reduce(d[take],axis=0);sellog[pool,valid]+=np.add.reduce(l[take],axis=0);seldet[pool,valid]+=(d[take]>0).sum(axis=0)
            md=models[vi];unscaled=np.column_stack([np.log1p(RNA[ids]),np.log1p(det[ids]),np.full(len(ids),p.full_view_cells/p.all_material_cells),np.full(len(ids),np.log1p(p.full_view_cells))]);q=(unscaled-md['mu'])/md['sd'];quality[ids]=q
            design=np.column_stack([np.ones(len(ids)),q[:,:2],np.full(len(ids),float(p.dataset=='GSE294109')),q[:,2:]])
            for k,model in md['pieces'].items():
                if not getattr(p,'MCP'+str(k//2+1)):continue
                ui=columns[si,model['gene_index']];assert (ui>=0).all();scores[ids,k]=np.sum(log[local][:,ui]*model['b'],axis=1)-model['center']-design@model['coef']
            if vi<2:
                ci=np.searchsorted(cache['global_original_T_row_index'],rows[ids]);assert np.array_equal(cache['global_original_T_row_index'][ci],rows[ids]);assert np.array_equal(dense[local][:,columns[si,[genes.index(g) for g in markers]]],raw[ci]);assert np.array_equal(RNA[ids],cache['full_native_total_counts'][ci]) and np.array_equal(det[ids],cache['full_native_detected_genes'][ci])
        if lo//1024%16==0 or hi==n:(out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_FULL_NATIVE_REPLAY',captured_cells_replayed=hi,captured_cells=n,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    original_digest=json.loads((producer/'private_canonical_full_native_digests.json').read_text());assert original_digest['captured_cells']==n and original_digest['canonical_nonzero_values']==nnz and original_digest['digests']=={k:v.hexdigest() for k,v in hashes.items()}
    saved=np.load(producer/'private_selected_cells_quality_scores.npz');assert np.array_equal(saved['global_row'],rows[chosen]) and np.array_equal(saved['pool_index'],pi[chosen]) and np.array_equal(saved['fullRNA'],RNA[chosen]) and np.array_equal(saved['detected'],det[chosen]);close(saved['quality'],quality[chosen]);close(saved['scores'],scores[chosen]);score_error=float(np.nanmax(np.abs(saved['scores']-scores[chosen])))
    stats=np.load(producer/'private_all_gene_measurement_sufficient_statistics.npz');assert stats['genes'].tolist()==genes
    log_error=0.
    for pool,p in pools.iterrows():
        invalid=gm[datasets.index(p.dataset)]<0
        for a in [allraw,alldet,selraw,seldet]:a[pool,invalid]=-1
        for a in [alllog,sellog]:a[pool,invalid]=np.nan
    for name,a in [('all_raw_sum',allraw),('all_detected_cells',alldet),('selected_raw_sum',selraw),('selected_detected_cells',seldet)]:assert np.array_equal(stats[name],a)
    for name,a,denom in [('all_logCP10K_sum',alllog,pools.full_view_cells.to_numpy()),('selected_logCP10K_sum',sellog,pools.selected_cells.to_numpy())]:
        close(stats[name]/denom[:,None],a/denom[:,None],1e-10);log_error=max(log_error,float(np.nanmax(np.abs(stats[name]/denom[:,None]-a/denom[:,None]))))
    case=pd.read_csv(producer/'private_all3488_programme_measurements.tsv',sep='\t');assert len(case)==3488 and not case.duplicated(KEYS).any();caseidx=case.set_index(KEYS);domain=pd.read_csv(producer/'private_all_capture_training_domain_descriptions.tsv',sep='\t').set_index(['pool_index','covariate']);assert len(domain)==nu*4
    for pool,p in pools.iterrows():
        ids=np.flatnonzero(pi==pool);sel=ids[chosen[ids]]
        for j,f in enumerate(QFIELDS):
            t=models[int(p.view_index)]['training'][f];v=quality[ids,j];d=domain.loc[(pool,f)];close([d.training_min,d.training_max,d.new_min,d.new_max],[t.min(),t.max(),v.min(),v.max()]);assert d.outside_training_cells==int(((v<t.min())|(v>t.max())).sum()) and d.captured_cells==len(v) and d['view']==p['view']
        for pr,program in enumerate(['MCP1','MCP2']):
            if not getattr(p,program):continue
            for st,variant in enumerate(['full','fixed_stress_complement']):
                r=caseidx.loc[(p.sampleID,p['view'],program,variant)];k=2*pr+st;close([r.score,r.all_captured_mean_sensitivity],[np.mean(scores[sel,k]),np.mean(scores[ids,k])]);close(r[QFIELDS].to_numpy(float),np.mean(quality[sel],axis=0));close(r[['all_'+f for f in QFIELDS]].to_numpy(float),np.mean(quality[ids],axis=0));assert r.pool_index==pool and r.dataset==p.dataset and r.full_view_cells==len(ids) and r.selected_cells==len(sel) and r.all_material_cells==p.all_material_cells
    partition=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_metadata_only_case_source_partition.tsv',sep='\t').set_index('fit_sample').sampleID
    official=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t');official=official[official.fit_id.eq('K2_seed20261018_source_adjusted')].copy();assert len(official)==396
    for r in official.itertuples(index=False):close(caseidx.loc[(partition.loc[r.sample],r.view,r.program,'full'),'score'],r.score)
    paired=pd.read_csv(root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1/private_all624_material_programme_measurements.tsv',sep='\t');assert len(paired)==624
    for r in paired.itertuples(index=False):close(caseidx.loc[(r.sampleID,r.view,r.program,r.variant),['score','all_captured_mean_sensitivity']].to_numpy(float),[r.score,r.all_captured_mean_sensitivity])
    savedofficial=pd.read_csv(producer/'private_all396_official_anchor_checks.tsv',sep='\t');savedpaired=pd.read_csv(producer/'private_all624_paired_anchor_checks.tsv',sep='\t');assert len(savedofficial)==396 and len(savedpaired)==624
    for f in ['score']:close(savedofficial[f+'_new'],savedofficial[f+'_accepted'])
    for f in ['score','all_captured_mean_sensitivity']:close(savedpaired[f+'_new'],savedpaired[f+'_accepted'])
    terminal=dict(status='ACCEPTED_COMPLETE_ALL_SELECTED_ATLAS_NATIVE_COUNTS_MODELS_SELECTED_SCORES_ALL_SUMMARIES_AND_ANCHORS_NO_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(producer/'private_terminal.json'),captured_cells=n,selected_cells=int(chosen.sum()),capture_pools=nu,complete_axes=1744,all_gene_pool_positions=nu*ng,canonical_native_nonzero_values=int(nnz),all_official_anchors=396,all_paired_anchors=624,all_programme_measurements=3488,maximum_selected_score_error=score_error,maximum_mean_log_gene_error=log_error,no_new_association_or_clinical_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),captured_cells=n,maximum_selected_score_error=score_error)))


if __name__=='__main__':main()
