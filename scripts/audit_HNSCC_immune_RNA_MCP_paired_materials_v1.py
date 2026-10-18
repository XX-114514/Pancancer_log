"""Independent row batches, dense sums/QR and all paired-material assay values."""
import gzip
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy.io import mmread
from scipy.sparse import csr_matrix, load_npz

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

VIEWS=['R_FOXP3_capture','E_effector_capture','Myeloid']
STRESS={'FOS','FOSB','JUN','JUNB','JUND','ATF3','EGR1','DUSP1','HSPA1A','HSPA1B','DNAJB1','IER2','IER3','NR4A1','BTG2'}
FIELDS=['cellQ','geneQ','fraction','logN']


def close(a,b):
    assert np.allclose(a,b,rtol=0,atol=1e-8,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_20261019_v1'
    launchpath=root/'HNSCC_receiver_followup_job_control_20261016_v2/measure_HNSCC_immune_RNA_MCP_paired_materials_v1/launch_receipt.json'
    launch=json.loads(launchpath.read_text());assert launch['source_sha256']=='21c7010b54590a7b86c3a46a6d55cd2ef7a066b6ace6887012da3c04cca4462d'
    out=root/'HNSCC_immune_RNA_MCP_paired_native_measurement_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_native_paired_material_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,producer_launch_sha256=sha(launchpath),
        scope='Reconstruct every captured global-row membership/primary hash selection; different1024captured-row reads, all fullRNA/detected and selected physical counts, dense expression sums and independent QR coefficients, all624 material measurements/12domains/6scales and original official anchors.',
        no_paired_difference_correlation_or_clinical_effects=True),indent=2)+'\n')
    deadline=time.monotonic()+21600
    while not (source/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'PRODUCER_MISSING_WITHOUT_TERMINAL'
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline
        (out/'private_progress.json').write_text(json.dumps(dict(phase='WAITING_SAME_LIVE_PAIRED_NATIVE_PRODUCER',observed_UTC=datetime.now(timezone.utc).isoformat(),pid=launch['pid'],startticks=st[19]))+'\n');time.sleep(15)
    term=json.loads((source/'private_terminal.json').read_text());assert term['status']=='COMPLETE_ALL52_PAIRED_NATIVE_FULL_FROZEN_RNA_MEASUREMENTS_NO_EFFECTS_AUDIT_PENDING'
    for n,v in term['files'].items():assert sha(source/n)==v['sha256']
    paired=root/'HNSCC_immune_RNA_MCP_paired_materials_20261019_v1';assert sha(paired/'private_independent_metadata_acceptance_20261019_v1.json')=='8dc0f5b74ea8cf41ff10a62eb5c804064a0dc558ec90e787f4e9cbac5b4461c9'
    pp=pd.read_csv(paired/'private_all_exact_pre_post_patient_material_pairs.tsv',sep='\t');pp=pp[~pp.post_used_in_original_discovery]
    material_rows=[(r.patient_registry,phase,getattr(r,phase+'_sampleID')) for r in pp.itertuples(index=False) for phase in ['pre','post']]
    materials=pd.read_csv(source/'private_all52_paired_materials.tsv',sep='\t');assert list(materials[['patient_registry','phase','sampleID']].itertuples(index=False,name=None))==material_rows
    ml={key:i for i,(_,_,key) in enumerate(material_rows)};assert len(ml)==52
    scope=root/'immune_RNA_MCP_discovery_scope_20261018_v1';overlay=pd.read_csv(scope/'private_all236_current_case_phase_overlay.tsv',sep='\t');overlay=overlay[overlay.threshold.eq(1)].set_index('sampleID')
    assert np.array_equal(materials.all_material_cells,overlay.loc[materials.sampleID,'annotated_scope_cells'])
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');T=pd.read_csv(old/'private_original_T_unit_scope.tsv',sep='\t')
    cu=cache['original_T_unit_index'];tc=np.array([ml.get(str(r.dataset)+'::'+str(r.library_id),-1) for r in T.itertuples(index=False)])[cu]
    marker_names=cache['genes'].tolist();raw=cache['raw_gene_counts'];fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    def gate(names,k):return np.count_nonzero(raw[:,[marker_names.index(g) for g in names]]>=1,axis=1)>=k
    core=gate(panels['pan_T'],2);reg=gate(panels['regulatory_associated'],2)&gate(['FOXP3'],1);eff=gate(['GZMB','PRF1','CTSW','NKG7'],2)
    oldscope=root/'existing_immune_fine_native_measurement_scope_20261016_v1';units=pd.read_csv(oldscope/'private_all_frozen_T_myeloid_node_units.tsv',sep='\t');mapping=np.load(oldscope/'private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    um=np.array([ml.get(str(r.dataset)+'::'+str(r.library_id),-1) if r.major=='Myeloid' else -1 for r in units.itertuples(index=False)])
    rr=np.flatnonzero(mapping>=0);mm=um[mapping[rr]]
    rows=[]
    for vi,mask in enumerate([core&reg&~eff&(tc>=0),core&eff&~reg&(tc>=0)]):
        rows.extend(zip(cache['global_original_T_row_index'][mask].astype(int).tolist(),tc[mask].astype(int).tolist(),[vi]*int(mask.sum())))
    rows.extend(zip(rr[mm>=0].astype(int).tolist(),mm[mm>=0].astype(int).tolist(),[2]*int((mm>=0).sum())))
    expected=pd.DataFrame(rows,columns=['global_row','material_index','view_index']).sort_values('global_row').reset_index(drop=True);assert expected.global_row.is_unique
    expected['digest']=[hashlib.sha256(('2026101802|'+VIEWS[vi]+'|'+str(r)).encode()).hexdigest() for r,vi in expected[['global_row','view_index']].itertuples(index=False,name=None)]
    primary=expected.sort_values('digest').groupby(['material_index','view_index'],sort=False).head(100).sort_values('global_row').reset_index(drop=True)
    selected=pd.read_csv(source/'private_selected_cells_and_all_scores.tsv.gz',sep='\t');assert list(primary[['global_row','material_index','view_index']].itertuples(index=False,name=None))==list(selected[['global_row','material_index','view_index']].itertuples(index=False,name=None))
    cover=pd.read_csv(source/'private_all156_capture_scope.tsv',sep='\t').set_index(['material_index','view']);assert len(cover)==156
    for (m,vi),g in expected.groupby(['material_index','view_index']):
        c=cover.loc[(m,VIEWS[vi])];assert c.full_view_cells==len(g)>=20 and c.selected_cells==min(100,len(g)) and c.all_material_cells==materials.iloc[m].all_material_cells
    physical=np.load(source/'private_union_native_columns.npy');native_to_union={int(c):i for i,c in enumerate(physical)}
    gm=pd.read_csv(scope/'private_marker_excluded_native_gene_map.tsv',sep='\t').set_index('gene')
    wpath=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wpath)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7'
    weights=pd.read_csv(wpath,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    contract=next(r for r in json.loads((old/'private_source_gene_measurement_contract.json').read_text()) if r['dataset']=='GSE286827')
    markers=np.array(contract['physical_columns'],int);assert np.array_equal(physical,np.array(sorted(set(gm.loc[weights.gene,'GSE286827'].to_numpy(int))|set(markers)),int))
    params=pd.read_csv(source/'private_all12_original_frozen_models.tsv',sep='\t').set_index(['view','program','variant'])
    inputs=root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs';models={};training={};maxerr=0.
    for vi,view in enumerate(VIEWS):
        md=pd.read_csv(inputs/(view+'_discovery_metadata.tsv'),sep='\t');training[vi]=md
        D=np.column_stack([np.ones(len(md)),md[['cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float)]);Q,R=np.linalg.qr(D,mode='reduced')
        names=Path(str(inputs/(view+'_discovery'))+'_genes.tsv').read_text().splitlines()
        with gzip.open(str(inputs/(view+'_discovery'))+'_tpm.mtx.gz','rb') as h:tr=mmread(h).tocsr()
        frozen=np.load(inputs/(view+'_frozen_projection.npz'));models[vi]=dict(mu=frozen['quality_mean'],sd=frozen['quality_sd'],pieces={})
        for pi,program in enumerate(['MCP1','MCP2']):
            w=weights[weights.view.eq(view)&weights.program.eq(program)];b=w.weight.to_numpy()/w.discovery_sd.to_numpy();center=float(np.sum(b*w.discovery_mean.to_numpy()))
            tidx=pd.Index(names).get_indexer(w.gene);native=gm.loc[w.gene,'GSE286827'].to_numpy(int);ui=np.array([native_to_union[int(c)] for c in native]);sg=w.gene.str.upper().isin(STRESS).to_numpy()
            fulltrain=np.asarray(tr[tidx].T.multiply(b).sum(axis=1)).ravel()-center
            stresstrain=np.asarray(tr[tidx[sg]].T.multiply(b[sg]).sum(axis=1)).ravel()-float(np.sum(b[sg]*w.discovery_mean.to_numpy()[sg])) if sg.any() else np.zeros(len(md))
            for si,variant in enumerate(['full','fixed_stress_complement']):
                y=fulltrain if si==0 else fulltrain-stresstrain;coef=np.linalg.solve(R,Q.T@y);pr=params.loc[(view,program,variant)]
                close(pr[['intercept','cellQ','geneQ','sourceCC','fraction','logN']].to_numpy(float),coef)
                ix=np.ones(len(w),bool) if si==0 else ~sg
                close(pr.center,float(np.sum(b[ix]*w.discovery_mean.to_numpy()[ix])));assert pr.active_genes==int(ix.sum())
                models[vi]['pieces'][2*pi+si]=dict(u=ui[ix],b=b[ix],center=float(np.sum(b[ix]*w.discovery_mean.to_numpy()[ix])),coef=coef)
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r');native=z['layers/counts'];nc,nf=map(int,native.attrs['shape']);assert len(mapping)==nc
    erows=expected.global_row.to_numpy(int);n=len(erows);RNA=np.zeros(n,np.int64);det=np.zeros(n,np.int32);qs=np.empty((n,4));values=np.empty((n,4));selected_raw=load_npz(source/'private_selected_raw_union_counts.npz');assert selected_raw.shape==(len(selected),len(physical))
    fullcounts=np.zeros((52,3),int)
    for (m,vi),g in expected.groupby(['material_index','view_index']):fullcounts[m,vi]=len(g)
    selected_ptr=0;slot_count=0
    for lo in range(0,n,1024):
        hi=min(n,lo+1024);batch=erows[lo:hi];first=int(batch[0]);last=int(batch[-1])+1
        ptr=np.asarray(native['indptr'][first:last+1],np.int64);x=csr_matrix((np.asarray(native['data'][int(ptr[0]):int(ptr[-1])]),np.asarray(native['indices'][int(ptr[0]):int(ptr[-1])]),ptr-ptr[0]),shape=(last-first,nf))[batch-first].copy()
        x.sum_duplicates();x.eliminate_zeros();x.sort_indices();assert np.isfinite(x.data).all() and (x.data>0).all() and np.equal(x.data,np.rint(x.data)).all()
        RNA[lo:hi]=np.asarray(x.sum(axis=1)).ravel().astype(np.int64);det[lo:hi]=np.diff(x.indptr);assert (RNA[lo:hi]>0).all()
        dense=x[:,physical].toarray();anchor=(selected.global_row.to_numpy()>=batch[0])&(selected.global_row.to_numpy()<=batch[-1]);ind=np.searchsorted(batch,selected.loc[anchor,'global_row'].to_numpy(int));ns=len(ind)
        if ns:assert np.array_equal(selected_raw[selected_ptr:selected_ptr+ns].toarray(),dense[ind]);selected_ptr+=ns;slot_count+=ns*len(physical)
        logged=np.log1p(dense*10000/RNA[lo:hi,None])
        for vi,model in models.items():
            idx=np.flatnonzero(expected.iloc[lo:hi].view_index.to_numpy()==vi)
            if not len(idx):continue
            ids=lo+idx;mi=expected.iloc[ids].material_index.to_numpy(int);cnt=fullcounts[mi,vi]
            rawq=np.stack([np.log1p(RNA[ids]),np.log1p(det[ids]),cnt/materials.iloc[mi].all_material_cells.to_numpy(int),np.log1p(cnt)],axis=1);q=(rawq-model['mu'])/model['sd'];qs[ids]=q
            for k,p in model['pieces'].items():values[ids,k]=np.sum(logged[idx][:,p['u']]*p['b'],axis=1)-p['center']-np.column_stack([np.ones(len(ids)),q])@p['coef'][[0,1,2,4,5]]
            if vi!=2:
                tix=np.searchsorted(cache['global_original_T_row_index'],erows[ids]);assert np.array_equal(cache['global_original_T_row_index'][tix],erows[ids])
                assert np.array_equal(dense[idx][:,[native_to_union[int(c)] for c in markers]],raw[tix])
                assert np.array_equal(RNA[ids],cache['full_native_total_counts'][tix]) and np.array_equal(det[ids],cache['full_native_detected_genes'][tix])
        if lo//1024%16==0 or hi==n:(out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_ALL_CAPTURED_NATIVE_ROW_REPLAY',captured_cells_replayed=hi,captured_cells=n,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert selected_ptr==len(selected) and np.isfinite(values).all()
    si=pd.Index(erows).get_indexer(selected.global_row);assert (si>=0).all()
    assert np.array_equal(selected.fullRNA,RNA[si]) and np.array_equal(selected.full_detected_genes,det[si]);close(selected[FIELDS].to_numpy(),qs[si])
    for k,f in enumerate(['MCP1_full','MCP1_fixed_stress_complement','MCP2_full','MCP2_fixed_stress_complement']):close(selected[f],values[si,k]);maxerr=max(maxerr,float(np.max(np.abs(selected[f].to_numpy()-values[si,k]))))
    measures=pd.read_csv(source/'private_all624_material_programme_measurements.tsv',sep='\t').set_index(['material_index','view','program','variant']);assert len(measures)==624
    domain=pd.read_csv(source/'private_all12_training_domain_descriptions.tsv',sep='\t').set_index(['view','covariate']);assert len(domain)==12
    for vi,view in enumerate(VIEWS):
        for j,f in enumerate(FIELDS):
            vv=qs[expected.view_index.eq(vi),j];a=float(training[vi][f].min());b=float(training[vi][f].max());r=domain.loc[(view,f)]
            close([r.training_min,r.training_max,r.new_min,r.new_max],[a,b,vv.min(),vv.max()]);assert r.captured_cells==len(vv) and r.outside_training_cells==int(((vv<a)|(vv>b)).sum())
        for m in range(52):
            ix=np.flatnonzero(expected.view_index.eq(vi)&expected.material_index.eq(m));sub=si[(selected.view_index.to_numpy()==vi)&(selected.material_index.to_numpy()==m)]
            for pi,program in enumerate(['MCP1','MCP2']):
                for st,variant in enumerate(['full','fixed_stress_complement']):
                    r=measures.loc[(m,view,program,variant)];k=pi*2+st;close([r.score,r.all_captured_mean_sensitivity],[values[sub,k].mean(),values[ix,k].mean()]);close(r[FIELDS].to_numpy(float),qs[sub].mean(axis=0))
                    assert r.sampleID==materials.iloc[m].sampleID and r.patient_registry==materials.iloc[m].patient_registry and r.phase==materials.iloc[m].phase and r.selected_cells==len(sub) and r.full_view_cells==len(ix) and r.all_material_cells==materials.iloc[m].all_material_cells
    official=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_case_mean_scores.tsv',sep='\t');official=official[official.fit_id.eq('K2_seed20261018_source_adjusted')&official.partition.eq('discovery')]
    partition=pd.read_csv(scope/'private_metadata_only_case_source_partition.tsv',sep='\t').set_index('fit_sample').sampleID;official['sampleID']=official['sample'].map(partition)
    scales=pd.read_csv(source/'private_all6_discovery_case_scales.tsv',sep='\t').set_index(['view','program'])
    for (v,p),g in official.groupby(['view','program']):close(scales.loc[(v,p),'discovery_case_sd'],np.std(g.score.to_numpy(),ddof=0))
    anchor=pd.read_csv(source/'private_all_original_discovery_anchor_checks.tsv',sep='\t')
    baseline=official.set_index(['sampleID','view','program']).score;actual=measures.reset_index();actual=actual[actual.variant.eq('full')]
    checks=0
    for r in actual.itertuples(index=False):
        key=(r.sampleID,r.view,r.program)
        if key in baseline.index:close(r.score,baseline.loc[key]);checks+=1
    assert checks==len(anchor)==term['original_anchor_measurements'];close(anchor.score_new,anchor.score_accepted)
    t=dict(status='ACCEPTED_ALL52_PAIRED_NATIVE_CAPTURE_COUNTS_FROZEN_FULL_COMPLEMENT_SCORES_AND_OFFICIAL_ANCHORS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),captured_cells_replayed=n,selected_cells_replayed=len(selected),selected_physical_count_slots_replayed=slot_count,all_material_measurements_replayed=624,all_original_anchor_scores_replayed=checks,maximum_selected_score_error=maxerr,no_paired_or_clinical_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
